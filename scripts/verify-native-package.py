#!/usr/bin/env python3
"""Verify an Axiolid native archive before installation or publication."""

from __future__ import annotations

import argparse
import hashlib
import json
import re
import stat
import sys
import tarfile
import zipfile
from pathlib import Path, PurePosixPath

MAX_MEMBER_BYTES = 128 * 1024 * 1024
MAX_ARCHIVE_BYTES = 256 * 1024 * 1024
SUPPORTED_TARGETS = {
    "x86_64-unknown-linux-gnu",
    "aarch64-unknown-linux-gnu",
    "x86_64-apple-darwin",
    "x86_64-pc-windows-msvc",
    "aarch64-pc-windows-msvc",
}


def digest(data: bytes) -> str:
    return hashlib.sha256(data).hexdigest()


def digest_file(path: Path) -> str:
    if path.stat().st_size > MAX_ARCHIVE_BYTES:
        raise ValueError("archive file exceeds verification size budget")
    value = hashlib.sha256()
    with path.open("rb") as stream:
        for chunk in iter(lambda: stream.read(1024 * 1024), b""):
            value.update(chunk)
    return value.hexdigest()


def safe_name(name: str) -> PurePosixPath:
    path = PurePosixPath(name)
    if (
        path.is_absolute()
        or ".." in path.parts
        or not path.parts
        or "\\" in name
        or ":" in path.parts[0]
    ):
        raise ValueError(f"unsafe archive path: {name!r}")
    return path


def read_archive(path: Path) -> dict[str, bytes]:
    files: dict[str, bytes] = {}
    total_size = 0
    if path.name.endswith(".zip"):
        with zipfile.ZipFile(path) as bundle:
            for info in bundle.infolist():
                name = safe_name(info.filename).as_posix()
                if info.is_dir():
                    continue
                if stat.S_ISLNK(info.external_attr >> 16):
                    raise ValueError(f"symbolic link archive member: {name}")
                if (
                    info.file_size > MAX_MEMBER_BYTES
                    or total_size + info.file_size > MAX_ARCHIVE_BYTES
                ):
                    raise ValueError("archive exceeds verification size budget")
                if name in files:
                    raise ValueError(f"duplicate archive member: {name}")
                files[name] = bundle.read(info)
                total_size += info.file_size
    elif path.name.endswith(".tar.gz"):
        with tarfile.open(path, "r:gz") as bundle:
            for member in bundle.getmembers():
                name = safe_name(member.name).as_posix()
                if member.isdir():
                    continue
                if not member.isfile():
                    raise ValueError(f"non-regular archive member: {name}")
                if (
                    member.size > MAX_MEMBER_BYTES
                    or total_size + member.size > MAX_ARCHIVE_BYTES
                ):
                    raise ValueError("archive exceeds verification size budget")
                if name in files:
                    raise ValueError(f"duplicate archive member: {name}")
                stream = bundle.extractfile(member)
                if stream is None:
                    raise ValueError(f"unreadable archive member: {name}")
                files[name] = stream.read()
                total_size += member.size
    else:
        raise ValueError("archive must end in .tar.gz or .zip")
    return files


def has_pe_image_header(data: bytes) -> bool:
    if not data.startswith(b"MZ") or len(data) < 64:
        return False
    offset = int.from_bytes(data[60:64], "little")
    return offset + 6 <= len(data) and data[offset : offset + 4] == bytes(
        (0x50, 0x45, 0, 0)
    )


def machine(data: bytes) -> str | None:
    if data.startswith(b"\x7fELF") and len(data) >= 20:
        endian = "little" if data[5] == 1 else "big"
        return {62: "x86_64", 183: "aarch64"}.get(int.from_bytes(data[18:20], endian))
    if data[:4] in (b"\xcf\xfa\xed\xfe", b"\xfe\xed\xfa\xcf") and len(data) >= 8:
        endian = "little" if data[:4] == b"\xcf\xfa\xed\xfe" else "big"
        return {0x01000007: "x86_64", 0x0100000C: "aarch64"}.get(
            int.from_bytes(data[4:8], endian)
        )
    if has_pe_image_header(data):
        offset = int.from_bytes(data[60:64], "little")
        return {0x8664: "x86_64", 0xAA64: "aarch64"}.get(
            int.from_bytes(data[offset + 4 : offset + 6], "little")
        )
    if len(data) >= 8 and data[:4] == bytes((0, 0, 0xFF, 0xFF)):
        return {0x8664: "x86_64", 0xAA64: "aarch64"}.get(
            int.from_bytes(data[6:8], "little")
        )
    if has_coff_object_header(data):
        return {0x8664: "x86_64", 0xAA64: "aarch64"}.get(
            int.from_bytes(data[:2], "little")
        )
    return None


def has_strong_object_magic(data: bytes) -> bool:
    return (
        data.startswith(b"\x7fELF")
        or data[:4]
        in {
            b"\xcf\xfa\xed\xfe",
            b"\xfe\xed\xfa\xcf",
            b"\xce\xfa\xed\xfe",
            b"\xfe\xed\xfa\xce",
        }
        or (len(data) >= 8 and data[:4] == bytes((0, 0, 0xFF, 0xFF)))
    )


def has_coff_object_header(data: bytes) -> bool:
    if len(data) < 20 or data[:4] == bytes((0, 0, 0xFF, 0xFF)):
        return False
    section_count = int.from_bytes(data[2:4], "little")
    symbol_table = int.from_bytes(data[8:12], "little")
    symbol_count = int.from_bytes(data[12:16], "little")
    optional_header = int.from_bytes(data[16:18], "little")
    if not 0 < section_count <= 0xFEFF or optional_header != 0:
        return False
    if len(data) < 20 + 40 * section_count:
        return False
    return symbol_table == 0 or symbol_table + 18 * symbol_count <= len(data)


def archive_machines(data: bytes) -> set[str]:
    if not data.startswith(b"!<arch>\n"):
        raise ValueError("static library is not an ar archive")
    found: set[str] = set()
    long_names = b""
    position = 8
    while position + 60 <= len(data):
        header = data[position : position + 60]
        if header[58:60] != b"`\n":
            raise ValueError("malformed ar member header")
        try:
            size = int(header[48:58].decode("ascii").strip())
        except ValueError as error:
            raise ValueError("malformed ar member size") from error
        payload_start = position + 60
        payload_end = payload_start + size
        if payload_end > len(data):
            raise ValueError("ar member exceeds archive length")
        payload = data[payload_start:payload_end]
        raw_name = header[:16].decode("ascii", errors="replace").strip()
        name = raw_name.rstrip("/")
        if raw_name == "//":
            long_names = payload
        elif raw_name.startswith("#1/"):
            name_length = int(raw_name[3:])
            if name_length > len(payload):
                raise ValueError("BSD ar filename exceeds member length")
            name = (
                payload[:name_length].decode("utf-8", errors="replace").rstrip(chr(0))
            )
            payload = payload[name_length:]
        elif raw_name.startswith("/") and raw_name[1:].isdigit():
            offset = int(raw_name[1:])
            if offset >= len(long_names):
                raise ValueError("GNU ar filename offset exceeds string table")
            name = (
                long_names[offset:]
                .split(b"/\n", 1)[0]
                .decode("utf-8", errors="replace")
            )
        detected = machine(payload)
        special_member = raw_name in {"/", "//", "/SYM64/"}
        object_member = not special_member and (
            has_strong_object_magic(payload) or has_coff_object_header(payload)
        )
        if object_member and detected is None:
            raise ValueError(
                f"archive member has unsupported object architecture: {name}"
            )
        if object_member:
            found.add(detected)
        position += 60 + size + (size % 2)
    if position != len(data):
        raise ValueError("trailing or truncated ar data")
    if not found:
        raise ValueError("static archive contains no recognized native object")
    return found


def expected_architecture(target: str) -> str:
    if target.startswith("x86_64-"):
        return "x86_64"
    if target.startswith("aarch64-"):
        return "aarch64"
    raise ValueError(f"unsupported architecture in target {target}")


def verify_checksum(archive: Path) -> None:
    sidecar = archive.with_suffix(archive.suffix + ".sha256")
    if not sidecar.is_file():
        raise ValueError(f"missing checksum sidecar: {sidecar.name}")
    fields = sidecar.read_text().strip().split()
    actual = digest_file(archive)
    if len(fields) != 2 or fields[0] != actual or fields[1] != archive.name:
        raise ValueError("checksum sidecar does not match archive")


def verify_files(files: dict[str, bytes], allow_dirty: bool) -> tuple[str, dict]:
    roots = {PurePosixPath(name).parts[0] for name in files}
    if len(roots) != 1:
        raise ValueError("archive must contain exactly one top-level directory")
    root = roots.pop()
    manifest_name = f"{root}/manifest.json"
    try:
        manifest = json.loads(files[manifest_name])
    except KeyError as error:
        raise ValueError("archive has no manifest.json") from error
    if not isinstance(manifest, dict):
        raise ValueError("manifest root must be an object")  # noqa: TRY004
    expected_keys = {
        "schema_version",
        "package_version",
        "abi_version",
        "source_commit",
        "source_dirty",
        "target",
        "profile",
        "rustc",
        "files",
    }
    if set(manifest) != expected_keys:
        raise ValueError("manifest keys differ from schema v1")
    if manifest.get("schema_version") != 1 or manifest.get("abi_version") != "0.4":
        raise ValueError("unsupported manifest or ABI version")
    if manifest.get("package_version") != "0.1.1":
        raise ValueError("unsupported native package version")
    target = manifest.get("target")
    if not isinstance(target, str) or target not in SUPPORTED_TARGETS:
        raise ValueError("unsupported target in manifest")
    profile = manifest.get("profile")
    if not isinstance(profile, str) or profile not in {"debug", "release"}:
        raise ValueError("unsupported build profile in manifest")
    source_commit = manifest.get("source_commit")
    if not isinstance(source_commit, str) or not re.fullmatch(
        r"[0-9a-f]{40}", source_commit
    ):
        raise ValueError("manifest source commit is not immutable")
    rustc = manifest.get("rustc")
    if (
        not isinstance(rustc, str)
        or not rustc.splitlines()
        or rustc.splitlines()[0] != "rustc 1.88.0 (6b00bc388 2025-06-23)"
    ):
        raise ValueError("manifest was not built with the supported Rust toolchain")
    if not isinstance(manifest.get("source_dirty"), bool):
        raise ValueError("manifest source_dirty must be a boolean")  # noqa: TRY004
    if manifest["source_dirty"] and not allow_dirty:
        raise ValueError("release archive records a dirty source tree")
    expected_root = (
        f"axiolid-native-v{manifest['package_version']}-{manifest['target']}"
    )
    if root != expected_root:
        raise ValueError(
            f"archive root {root!r} does not match manifest target/version"
        )
    required = {
        "LICENSE",
        "include/axiolid.h",
        "lib/cmake/Axiolid/AxiolidConfig.cmake",
        "lib/cmake/Axiolid/AxiolidConfigVersion.cmake",
        "lib/cmake/Axiolid/AxiolidFetch.cmake",
        "lib/cmake/Axiolid/AxiolidTargets.cmake",
    }
    relative_files = {
        str(PurePosixPath(name).relative_to(root))
        for name in files
        if name != manifest_name
    }
    missing = required - relative_files
    if missing:
        raise ValueError(f"archive is missing required files: {sorted(missing)}")
    recorded = manifest.get("files")
    if not isinstance(recorded, dict):
        raise ValueError("manifest files must be an object")  # noqa: TRY004
    if set(recorded) != relative_files:
        raise ValueError("manifest file inventory differs from archive")
    for relative, metadata in recorded.items():
        payload = files[f"{root}/{relative}"]
        if metadata != {"sha256": digest(payload), "size": len(payload)}:
            raise ValueError(f"manifest digest/size mismatch: {relative}")
    return root, manifest


def verify_binaries(files: dict[str, bytes], root: str, manifest: dict) -> None:
    target = manifest.get("target", "")
    expected = expected_architecture(target)
    if "windows-msvc" in target:
        shared, static, implib = (
            "bin/axiolid_capi.dll",
            "lib/axiolid_capi.lib",
            "lib/axiolid_capi.dll.lib",
        )
    elif "apple-darwin" in target:
        shared, static, implib = (
            "lib/libaxiolid_capi.dylib",
            "lib/libaxiolid_capi.a",
            None,
        )
    elif "linux" in target:
        shared, static, implib = "lib/libaxiolid_capi.so", "lib/libaxiolid_capi.a", None
    else:
        raise ValueError(f"unsupported target in manifest: {target}")
    required = [shared, static, *([implib] if implib else [])]
    for relative in required:
        if f"{root}/{relative}" not in files:
            raise ValueError(f"missing native artifact: {relative}")
    shared_payload = files[f"{root}/{shared}"]
    if "windows-msvc" in target:
        shared_magic_valid = has_pe_image_header(shared_payload)
    elif "apple-darwin" in target:
        shared_magic_valid = shared_payload[:4] in {
            b"\xcf\xfa\xed\xfe",
            b"\xfe\xed\xfa\xcf",
            b"\xce\xfa\xed\xfe",
            b"\xfe\xed\xfa\xce",
        }
    else:
        shared_magic_valid = shared_payload.startswith(b"\x7fELF")
    if not shared_magic_valid:
        raise ValueError(
            "shared library does not have the target platform's image format"
        )
    actual_shared = machine(shared_payload)
    if actual_shared != expected:
        raise ValueError(
            f"shared library machine is {actual_shared}, expected {expected}"
        )
    static_machines = archive_machines(files[f"{root}/{static}"])
    if static_machines != {expected}:
        raise ValueError(
            f"static library machines are {sorted(static_machines)}, expected {expected}"
        )
    if implib:
        import_machines = archive_machines(files[f"{root}/{implib}"])
        if import_machines != {expected}:
            raise ValueError(
                f"import library machines are {sorted(import_machines)}, expected {expected}"
            )
    header = files[f"{root}/include/axiolid.h"]
    if (
        b"axiolid_v0_4_version" not in header
        or b"Generated by `cargo xtask ffi header`" not in header
    ):
        raise ValueError("public header lacks the v0.4 API or generation marker")


def extract(files: dict[str, bytes], destination: Path) -> Path:
    if destination.exists():
        raise ValueError(f"extraction destination already exists: {destination}")
    destination.mkdir(parents=True)
    for name, payload in files.items():
        path = destination.joinpath(*safe_name(name).parts)
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_bytes(payload)
    root = next(destination.iterdir())
    return root


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("archive", type=Path)
    parser.add_argument("--extract-to", type=Path)
    parser.add_argument(
        "--allow-dirty", action="store_true", help="accept local-test manifests"
    )
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    try:
        verify_checksum(args.archive)
        files = read_archive(args.archive)
        root, manifest = verify_files(files, args.allow_dirty)
        extension = ".zip" if "windows-msvc" in manifest["target"] else ".tar.gz"
        if args.archive.name != f"{root}{extension}":
            raise ValueError("archive filename does not match manifest target/version")
        verify_binaries(files, root, manifest)
        extracted = (
            extract(files, args.extract_to.resolve()) if args.extract_to else None
        )
    except (
        OSError,
        ValueError,
        json.JSONDecodeError,
        tarfile.TarError,
        zipfile.BadZipFile,
    ) as error:
        print(f"native package verification failed: {error}", file=sys.stderr)
        return 1
    print(
        f"verified {args.archive.name}: target={manifest['target']} commit={manifest['source_commit']}"
    )
    if extracted:
        print(extracted)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
