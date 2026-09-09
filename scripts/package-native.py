#!/usr/bin/env python3
"""Build and reproducibly package the versioned Axiolid C ABI."""

from __future__ import annotations

import argparse
import datetime as dt
import gzip
import hashlib
import json
import os
import shutil
import subprocess
import sys
import tarfile
import tempfile
import zipfile
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
VERSION = "0.1.1"
SUPPORTED_TARGETS = {
    "x86_64-unknown-linux-gnu",
    "aarch64-unknown-linux-gnu",
    "x86_64-apple-darwin",
    "x86_64-pc-windows-msvc",
    "aarch64-pc-windows-msvc",
}


def run(*args: str, env: dict[str, str] | None = None) -> str:
    return subprocess.run(
        args,
        cwd=ROOT,
        env=env,
        check=True,
        text=True,
        stdout=subprocess.PIPE,
        stderr=subprocess.STDOUT,
    ).stdout.strip()


def host_triple() -> str:
    for line in run("rustc", "+1.88.0", "-vV").splitlines():
        if line.startswith("host: "):
            return line.removeprefix("host: ")
    raise RuntimeError("rustc did not report a host triple")


def layout(target: str) -> dict[str, str]:
    if target not in SUPPORTED_TARGETS:
        raise ValueError(f"unsupported native target: {target}")
    if "windows-msvc" in target:
        processor = (
            "^(aarch64|arm64|ARM64)$"
            if target.startswith("aarch64-")
            else "^(x86_64|AMD64|amd64)$"
        )
        return {
            "shared": "axiolid_capi.dll",
            "shared_location": "bin/axiolid_capi.dll",
            "implib": "axiolid_capi.dll.lib",
            "static": "axiolid_capi.lib",
            "format": "zip",
            # `rustc --print=native-static-libs` for `axiolid-capi` on
            # x86_64-pc-windows-msvc (verified via cross-compile, both
            # dev and release profiles, 2026-09-03): std pulls these in
            # for env/threading/RNG/error-reporting. msvcrt is handled by
            # CMake's own runtime-library selection, not listed here.
            "static_links": "kernel32.lib;ntdll.lib;userenv.lib;ws2_32.lib;dbghelp.lib",
            "system": "Windows",
            "processor_regex": processor,
        }
    if "apple-darwin" in target:
        return {
            "shared": "libaxiolid_capi.dylib",
            "shared_location": "lib/libaxiolid_capi.dylib",
            "implib": "",
            "static": "libaxiolid_capi.a",
            "format": "tar.gz",
            "static_links": "-framework CoreFoundation;-framework Security",
            "system": "Darwin",
            "processor_regex": "^(x86_64|AMD64|amd64)$",
        }
    if "linux" in target:
        processor = (
            "^(aarch64|arm64|ARM64)$"
            if target.startswith("aarch64-")
            else "^(x86_64|AMD64|amd64)$"
        )
        return {
            "shared": "libaxiolid_capi.so",
            "shared_location": "lib/libaxiolid_capi.so",
            "implib": "",
            "static": "libaxiolid_capi.a",
            "format": "tar.gz",
            "static_links": "dl;pthread;m",
            "system": "Linux",
            "processor_regex": processor,
        }
    raise ValueError(f"unsupported native target: {target}")


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for chunk in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def git_metadata() -> tuple[str, int, bool]:
    commit = run("git", "rev-parse", "HEAD")
    epoch = int(run("git", "show", "-s", "--format=%ct", "HEAD"))
    dirty = bool(run("git", "status", "--porcelain", "--untracked-files=normal"))
    return commit, epoch, dirty


def build(target: str, profile: str, target_dir: Path) -> Path:
    env = os.environ.copy()
    env["CARGO_TARGET_DIR"] = str(target_dir)
    run("cargo", "+1.88.0", "xtask", "ffi", "check", env=env)
    command = [
        "cargo",
        "+1.88.0",
        "build",
        "--package",
        "axiolid-capi",
        "--target",
        target,
    ]
    if profile == "release":
        command.append("--release")
    run(*command, env=env)
    return target_dir / target / ("release" if profile == "release" else "debug")


def render_targets(values: dict[str, str]) -> str:
    template = (ROOT / "native/cmake/AxiolidTargets.cmake.in").read_text()
    if values["implib"]:
        implib = (
            "set_target_properties(Axiolid::axiolid_shared PROPERTIES\n"
            f'  IMPORTED_IMPLIB "${{PACKAGE_PREFIX_DIR}}/lib/{values["implib"]}"\n)'
        )
    else:
        implib = ""
    replacements = {
        "@SHARED_LOCATION@": values["shared_location"],
        "@SHARED_IMPLIB_PROPERTY@": implib,
        "@STATIC_LOCATION@": f'lib/{values["static"]}',
        "@STATIC_LINK_LIBRARIES@": values["static_links"],
        "@SYSTEM_NAME@": values["system"],
        "@PROCESSOR_REGEX@": values["processor_regex"],
    }
    for marker, value in replacements.items():
        template = template.replace(marker, value)
    if "@" in template:
        raise RuntimeError("unresolved CMake template marker")
    return template


def copy_payload(stage: Path, artifact_dir: Path, values: dict[str, str]) -> None:
    (stage / "include").mkdir(parents=True)
    (stage / "lib/cmake/Axiolid").mkdir(parents=True)
    (stage / "bin").mkdir(parents=True)
    shutil.copy2(
        ROOT / "crates/facade/axiolid-capi/include/axiolid.h",
        stage / "include/axiolid.h",
    )
    shutil.copy2(ROOT / "LICENSE", stage / "LICENSE")
    shared_destination = stage / values["shared_location"]
    shared_destination.parent.mkdir(parents=True, exist_ok=True)
    shutil.copy2(artifact_dir / values["shared"], shared_destination)
    shutil.copy2(artifact_dir / values["static"], stage / "lib" / values["static"])
    if values["implib"]:
        shutil.copy2(artifact_dir / values["implib"], stage / "lib" / values["implib"])
    cmake_dir = stage / "lib/cmake/Axiolid"
    shutil.copy2(ROOT / "native/cmake/AxiolidConfig.cmake", cmake_dir)
    shutil.copy2(ROOT / "native/cmake/AxiolidConfigVersion.cmake", cmake_dir)
    shutil.copy2(ROOT / "native/cmake/AxiolidFetch.cmake", cmake_dir)
    (cmake_dir / "AxiolidTargets.cmake").write_text(render_targets(values))


def write_manifest(
    stage: Path, target: str, profile: str, commit: str, dirty: bool
) -> None:
    files = {}
    for path in sorted(item for item in stage.rglob("*") if item.is_file()):
        relative = path.relative_to(stage).as_posix()
        files[relative] = {"sha256": sha256(path), "size": path.stat().st_size}
    manifest = {
        "schema_version": 1,
        "package_version": VERSION,
        "abi_version": "0.4",
        "target": target,
        "profile": profile,
        "source_commit": commit,
        "source_dirty": dirty,
        "rustc": run("rustc", "+1.88.0", "-vV"),
        "files": files,
    }
    (stage / "manifest.json").write_text(
        json.dumps(manifest, indent=2, sort_keys=True) + "\n"
    )


def tar_archive(stage: Path, archive: Path, epoch: int) -> None:
    with archive.open("wb") as raw:
        with gzip.GzipFile(
            filename="", mode="wb", fileobj=raw, mtime=epoch
        ) as compressed:
            with tarfile.open(
                fileobj=compressed, mode="w", format=tarfile.USTAR_FORMAT
            ) as tar:
                for path in sorted(
                    [stage, *stage.rglob("*")], key=lambda item: item.as_posix()
                ):
                    relative = Path(stage.name) / path.relative_to(stage)
                    info = tar.gettarinfo(str(path), arcname=relative.as_posix())
                    info.uid = info.gid = 0
                    info.uname = info.gname = ""
                    info.mtime = epoch
                    info.mode = 0o755 if path.is_dir() else 0o644
                    if path.is_file():
                        with path.open("rb") as stream:
                            tar.addfile(info, stream)
                    else:
                        tar.addfile(info)


def zip_archive(stage: Path, archive: Path, epoch: int) -> None:
    stamp = dt.datetime.fromtimestamp(max(epoch, 315532800), tz=dt.UTC)
    date_time = (
        stamp.year,
        stamp.month,
        stamp.day,
        stamp.hour,
        stamp.minute,
        stamp.second,
    )
    with zipfile.ZipFile(
        archive, "w", compression=zipfile.ZIP_DEFLATED, compresslevel=9
    ) as bundle:
        for path in sorted(item for item in stage.rglob("*") if item.is_file()):
            name = (Path(stage.name) / path.relative_to(stage)).as_posix()
            info = zipfile.ZipInfo(name, date_time=date_time)
            info.compress_type = zipfile.ZIP_DEFLATED
            info.external_attr = 0o100644 << 16
            bundle.writestr(info, path.read_bytes(), compresslevel=9)


def package(args: argparse.Namespace) -> Path:
    target = args.target or host_triple()
    values = layout(target)
    commit, git_epoch, dirty = git_metadata()
    if dirty and not args.allow_dirty:
        raise RuntimeError(
            "refusing to package a dirty tree; commit it or pass --allow-dirty for local testing"
        )
    epoch = int(os.environ.get("SOURCE_DATE_EPOCH", git_epoch))
    target_dir = Path(os.environ.get("CARGO_TARGET_DIR", ROOT / "target")).resolve()
    artifact_dir = build(target, args.profile, target_dir)
    output_dir = args.output_dir.resolve()
    output_dir.mkdir(parents=True, exist_ok=True)
    stem = f"axiolid-native-v{VERSION}-{target}"
    extension = values["format"]
    archive = output_dir / f"{stem}.{extension}"
    with tempfile.TemporaryDirectory(prefix="axiolid-native-") as temporary:
        stage = Path(temporary) / stem
        stage.mkdir()
        copy_payload(stage, artifact_dir, values)
        write_manifest(stage, target, args.profile, commit, dirty)
        if extension == "zip":
            zip_archive(stage, archive, epoch)
        else:
            tar_archive(stage, archive, epoch)
    checksum = sha256(archive)
    archive.with_suffix(archive.suffix + ".sha256").write_text(
        f"{checksum}  {archive.name}\n"
    )
    return archive


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--target", help="Rust target triple; defaults to rustc host")
    parser.add_argument("--profile", choices=("debug", "release"), default="release")
    parser.add_argument("--output-dir", type=Path, default=ROOT / "dist/native")
    parser.add_argument(
        "--allow-dirty",
        action="store_true",
        help="local testing only; recorded in manifest",
    )
    return parser.parse_args()


def main() -> int:
    try:
        archive = package(parse_args())
    except (OSError, ValueError, RuntimeError, subprocess.CalledProcessError) as error:
        print(f"native package failed: {error}", file=sys.stderr)
        return 1
    print(archive)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
