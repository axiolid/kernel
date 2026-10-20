#!/usr/bin/env python3
"""Bootstrap-preflight every publishable workspace crate without uploading.

First releases cannot package dependent crates unpatched before their internal
names exist on crates.io. This preflight therefore compiles normalized,
lock-free source archives with command-scoped local patches. It rejects any
Cargo.lock (and thus patch resolution metadata) in those archives, and proves
byte identity against unpatched archives for every internal-dependency root.

It then checks the crates about to be released -- every publishable crate
whose version is not on crates.io yet -- against their minimal internal
dependency versions (kernel#337). Each one is built in its own scratch
workspace from its packaged archive, against crates.io plus the other
to-be-released archives through `[patch.crates-io]`, with every `axiolid-*`
crate in its graph lowered to the oldest version the published requirements
admit. A consumer whose lockfile holds those versions is a valid consumer, so
a crate that does not build there was published with a floor that is too low;
the failure names that crate.

This is deliberately not the upload-byte gate. ``publish-workspace.py`` later
packages and compiles each exact unpatched, lockful archive after its children
are visible, reproduces it through ``cargo publish --dry-run``, and binds the
uploaded bytes to the crates.io checksum.
"""

from __future__ import annotations

import argparse
from concurrent.futures import ThreadPoolExecutor
import hashlib
import json
import os
from pathlib import Path
import re
import subprocess
import tarfile
import tempfile
import time
import tomllib
import urllib.error
import urllib.request

ROOT = Path(__file__).resolve().parents[1]
CARGO = os.environ.get("CARGO", "cargo")
INTERNAL_PREFIX = "axiolid"
# Pinned so the resolution is reproducible: `-Z minimal-versions` is a
# nightly-only resolver mode. Only the lockfile is made with it; the build
# uses the workspace toolchain (rust-toolchain.toml).
MINIMAL_VERSIONS_TOOLCHAIN = "nightly-2026-09-25"
SPARSE_INDEX = "https://index.crates.io"
USER_AGENT = "axiolid-release-workflow/1.0"
# The facade compiles nothing without a capability feature by design
# (kernel#9); `standard` is the migration target its diagnostic recommends.
DEFAULT_FEATURES = {"axiolid": ["standard"]}
PROBE_PACKAGE = "axiolid-minimal-versions-probe"


def cargo_metadata() -> dict:
    result = subprocess.run(
        [CARGO, "metadata", "--format-version", "1", "--no-deps", "--locked"],
        cwd=ROOT,
        check=True,
        text=True,
        stdout=subprocess.PIPE,
    )
    return json.loads(result.stdout)


def archive_checksum(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for chunk in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def assert_no_packaged_lockfile(archive: Path) -> None:
    with tarfile.open(archive, "r:gz") as package:
        lockfiles = [name for name in package.getnames() if name.endswith("/Cargo.lock")]
    if lockfiles:
        raise SystemExit(
            f"bootstrap archive {archive.name} contains Cargo.lock/patch metadata: {lockfiles}"
        )


def assert_archives_identical(patched: Path, unpatched: Path) -> None:
    patched_checksum = archive_checksum(patched)
    unpatched_checksum = archive_checksum(unpatched)
    if patched_checksum != unpatched_checksum or patched.read_bytes() != unpatched.read_bytes():
        raise SystemExit(
            f"patched/unpatched normalized archive mismatch for {patched.name}: "
            f"patched={patched_checksum} unpatched={unpatched_checksum}"
        )


def workspace_packages(data: dict) -> tuple[list[dict], list[dict]]:
    workspace_ids = set(data["workspace_members"])
    members = sorted(
        (package for package in data["packages"] if package["id"] in workspace_ids),
        key=lambda package: package["name"],
    )
    return (
        [package for package in members if package.get("publish") != []],
        [package for package in members if package.get("publish") == []],
    )


def internal_roots(packages: list[dict]) -> list[dict]:
    by_path = {
        str(Path(package["manifest_path"]).parent.resolve()): package["name"]
        for package in packages
    }
    roots = []
    for package in packages:
        internal = {
            by_path[str(Path(dependency["path"]).resolve())]
            for dependency in package["dependencies"]
            if dependency.get("path")
            and str(Path(dependency["path"]).resolve()) in by_path
        }
        if not internal:
            roots.append(package)
    return roots


def write_patch_config(publishable: list[dict], path: Path) -> None:
    rows = ["[patch.crates-io]"]
    for package in publishable:
        package_path = Path(package["manifest_path"]).parent
        rows.append(f'"{package["name"]}" = {{ path = "{package_path}" }}')
    path.write_text("\n".join(rows) + "\n", encoding="utf-8")


def bootstrap_preflight(
    publishable: list[dict],
    excluded: list[dict],
    temporary_path: Path,
    patch_config: Path,
    allow_dirty: bool,
) -> tuple[Path, int]:
    """Package and verify every publishable archive; return the archive
    directory and the number of byte-identical internal roots."""
    patched_target = temporary_path / "patched-target"
    command = [
        CARGO,
        "package",
        "--workspace",
        "--locked",
        "--quiet",
        "--exclude-lockfile",
    ]
    # The facade compiles nothing without a capability feature by design
    # (kernel#9), so verifying its tarball with default features would hit
    # its own `compile_error!`. It is packaged separately with `standard`
    # -- the migration target its diagnostic recommends -- so the tarball
    # is still verified rather than merely skipped.
    FACADE = "axiolid"
    facade_publishable = any(p["name"] == FACADE for p in publishable)
    for package in excluded:
        command.extend(["--exclude", package["name"]])
    if facade_publishable:
        command.extend(["--exclude", FACADE])
    if allow_dirty:
        command.append("--allow-dirty")
    command.extend(["--config", str(patch_config)])
    environment = os.environ.copy()
    environment["CARGO_TARGET_DIR"] = str(patched_target)
    subprocess.run(command, cwd=ROOT, env=environment, check=True)

    if facade_publishable:
        facade_command = [
            CARGO,
            "package",
            "--locked",
            "--quiet",
            "--exclude-lockfile",
            "-p",
            FACADE,
            "--features",
            "standard",
        ]
        if allow_dirty:
            facade_command.append("--allow-dirty")
        facade_command.extend(["--config", str(patch_config)])
        subprocess.run(facade_command, cwd=ROOT, env=environment, check=True)

    expected = {
        f'{package["name"]}-{package["version"]}.crate' for package in publishable
    }
    actual = {archive.name for archive in (patched_target / "package").glob("*.crate")}
    if actual != expected:
        raise SystemExit(
            "package archive mismatch: "
            f"missing={sorted(expected - actual)}, unexpected={sorted(actual - expected)}"
        )
    for archive in sorted((patched_target / "package").glob("*.crate")):
        assert_no_packaged_lockfile(archive)

    roots = internal_roots(publishable)
    for package in roots:
        root_target = temporary_path / f'unpatched-{package["name"]}'
        root_environment = os.environ.copy()
        root_environment["CARGO_TARGET_DIR"] = str(root_target)
        root_command = [
            CARGO,
            "package",
            "-p",
            package["name"],
            "--locked",
            "--quiet",
            "--exclude-lockfile",
        ]
        if allow_dirty:
            root_command.append("--allow-dirty")
        subprocess.run(root_command, cwd=ROOT, env=root_environment, check=True)
        filename = f'{package["name"]}-{package["version"]}.crate'
        assert_archives_identical(
            patched_target / "package" / filename,
            root_target / "package" / filename,
        )
    return patched_target / "package", len(roots)


# --- minimal internal versions (kernel#337) ---------------------------------


def is_internal(name: str) -> bool:
    if name == PROBE_PACKAGE:
        return False
    return name == INTERNAL_PREFIX or name.startswith(f"{INTERNAL_PREFIX}-")


def sparse_index_path(name: str) -> str:
    """The crates.io sparse-index path of a crate name (Cargo's layout)."""
    lowered = name.lower()
    if len(lowered) <= 2:
        return f"{len(lowered)}/{lowered}"
    if len(lowered) == 3:
        return f"3/{lowered[0]}/{lowered}"
    return f"{lowered[:2]}/{lowered[2:4]}/{lowered}"


def parse_index_versions(text: str) -> set[str]:
    """Every version a sparse-index file lists, yanked ones included: a
    yanked version still occupies its number and cannot be released again."""
    return {json.loads(line)["vers"] for line in text.splitlines() if line.strip()}


def published_versions(name: str, attempts: int = 3) -> set[str]:
    """The versions of `name` on crates.io; an unclear answer is an error.

    Guessing "not published" would put every crate into the check when the
    network is down, and guessing "published" would drop a crate from it.
    """
    request = urllib.request.Request(
        f"{SPARSE_INDEX}/{sparse_index_path(name)}", headers={"User-Agent": USER_AGENT}
    )
    for attempt in range(1, attempts + 1):
        try:
            with urllib.request.urlopen(request, timeout=30) as response:
                return parse_index_versions(response.read().decode("utf-8"))
        except urllib.error.HTTPError as error:
            if error.code in (403, 404, 410):
                return set()
            failure = f"HTTP {error.code}"
        except (OSError, TimeoutError) as error:
            failure = str(error)
        if attempt < attempts:
            time.sleep(2 * attempt)
    raise SystemExit(
        f"cannot tell whether {name} is on crates.io ({failure}); the minimal-versions "
        "check needs the registry index to know which crates are to be released"
    )


def to_be_released(publishable: list[dict], lookup=published_versions) -> list[dict]:
    """The publishable crates whose current version is not on crates.io --
    the same rule publish-workspace.py uses to decide what it uploads."""
    with ThreadPoolExecutor(max_workers=8) as pool:
        versions = list(pool.map(lambda package: lookup(package["name"]), publishable))
    return [
        package
        for package, published in zip(publishable, versions)
        if package["version"] not in published
    ]


def select_packages(publishable: list[dict], names: list[str]) -> list[dict]:
    by_name = {package["name"]: package for package in publishable}
    unknown = sorted(set(names) - set(by_name))
    if unknown:
        raise SystemExit(f"--crates names no publishable workspace crate: {unknown}")
    return [by_name[name] for name in sorted(set(names))]


def feature_sets(package: dict) -> list[tuple[str, list[str]]]:
    """Default features (the facade's `standard`), then every feature.

    Both: an optional dependency can raise a shared floor and so hide a
    default build that does not compile, and a feature-only dependency is
    invisible to a default build.
    """
    default = DEFAULT_FEATURES.get(package["name"], [])
    every = sorted(package.get("features", {}))
    sets = [("default", default)]
    if set(every) - set(default) - {"default"}:
        sets.append(("all-features", every))
    return sets


def probe_workspace_files(
    name: str, version: str, features: list[str], patches: dict[str, Path]
) -> dict[str, str]:
    """A one-member scratch workspace that depends on exactly `name@version`."""
    feature_list = ", ".join(json.dumps(feature) for feature in features)
    rows = ["[workspace]", 'resolver = "2"', 'members = ["probe"]', "", "[patch.crates-io]"]
    for patch_name, path in sorted(patches.items()):
        rows.append(f"{json.dumps(patch_name)} = {{ path = {json.dumps(str(path))} }}")
    probe = [
        "[package]",
        f'name = "{PROBE_PACKAGE}"',
        'version = "0.0.0"',
        'edition = "2021"',
        "publish = false",
        "",
        "[dependencies]",
        f'{json.dumps(name)} = {{ version = "={version}", features = [{feature_list}] }}',
    ]
    return {
        "Cargo.toml": "\n".join(rows) + "\n",
        "probe/Cargo.toml": "\n".join(probe) + "\n",
        "probe/src/lib.rs": "",
    }


def internal_update_specs(lock_text: str) -> list[str]:
    """`name@version` for every internal crate in a lockfile, for the
    `cargo update -p` that lowers exactly those to their minimum."""
    lock = tomllib.loads(lock_text)
    return sorted(
        f'{row["name"]}@{row["version"]}'
        for row in lock.get("package", [])
        if is_internal(row["name"])
    )


def locked_internal_versions(lock_text: str) -> list[tuple[str, str, bool]]:
    """(name, version, from_registry) for every internal crate in a lockfile."""
    lock = tomllib.loads(lock_text)
    return sorted(
        (row["name"], row["version"], row.get("source") is not None)
        for row in lock.get("package", [])
        if is_internal(row["name"])
    )


COULD_NOT_COMPILE = re.compile(r"could not compile `([^`]+)`")


def failing_crates(output: str) -> list[str]:
    """The crates cargo reports it could not compile, in first-seen order."""
    seen: list[str] = []
    for match in COULD_NOT_COMPILE.finditer(output):
        if match.group(1) not in seen and match.group(1) != PROBE_PACKAGE:
            seen.append(match.group(1))
    return seen


def minimal_versions_failure(
    name: str, version: str, label: str, output: str, locked: list[tuple[str, str, bool]]
) -> str:
    failing = failing_crates(output) or [name]
    versions = {row[0]: row for row in locked}
    culprits = []
    for crate in failing:
        row = versions.get(crate)
        if row is None:
            culprits.append(crate)
        else:
            origin = "crates.io" if row[2] else "to be released"
            culprits.append(f"{crate} {row[1]} ({origin})")
    floors = ", ".join(f"{row[0]} {row[1]}" for row in locked if row[0] != name)
    tail = "\n".join(output.strip().splitlines()[-40:])
    return (
        f"MINIMAL_VERSIONS=FAIL {name} {version} ({label}): "
        f"{', '.join(culprits)} does not build with its minimal internal dependency "
        "versions, so its axiolid-* requirements admit versions it cannot use. Raise "
        "the floor it names, or that of the crate re-exporting the missing API, and "
        "release that crate.\n"
        f"minimal internal versions: {floors or '-'}\n{tail}"
    )


def minimal_update_command(manifest: Path, specs: list[str]) -> list[str]:
    """Lower exactly `specs` (the internal crates) to their minimal versions.

    `-Z minimal-versions`, not `-Z direct-minimal-versions`: the stale floor
    this exists for is transitive. `axiolid-mesh-compile` reached
    `axiolid-evaluate` only through `axiolid-reference`, and the direct mode
    keeps every indirect dependency at its newest version, which hides it.
    `-p` confines the mode to the internal crates; everything else stays
    where the normal resolution put it.
    """
    command = [
        "cargo",
        f"+{MINIMAL_VERSIONS_TOOLCHAIN}",
        "-Z",
        "minimal-versions",
        "update",
        "--manifest-path",
        str(manifest),
    ]
    for spec in specs:
        command.extend(["-p", spec])
    return command


def run_logged(command: list[str], *, cwd: Path, env: dict[str, str]) -> subprocess.CompletedProcess:
    return subprocess.run(
        command, cwd=cwd, env=env, text=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT
    )


def require_toolchain() -> None:
    result = subprocess.run(
        ["rustup", "run", MINIMAL_VERSIONS_TOOLCHAIN, "cargo", "--version"],
        text=True,
        stdout=subprocess.PIPE,
        stderr=subprocess.STDOUT,
    )
    if result.returncode != 0:
        raise SystemExit(
            f"the minimal-versions check needs the pinned {MINIMAL_VERSIONS_TOOLCHAIN} "
            f"toolchain: rustup toolchain install {MINIMAL_VERSIONS_TOOLCHAIN} --profile minimal"
        )


def unpack_archives(archives: Path, packages: list[dict], destination: Path) -> dict[str, Path]:
    destination.mkdir(parents=True, exist_ok=True)
    unpacked: dict[str, Path] = {}
    for package in packages:
        stem = f'{package["name"]}-{package["version"]}'
        archive = archives / f"{stem}.crate"
        if not archive.is_file():
            raise SystemExit(f"no packaged archive {archive.name} for the minimal-versions check")
        with tarfile.open(archive, "r:gz") as tar:
            tar.extractall(destination, filter="data")
        unpacked[package["name"]] = destination / stem
    return unpacked


def check_minimal_versions(
    released: list[dict], archives: Path, temporary_path: Path
) -> list[str]:
    """Build each to-be-released crate against its minimal internal versions;
    return one line per checked (crate, feature set)."""
    require_toolchain()
    # Outside every scratch workspace: a path dependency below a workspace
    # root would become a member and drag its dev-dependencies into the
    # resolution.
    patches = unpack_archives(archives, released, temporary_path / "minimal-packages")
    environment = os.environ.copy()
    environment["CARGO_TARGET_DIR"] = str(temporary_path / "minimal-target")
    environment.pop("RUSTUP_TOOLCHAIN", None)
    lines: list[str] = []
    failures: list[str] = []
    for package in released:
        name, version = package["name"], package["version"]
        for label, features in feature_sets(package):
            started = time.monotonic()
            workspace = temporary_path / f"minimal-{name}-{label}"
            for relative, text in probe_workspace_files(name, version, features, patches).items():
                path = workspace / relative
                path.parent.mkdir(parents=True, exist_ok=True)
                path.write_text(text, encoding="utf-8")
            manifest = ["--manifest-path", str(workspace / "Cargo.toml")]
            # Everything at its newest first, then only the internal crates
            # down to the oldest version the requirements admit. External
            # crates keep resolving normally: their floors are not ours.
            result = run_logged([CARGO, "generate-lockfile", *manifest], cwd=ROOT, env=environment)
            if result.returncode != 0:
                raise SystemExit(f"cannot resolve {name} {version}:\n{result.stdout}")
            specs = internal_update_specs((workspace / "Cargo.lock").read_text(encoding="utf-8"))
            result = run_logged(
                minimal_update_command(workspace / "Cargo.toml", specs), cwd=ROOT, env=environment
            )
            if result.returncode != 0:
                raise SystemExit(
                    f"cannot resolve {name} {version} at minimal internal versions:\n{result.stdout}"
                )
            locked = locked_internal_versions((workspace / "Cargo.lock").read_text(encoding="utf-8"))
            # The workspace toolchain builds: the nightly only resolved.
            result = run_logged(
                [CARGO, "check", "--locked", "--quiet", *manifest], cwd=ROOT, env=environment
            )
            if result.returncode != 0:
                failures.append(minimal_versions_failure(name, version, label, result.stdout, locked))
                print(failures[-1].splitlines()[0], flush=True)
                continue
            lowered = sum(1 for row in locked if row[2] and row[0] != name)
            seconds = time.monotonic() - started
            lines.append(
                f"MINIMAL_VERSIONS=PASS {name} {version} ({label}) "
                f"internal={len(locked)} from_registry={lowered} seconds={seconds:.0f}"
            )
            print(lines[-1], flush=True)
    if failures:
        raise SystemExit("\n\n".join(failures))
    return lines


def package_selected(
    packages: list[dict], temporary_path: Path, patch_config: Path, allow_dirty: bool
) -> Path:
    """Package only `packages`, unverified: the minimal-versions build that
    follows compiles them anyway."""
    target = temporary_path / "selected-target"
    command = [CARGO, "package", "--locked", "--quiet", "--exclude-lockfile", "--no-verify"]
    for package in packages:
        command.extend(["-p", package["name"]])
    if allow_dirty:
        command.append("--allow-dirty")
    command.extend(["--config", str(patch_config)])
    environment = os.environ.copy()
    environment["CARGO_TARGET_DIR"] = str(target)
    subprocess.run(command, cwd=ROOT, env=environment, check=True)
    return target / "package"


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--allow-dirty",
        action="store_true",
        help="forward --allow-dirty for local pre-commit verification",
    )
    parser.add_argument(
        "--crates",
        help="comma-separated crates to check at minimal internal versions, instead of "
        "the crates whose version is not on crates.io yet",
    )
    parser.add_argument(
        "--minimal-versions-only",
        action="store_true",
        help="skip the archive preflight; package and check only the minimal-versions set",
    )
    args = parser.parse_args()

    publishable, excluded = workspace_packages(cargo_metadata())
    lookup_started = time.monotonic()
    if args.crates is not None:
        released = select_packages(publishable, [n for n in args.crates.split(",") if n])
    else:
        released = to_be_released(publishable)
    print(
        "MINIMAL_VERSIONS_SET="
        + (",".join(f'{p["name"]}@{p["version"]}' for p in released) or "-")
        + f" lookup_seconds={time.monotonic() - lookup_started:.1f}",
        flush=True,
    )

    with tempfile.TemporaryDirectory(prefix="axiolid-package-preflight-") as temporary:
        temporary_path = Path(temporary)
        patch_config = temporary_path / "workspace-patches.toml"
        write_patch_config(publishable, patch_config)

        roots = None
        if args.minimal_versions_only:
            if not released:
                print("MINIMAL_VERSIONS=PASS crates=0 builds=0 seconds=0")
                return
            archives = package_selected(released, temporary_path, patch_config, args.allow_dirty)
        else:
            archives, roots = bootstrap_preflight(
                publishable, excluded, temporary_path, patch_config, args.allow_dirty
            )
        started = time.monotonic()
        checked = check_minimal_versions(released, archives, temporary_path) if released else []
        minimal_seconds = time.monotonic() - started

    if roots is not None:
        print(
            f"PACKAGE_SOURCE_PREFLIGHT=PASS archives={len(publishable)} "
            f"lockfiles=absent roots_byte_identical={roots} "
            f"excluded={','.join(package['name'] for package in excluded) or '-'}"
        )
    print(
        f"MINIMAL_VERSIONS=PASS crates={len(released)} builds={len(checked)} "
        f"seconds={minimal_seconds:.0f}"
    )


if __name__ == "__main__":
    main()
