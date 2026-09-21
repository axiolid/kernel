#!/usr/bin/env python3
"""Bootstrap-preflight every publishable workspace crate without uploading.

First releases cannot package dependent crates unpatched before their internal
names exist on crates.io. This preflight therefore compiles normalized,
lock-free source archives with command-scoped local patches. It rejects any
Cargo.lock (and thus patch resolution metadata) in those archives, and proves
byte identity against unpatched archives for every internal-dependency root.

This is deliberately not the upload-byte gate. ``publish-workspace.py`` later
packages and compiles each exact unpatched, lockful archive after its children
are visible, reproduces it through ``cargo publish --dry-run``, and binds the
uploaded bytes to the crates.io checksum.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import os
from pathlib import Path
import subprocess
import tarfile
import tempfile

ROOT = Path(__file__).resolve().parents[1]
CARGO = os.environ.get("CARGO", "cargo")


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


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--allow-dirty",
        action="store_true",
        help="forward --allow-dirty for local pre-commit verification",
    )
    args = parser.parse_args()

    publishable, excluded = workspace_packages(cargo_metadata())
    with tempfile.TemporaryDirectory(prefix="axiolid-package-preflight-") as temporary:
        temporary_path = Path(temporary)
        patch_config = temporary_path / "workspace-patches.toml"
        patched_target = temporary_path / "patched-target"
        rows = ["[patch.crates-io]"]
        for package in publishable:
            package_path = Path(package["manifest_path"]).parent
            rows.append(f'"{package["name"]}" = {{ path = "{package_path}" }}')
        patch_config.write_text("\n".join(rows) + "\n", encoding="utf-8")

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
        if args.allow_dirty:
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
            if args.allow_dirty:
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
            if args.allow_dirty:
                root_command.append("--allow-dirty")
            subprocess.run(root_command, cwd=ROOT, env=root_environment, check=True)
            filename = f'{package["name"]}-{package["version"]}.crate'
            assert_archives_identical(
                patched_target / "package" / filename,
                root_target / "package" / filename,
            )

    print(
        f"PACKAGE_SOURCE_PREFLIGHT=PASS archives={len(publishable)} "
        f"lockfiles=absent roots_byte_identical={len(roots)} "
        f"excluded={','.join(package['name'] for package in excluded) or '-'}"
    )


if __name__ == "__main__":
    main()
