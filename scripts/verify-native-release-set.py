#!/usr/bin/env python3
"""Verify the complete supported native release set and emit SHA256SUMS."""

from __future__ import annotations

import argparse
import importlib.util
import json
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
EXPECTED_TARGETS = {
    "x86_64-unknown-linux-gnu",
    "aarch64-unknown-linux-gnu",
    "x86_64-apple-darwin",
    "x86_64-pc-windows-msvc",
    "aarch64-pc-windows-msvc",
}


def load_verifier():
    path = ROOT / "scripts/verify-native-package.py"
    spec = importlib.util.spec_from_file_location("verify_native_package", path)
    assert spec and spec.loader
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("directory", type=Path)
    parser.add_argument(
        "--source-commit", help="Required exact commit for every archive"
    )
    parser.add_argument(
        "--version",
        help="Required package version for every archive (the release tag without 'v')",
    )
    args = parser.parse_args()
    verify = load_verifier()
    archives = sorted(
        [
            *args.directory.glob("axiolid-native-*.tar.gz"),
            *args.directory.glob("axiolid-native-*.zip"),
        ]
    )
    found: dict[str, Path] = {}
    commits: set[str] = set()
    try:
        if EXPECTED_TARGETS != verify.SUPPORTED_TARGETS:
            raise ValueError(
                "release target policy differs between packager and set verifier"
            )
        for archive in archives:
            verify.verify_checksum(archive)
            files = verify.read_archive(archive)
            root, manifest = verify.verify_files(files, allow_dirty=False)
            verify.verify_binaries(files, root, manifest)
            target = manifest["target"]
            if manifest["profile"] != "release":
                raise ValueError(
                    f"release target {target} was built with {manifest['profile']} profile"
                )
            commits.add(manifest["source_commit"])
            if args.version and manifest["package_version"] != args.version:
                raise ValueError(
                    f"{archive.name} is package version {manifest['package_version']}, "
                    f"release is {args.version}"
                )
            if target in found:
                raise ValueError(f"duplicate release target: {target}")
            found[target] = archive
        if set(found) != EXPECTED_TARGETS:
            raise ValueError(
                f"release target set differs: missing={sorted(EXPECTED_TARGETS - set(found))}, "
                f"unexpected={sorted(set(found) - EXPECTED_TARGETS)}"
            )
        if len(commits) != 1:
            raise ValueError(
                f"release archives disagree on source commit: {sorted(commits)}"
            )
        if args.source_commit and commits != {args.source_commit}:
            raise ValueError(
                f"release archives were built from {sorted(commits)}, expected {args.source_commit}"
            )
        lines = [
            f"{verify.digest_file(path)}  {path.name}"
            for path in sorted(found.values())
        ]
        (args.directory / "SHA256SUMS").write_text("\n".join(lines) + "\n")
    except (OSError, ValueError, json.JSONDecodeError) as error:
        print(f"native release set verification failed: {error}", file=sys.stderr)
        return 1
    print(f"verified native release set: {len(found)} targets")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
