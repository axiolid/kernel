#!/usr/bin/env python3
"""What `scripts/gate.sh` needs to know about the change under test (#443).

Three questions, all answered against the merge base with `origin/main`:

    changed    the paths that differ (committed, staged, unstaged, untracked)
    infra      which of them are gate infrastructure; exit 0 when there is
               one, 1 when there is none, 2 when it cannot be decided
    packages   the workspace packages the branch tier must check: those that
               own a changed file, plus every workspace package that depends
               on them, directly or transitively (any dependency kind)

Fail closed: an unknown base, a failing `git` or `cargo metadata` is an error
(exit 2), never an empty answer. `gate.sh` reads exit 2 from `infra` as
"infrastructure changed", so a guess can only run more checks, never fewer.

`--files-from FILE` replaces the git query with one path per line, for tests.
"""
from __future__ import annotations

import argparse
import fnmatch
import json
import subprocess
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
BASE_REF = "origin/main"

# A change under one of these runs the full gate, and the gate's own
# mutation probes: the branch tier cannot vouch for it.
INFRA_PREFIXES = (
    "scripts/",
    "tools/xtask/",
    "crates/facade/axiolid-capi/",
    "native/",
    "tests/consumers/",
    "tests/downstream/",
    ".github/",
    ".cargo/",
)
INFRA_GLOBS = ("docs/architecture/*.toml",)
INFRA_NAMES = ("Cargo.toml", "Cargo.lock", "rust-toolchain.toml")


class GateAffectedError(Exception):
    pass


def _run(argv: list[str]) -> str:
    result = subprocess.run(argv, cwd=ROOT, capture_output=True, text=True)
    if result.returncode != 0:
        raise GateAffectedError(f"{' '.join(argv)}: {result.stderr.strip()}")
    return result.stdout


def merge_base(ref: str = BASE_REF) -> str:
    return _run(["git", "merge-base", "HEAD", ref]).strip()


def changed_files(base: str) -> list[str]:
    """Paths that differ from `base`, in the working tree (sorted, unique)."""
    tracked = _run(["git", "diff", "--name-only", "--no-renames", base]).splitlines()
    untracked = _run(["git", "ls-files", "--others", "--exclude-standard"]).splitlines()
    return sorted({path for path in tracked + untracked if path})


def is_infra(path: str) -> bool:
    if path.startswith(INFRA_PREFIXES):
        return True
    if path.rsplit("/", 1)[-1] in INFRA_NAMES:
        return True
    return any(fnmatch.fnmatchcase(path, glob) for glob in INFRA_GLOBS)


def infra_files(files: list[str]) -> list[str]:
    return [path for path in files if is_infra(path)]


def workspace_metadata() -> dict:
    return json.loads(
        _run(["cargo", "metadata", "--no-deps", "--locked", "--offline", "--format-version", "1"])
    )


def affected_packages(metadata: dict, files: list[str]) -> list[str]:
    """Owners of `files` plus all workspace packages depending on them."""
    root = Path(metadata["workspace_root"])
    owners: dict[str, str] = {}  # package directory (relative) -> name
    deps: dict[str, set[str]] = {}
    names = {package["name"] for package in metadata["packages"]}
    for package in metadata["packages"]:
        directory = Path(package["manifest_path"]).parent.relative_to(root).as_posix()
        owners[directory] = package["name"]
        deps[package["name"]] = {
            dep["name"] for dep in package["dependencies"] if dep["name"] in names
        }

    affected: set[str] = set()
    for path in files:
        # The longest directory wins, so a nested package owns its own files.
        best = None
        for directory in owners:
            inside = directory == "." or path == directory or path.startswith(directory + "/")
            if inside and (best is None or len(directory) > len(best)):
                best = directory
        if best is not None and best != ".":
            affected.add(owners[best])

    # Reverse dependents, to a fixed point.
    grew = True
    while grew:
        grew = False
        for name, wants in deps.items():
            if name not in affected and wants & affected:
                affected.add(name)
                grew = True
    return sorted(affected)


def main(argv: list[str]) -> int:
    parser = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    parser.add_argument("command", choices=("changed", "infra", "packages"))
    parser.add_argument("--base", default=None, help=f"merge-base ref (default {BASE_REF})")
    parser.add_argument("--files-from", type=Path, help="read the changed paths from FILE")
    args = parser.parse_args(argv)
    try:
        if args.files_from is not None:
            files = [p for p in args.files_from.read_text().splitlines() if p]
        else:
            files = changed_files(merge_base(args.base or BASE_REF))
        if args.command == "changed":
            print("\n".join(files))
        elif args.command == "infra":
            found = infra_files(files)
            print("\n".join(found))
            return 0 if found else 1
        else:
            print("\n".join(affected_packages(workspace_metadata(), files)))
    except (GateAffectedError, OSError, ValueError, KeyError) as error:
        print(f"gate-affected: {error}", file=sys.stderr)
        return 2
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
