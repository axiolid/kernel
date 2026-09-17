#!/usr/bin/env python3
"""Crates that must build standalone, derived from `cargo metadata`.

The gate used to carry this list by hand. A new workspace member was
not added automatically, so it silently escaped the isolated-build
check -- the failure kernel#38 was filed about, and which recurred
while adding `axiolid-curve-evaluate-contract`.

A published crate is compiled by whoever depends on it, with only the
features they ask for. Building it alone is the only check that its
own `Cargo.toml` is complete: a dependency reachable solely through a
sibling's feature unification is a defect the workspace build cannot
see.

Unpublished members are excluded deliberately. Nobody can consume
`xtask` or `axiolid-benchmark` standalone, so an isolated build proves
nothing a workspace build has not already proven. Publication status
is the rule rather than a name list, so the exclusion cannot silently
widen.
"""

from __future__ import annotations

import json
import subprocess
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]


def isolated_build_targets() -> list[str]:
    """Publishable workspace members, in dependency-insensitive sorted order."""
    result = subprocess.run(
        ["cargo", "metadata", "--format-version", "1", "--no-deps", "--locked"],
        cwd=ROOT,
        check=True,
        text=True,
        stdout=subprocess.PIPE,
    )
    packages = json.loads(result.stdout)["packages"]
    # `publish = []` means private. Any other value (absent, or a registry
    # allow-list) means the crate can be consumed standalone.
    return sorted(p["name"] for p in packages if p.get("publish") != [])


def main() -> int:
    for name in isolated_build_targets():
        print(name)
    return 0


if __name__ == "__main__":
    sys.exit(main())
