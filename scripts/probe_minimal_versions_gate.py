#!/usr/bin/env python3
"""Prove the minimal-internal-versions release check can fail (kernel#337).

The first mutant restates the `axiolid-mesh-compile` 0.3.28 situation: its
floors on internal crates sit below what it uses, so the build against the
lowest versions its manifest admits must fail, and
`scripts/verify-packages.py` must reject it. The floors to lower are chosen at
run time from `cargo metadata` and the crates.io sparse index (see
`stale_floor_edits`), so the probe does not depend on which crates are pending
release. When every dependency is pending release there is no published
floor to make stale, so it reports the probe as inapplicable and exits 0; any
other way of finding no floor fails loudly.
TODO(#449): make the probe hermetic with a fixture registry.

The others break the check itself while that stale floor is in place, and the
check must then still reject it: each is run through `REJECTS`, which passes
only when the check fails with `MINIMAL_VERSIONS=FAIL` naming
`axiolid-mesh-compile` as the crate that does not build. A mutant that lets the
stale floor through therefore fails the command and is killed. The pure
helpers have their own unit tests in `scripts/test_release_scripts.py`.

Needs the pinned nightly (`MINIMAL_VERSIONS_TOOLCHAIN`) and crates.io.
"""

from __future__ import annotations

import json
import os
from pathlib import Path
import subprocess
import sys
import urllib.request

import probe_lib

ROOT = Path(__file__).resolve().parent.parent
MANIFEST = ROOT / "crates/execution/compile/Cargo.toml"
SCRIPT = ROOT / "scripts/verify-packages.py"
PROBE_CRATE = "axiolid-mesh-compile"
USER_AGENT = "axiolid-release-probe (https://github.com/axiolid/kernel)"

CHECK = [
    sys.executable,
    str(SCRIPT),
    "--minimal-versions-only",
    "--allow-dirty",
    "--crates",
    PROBE_CRATE,
]
REJECTS = [
    sys.executable,
    "-c",
    "import re, subprocess, sys\n"
    "r = subprocess.run(sys.argv[1:], text=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT)\n"
    "print(r.stdout[-2000:])\n"
    "named = re.search(r'MINIMAL_VERSIONS=FAIL .*axiolid-mesh-compile \\S+ \\(to be released\\) "
    "does not build', r.stdout)\n"
    "sys.exit(0 if r.returncode != 0 and named else 1)\n",
    *CHECK,
]


def sparse_index_path(name: str) -> str:
    lowered = name.lower()
    if len(lowered) <= 2:
        return f"{len(lowered)}/{lowered}"
    if len(lowered) == 3:
        return f"3/{lowered[0]}/{lowered}"
    return f"{lowered[:2]}/{lowered[2:4]}/{lowered}"


def index_versions(name: str) -> dict[str, bool]:
    """Every version of `name` on crates.io, mapped to whether it is yanked."""
    request = urllib.request.Request(
        f"https://index.crates.io/{sparse_index_path(name)}",
        headers={"User-Agent": USER_AGENT},
    )
    with urllib.request.urlopen(request, timeout=30) as response:
        rows = [json.loads(line) for line in response.read().decode().splitlines() if line.strip()]
    return {row["vers"]: bool(row["yanked"]) for row in rows}


def version_key(version: str) -> tuple[int, ...]:
    return tuple(int(part) for part in version.split("."))


def caret_compatible(older: str, current: str) -> bool:
    a, b = version_key(older), version_key(current)
    return a[0] == b[0] and (a[0] != 0 or a[1] == b[1])


class Inapplicable(Exception):
    """Every candidate dependency is pending release; carries their names."""


def stale_floor_edits(metadata: dict, manifest: Path, lookup=index_versions) -> list[tuple[str, str]]:
    """The (old, new) manifest lines that lower every floor of the probe crate
    that can be lowered.

    A dependency qualifies when it is a workspace member published at its
    current version (a crate pending release is built from the tree, so its
    floor is not lowered) and has an older, unyanked, caret-compatible
    published version. Its `name.workspace = true` line becomes a path
    dependency requiring the oldest such version.

    Raises Inapplicable when every candidate dependency is pending release:
    the check then has no published floor to test. Any other way of finding no
    floor (no dependency, no older version, registry unreachable) is an error.
    """
    members = {package["name"]: package for package in metadata["packages"]}
    text = manifest.read_text()
    edits = []
    pending = []
    candidates = 0
    for dependency in members[PROBE_CRATE]["dependencies"]:
        name = dependency["name"]
        line = f"{name}.workspace = true\n"
        member = members.get(name)
        if member is None or dependency["kind"] is not None or line not in text:
            continue
        candidates += 1
        versions = lookup(name)
        current = member["version"]
        if current not in versions:
            pending.append(name)
            continue
        older = sorted(
            (
                version
                for version, yanked in versions.items()
                if not yanked
                and version_key(version) < version_key(current)
                and caret_compatible(version, current)
            ),
            key=version_key,
        )
        if not older:
            continue
        path = os.path.relpath(Path(member["manifest_path"]).parent, manifest.parent)
        edits.append((line, f'{name} = {{ path = "{path}", version = "{older[0]}" }}\n'))
    if not edits and candidates and len(pending) == candidates:
        raise Inapplicable(pending)
    if not edits:
        raise SystemExit(
            "no published dependency to probe with: no dependency of "
            f"{PROBE_CRATE} is published at its current version with an older "
            "published version"
        )
    return edits


def build_mutants(stale_floor: list[tuple[str, str]]) -> list:
    def broken_check(name: str, old: str, new: str) -> probe_lib.Mutant:
        edits = [(MANIFEST, a, b) for a, b in stale_floor] + [(SCRIPT, old, new)]
        return probe_lib.Mutant(name, edits, [REJECTS])

    return [
        ("stale re-exported floor is released", MANIFEST, stale_floor, [CHECK]),
        broken_check(
            "only direct dependencies are lowered",
            '        "minimal-versions",\n        "update",\n',
            '        "direct-minimal-versions",\n        "update",\n',
        ),
        broken_check(
            "published internal crates are not lowered",
            '        for row in lock.get("package", [])\n'
            '        if is_internal(row["name"])\n'
            "    )\n\n\ndef locked_internal_versions",
            '        for row in lock.get("package", [])\n'
            '        if is_internal(row["name"]) and row.get("source") is None\n'
            "    )\n\n\ndef locked_internal_versions",
        ),
        broken_check(
            "a failed build is not reported",
            "            if result.returncode != 0:\n                failures.append(",
            "            if False:\n                failures.append(",
        ),
        broken_check(
            "the failure does not name the crate that does not build",
            "    failing = failing_crates(output) or [name]\n",
            "    failing = [PROBE_PACKAGE]\n",
        ),
    ]


def mutants_to_probe() -> list:
    metadata = json.loads(
        subprocess.check_output(
            ["cargo", "metadata", "--no-deps", "--format-version", "1"], cwd=ROOT, text=True
        )
    )
    try:
        floor_edits = stale_floor_edits(metadata, MANIFEST)
    except Inapplicable as pending:
        print(
            f"minimal-versions probe inapplicable: every dependency of {PROBE_CRATE} is "
            f"pending release ({', '.join(pending.args[0])}); the check has no "
            "published floor to test"
        )
        sys.exit(0)
    mutants = build_mutants(floor_edits)

    # The clean tree must pass first, or a mutant could be "killed" by a broken
    # environment (no nightly, no registry). probe_lib's own baseline would also
    # run REJECTS, which is red on a clean tree by design.
    if not {"--list", "--check-anchors"} & set(sys.argv):
        code, output = probe_lib.run_command(CHECK, cwd=ROOT, timeout=1800, env=None)
        if code != 0:
            print(output[-3000:])
            sys.exit("baseline: FAIL (the clean tree must pass the check) -- fix it before probing")
        print("baseline: PASS", flush=True)

    return mutants


# Importing the module (the unit tests do) must not run the probe.
probe_lib.main(mutants_to_probe(), command=(), timeout=1800) if __name__ == "__main__" else None
