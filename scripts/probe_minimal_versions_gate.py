#!/usr/bin/env python3
"""Prove the minimal-internal-versions release check can fail (kernel#337).

The first mutant restates the `axiolid-mesh-compile` 0.3.28 situation on the
current tree: `axiolid-mesh-compile` requires `axiolid-reference` 0.3.8 (which
admits `axiolid-evaluate` 0.3.12) and `axiolid-brep-boolean` 0.1.10 (0.3.16),
while it uses `axiolid-evaluate` 0.3.17's frame laws through `axiolid-reference`.
`scripts/verify-packages.py` must reject it.

The others break the check itself while that stale floor is in place, and the
check must then still reject it: each is run through `REJECTS`, which passes
only when the check fails with `MINIMAL_VERSIONS=FAIL` naming
`axiolid-mesh-compile` as the crate that does not build. A mutant that lets the
stale floor through therefore fails the command and is killed. The pure
helpers have their own unit tests in `scripts/test_release_scripts.py`.

Needs the pinned nightly (`MINIMAL_VERSIONS_TOOLCHAIN`) and crates.io.
"""

from pathlib import Path
import sys

import probe_lib

ROOT = Path(__file__).resolve().parent.parent
MANIFEST = ROOT / "crates/execution/compile/Cargo.toml"
SCRIPT = ROOT / "scripts/verify-packages.py"

CHECK = [
    sys.executable,
    str(SCRIPT),
    "--minimal-versions-only",
    "--allow-dirty",
    "--crates",
    "axiolid-mesh-compile",
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

# Published versions, so they stay valid after later releases.
STALE_FLOOR = [
    (
        "axiolid-reference.workspace = true\n",
        'axiolid-reference = { path = "../../algorithms/reference", version = "0.3.8" }\n',
    ),
    (
        "axiolid-brep-boolean.workspace = true\n",
        'axiolid-brep-boolean = { path = "../../algorithms/construction/brep-boolean", '
        'version = "0.1.10" }\n',
    ),
]


def broken_check(name: str, old: str, new: str) -> probe_lib.Mutant:
    edits = [(MANIFEST, a, b) for a, b in STALE_FLOOR] + [(SCRIPT, old, new)]
    return probe_lib.Mutant(name, edits, [REJECTS])


MUTANTS = [
    ("stale re-exported floor is released", MANIFEST, STALE_FLOOR, [CHECK]),
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

# The clean tree must pass first, or a mutant could be "killed" by a broken
# environment (no nightly, no registry). probe_lib's own baseline would also
# run REJECTS, which is red on a clean tree by design.
if not {"--list", "--check-anchors"} & set(sys.argv):
    code, output = probe_lib.run_command(CHECK, cwd=ROOT, timeout=1800, env=None)
    if code != 0:
        print(output[-3000:])
        sys.exit("baseline: FAIL (the clean tree must pass the check) -- fix it before probing")
    print("baseline: PASS", flush=True)

probe_lib.main(MUTANTS, command=(), timeout=1800)
