"""Mutation probe for a boolean's deviation measured against its exact
result (#235).

Each mutant loosens a step the measured bound rests on: a flat region's
off-plane term or its covering test, and a trimmed face's domain
classification. Each must turn a test red: `boolean_deviation.rs` samples
the exact results of booleans against the compiled meshes and checks the
reported bound covers them (stopping once within the budget), the
`deviation::boolean` unit test checks a bound tightened to the sampled
maximum against the cut, and the `certify` unit tests probe flat regions
directly. The branch and bound's own steps are probed by
`probe_deviation_report_mutants.py`.
"""
import pathlib, subprocess, sys

ROOT = pathlib.Path(__file__).resolve().parents[1]
FLAT = "crates/execution/compile/src/certify/flat.rs"
BOOL = "crates/execution/compile/src/deviation/boolean.rs"
MEASURED = ["-p", "axiolid-mesh-compile", "--test", "boolean_deviation"]
UNIT = ["-p", "axiolid-mesh-compile", "--lib"]

MUTANTS = [
    ('flat region off-plane term dropped', FLAT,
     '            let bound = off + region.thickness;',
     '            let bound = region.thickness;', [UNIT, MEASURED]),
    ('flat region boundary edges ignored', FLAT,
     '            if inside && !touches(&quad, &region.boundary, region.margin) {',
     '            if inside {', [UNIT, MEASURED]),
    ('flat region joins folded neighbours', FLAT,
     '                    if nu.dot(normal) > 0.5 && level {',
     '                    if level {', [UNIT, MEASURED]),
    ('trimmed face drops its boundary cells', BOOL,
     '        if near {\n            return Coverage::Boundary;',
     '        if near {\n            return Coverage::Outside;', [UNIT, MEASURED]),
    ('trimmed face parity inverted', BOOL,
     '        if self.inside(centre) {',
     '        if !self.inside(centre) {', [UNIT, MEASURED]),
]

def run(target):
    return subprocess.run(
        ["cargo", "test", "-q", *target],
        cwd=ROOT, capture_output=True, text=True, timeout=1800,
    ).returncode

survivors = []
for name, rel, old, new, targets in MUTANTS:
    path = ROOT / rel
    original = path.read_text()
    assert original.count(old) == 1, f"anchor for '{name}' not unique/found"
    path.write_text(original.replace(old, new))
    code = 0
    try:
        for target in targets:
            try:
                code = run(target)
            except subprocess.TimeoutExpired:
                code = -1
            if code != 0:
                break
    finally:
        path.write_text(original)
    status = "killed" if code != 0 else "SURVIVED"
    print(f"{status:8} {name}", flush=True)
    if code == 0:
        survivors.append(name)
print(f"{len(MUTANTS) - len(survivors)}/{len(MUTANTS)} killed")
sys.exit(1 if survivors else 0)
