"""Mutation probe for the minimum-area rectangle (#182).

Each mutant skips a rotating-calipers step, an edge, the exact area
comparison or the tie rule, and must turn a test red.

Equivalent mutant, deliberately not listed: advancing a caliper only on a
strict increase instead of a non-decrease. On a tie the two vertices are
equally extreme, so the spans, and the area, are the same.
"""
import pathlib, subprocess, sys

ROOT = pathlib.Path(__file__).resolve().parents[1]
R = "crates/algorithms/planar/overlay/src/rectangle.rs"
TESTS = ["-p", "axiolid-overlay", "--test", "rectangle"]

MUTANTS = [
    ('far caliper not advanced', R, '            hi = advance(hi, d, false, true);', '', TESTS),
    ('top caliper not advanced', R, '            top = advance(top, d, true, true);', '', TESTS),
    ('near caliper not advanced', R, '            lo = advance(lo, d, false, false);', '', TESTS),
    ('calipers move one step at most', R, '        for _ in 0..n {\n            if !step(', '        for _ in 0..1 {\n            if !step(', TESTS),
    ('half the hull edges tried', R, '    for base in 0..n {', '    for base in 0..n.div_ceil(2) {', TESTS),
    ('larger area preferred', R, '            Some(Sign::Negative) => {\n                best = *c;', '            Some(Sign::Positive) => {\n                best = *c;', TESTS),
    ('ties keep the first found', R, '                if clockwise(canonical(best.d), axis) {', '                if false && clockwise(canonical(best.d), axis) {', TESTS),
    ('tied orientations not counted', R, '                    ties.push(axis);', '', TESTS),
    ('collinear hull vertices kept', R, '                if sign(&Orient { a, b, c }) == Sign::Positive {', '                if sign(&Orient { a, b, c }) != Sign::Negative {', TESTS),
    ('segment width left rounded', R, '            rectangle.half_extents[usize::from(axis == d || axis == -d)] = 0.0;', '', TESTS),
]

def run(target):
    return subprocess.run(
        ["cargo", "test", "-q", *target],
        cwd=ROOT, capture_output=True, text=True, timeout=1200,
    ).returncode

survivors = []
for name, rel, old, new, target in MUTANTS:
    path = ROOT / rel
    original = path.read_text()
    assert original.count(old) == 1, f"anchor for '{name}' not unique/found"
    path.write_text(original.replace(old, new))
    try:
        try:
            code = run(target)
        except subprocess.TimeoutExpired:
            code = -1
    finally:
        path.write_text(original)
    status = "killed" if code != 0 else "SURVIVED"
    print(f"{status:8} {name}", flush=True)
    if code == 0:
        survivors.append(name)
print(f"{len(MUTANTS) - len(survivors)}/{len(MUTANTS)} killed")
sys.exit(1 if survivors else 0)
