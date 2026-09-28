"""Mutation probe for the shortest walk forced through a region (#196).

Each mutant weakens the bracket -- the Lipschitz radius, the vertex-pair
bound, the proofs that a vertex is hidden from a cell, the floor of the
shortest walk -- or the handling of mismatched maps and unreached parts,
and must turn a test red.

Equivalent mutants, deliberately not listed:

- dropping the widening for the rounding of lengths, a few ulps of the
  length, far inside the tolerance;
- dropping the floor of the shortest walk overall. It is sound and cheap,
  but the vertex-pair bound reaches the shortest walk on every cell the
  walk passes through, so the floor never binds in a closed-form test.
"""
import pathlib, subprocess, sys

ROOT = pathlib.Path(__file__).resolve().parents[1]
F = "crates/algorithms/planar/route/src/forced.rs"
TESTS = ["-p", "axiolid-route", "--test", "forced_walk", "--lib"]

MUTANTS = [
    ('Lipschitz constant 1, not 2', F, '    f - 2.0 * radius\n', '    f - radius\n', TESTS),
    ('inherited anchor without radius', F, '(lipschitz_lower(f, radius(a)), (a, f))', '(f, (a, f))', TESTS),
    ('pair bound adds both terms', F, 'best = best.min(du + dv + (u - v).length().max(gu + gv));', 'best = best.min(du + dv + (u - v).length() + gu + gv);', TESTS),
    ('pair bound unused', F, '(bound.max(pairs) - slack(depth, radius(anchor.0)))', '(bound - slack(depth, radius(anchor.0)))', TESTS),
    ('pair search stops too early', F, '                if s + r >= best {', '                if s + r >= best - 1.0 {', TESTS),
    ('hidden when one corner is blocked', F, '            if crosses(v, t[0], p, q)? && crosses(v, t[1], p, q)? && crosses(v, t[2], p, q)? {', '            if crosses(v, t[0], p, q)? {', TESTS),
    ('no obstacle hides a vertex', F, '            if crosses(v, t[0], p, q)? && crosses(v, t[1], p, q)? && crosses(v, t[2], p, q)? {', '            if false {', TESTS),
    ('solid wedge on the wrong side', F, '(small_is_inside == hole).then_some', '(small_is_inside != hole).then_some', TESTS),
    ('solid wedge from one corner', F, '            if inside(t[0])? && inside(t[1])? && inside(t[2])? {\n                continue \'vertex;', '            if inside(t[0])? {\n                continue \'vertex;', TESTS),
    ('upper bound from points outside the polygon', F, '            if crate::map::in_polygon(self.through, a)? && self.upper', '            if self.upper', TESTS),
    ('mismatched maps accepted', F, ') -> Result<ForcedWalk, FarthestError> {\n    if !(tolerance.is_finite() && tolerance >= 0.0) {\n        return Err(FarthestError::InvalidTolerance);\n    }\n    if !from.same_space(to) {', ') -> Result<ForcedWalk, FarthestError> {\n    if !(tolerance.is_finite() && tolerance >= 0.0) {\n        return Err(FarthestError::InvalidTolerance);\n    }\n    if false {', TESTS),
    ('an unreached part is an error', F, '            // No walk reaches this part of the free space.\n            Anchored::Unreached => {}', '            Anchored::Unreached => return Err(FarthestError::Triangulation),', TESTS),
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
