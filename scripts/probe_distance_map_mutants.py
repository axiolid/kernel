"""Mutation probe for distance maps and the farthest point (#186).

Each mutant weakens the bracket -- drops the Lipschitz radius, the
anchors, the split, the proof that an anchor lies in its triangle -- or
the multi-source search, and must turn a test red.

Equivalent mutant, deliberately not listed: dropping the widening of the
interval for the rounding of lengths. It is a few ulps of the distance,
far inside the cell slack and the tolerance, so no closed-form test can
see it; it is kept for the proof, not for a test.
"""
import pathlib, subprocess, sys

ROOT = pathlib.Path(__file__).resolve().parents[1]
M = "crates/algorithms/planar/route/src/map.rs"
TESTS = ["-p", "axiolid-route", "--test", "distance_map", "--lib"]

MUTANTS = [
    ('no Lipschitz radius', M, '            let bound = d + radius(a);', '            let bound = d;', TESTS),
    ('parent anchor bound without radius', M, '(d + radius(a), (a, d))', '(d, (a, d))', TESTS),
    ('cells never split', M, '        for corners in [[a, ab, ca], [ab, b, bc], [ca, bc, c], [ab, bc, ca]] {', '        for corners in [[a, ab, ca], [ab, b, bc], [ca, bc, c], [ab, bc, ca]].into_iter().take(0) {', TESTS),
    ('middle child dropped', M, '[ca, bc, c], [ab, bc, ca]] {', '[ca, bc, c]] {', TESTS),
    ('anchors on barriers admitted', M, '        if side(p, q, a)? == Sign::Zero && within(p, q, a) {\n            return Ok(false);', '        if false && side(p, q, a)? == Sign::Zero && within(p, q, a) {\n            return Ok(false);', TESTS),
    ('anchors not proven in the triangle', M, '    if !in_triangle(root, a)? {\n        return Ok(false);\n    }', '', TESTS),
    ('lower bound from points outside the subregion', M, '            if in_polygon(self.subregion, a)? && self.lower', '            if self.lower', TESTS),
    ('unreachable anchors skipped', M, '            let Some(d) = self.distance(a)? else {\n                return Err(FarthestError::Unreachable {\n                    triangle: self.roots[root],\n                });\n            };', '            let Some(d) = self.distance(a)? else {\n                continue;\n            };', TESTS),
    ('subregion edges ignored', M, '            if meets_triangle(p, q, t)? {\n                return Ok(false);', '            if false && meets_triangle(p, q, t)? {\n                return Ok(false);', TESTS),
    ('only the first target seeded', M, '        if let Some(t) = targets.iter().position(|t| t == node) {', '        if let Some(t) = targets.iter().take(1).position(|t| t == node) {', TESTS),
    ('query ignores visibility', M, '                if !ok[k] {\n                    continue;\n                }', '', TESTS),
    ('crossing barriers accepted', M, '                return Err(FarthestError::CrossingObstacles);', '', TESTS),
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
