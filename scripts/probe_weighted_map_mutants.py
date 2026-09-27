"""Mutation probe for weighted distance maps (#195).

Each mutant makes a bound unsound -- a lower bound above the distance or
an upper bound below it -- or drops a refusal, and must turn a test red.
The tests that kill them compare brackets against closed forms, and two
maps of one scene at different spacings against each other.

Equivalent mutants, deliberately not listed:

- dropping the rounding margins (a few ulps, far inside every tolerance);
- the rules an optimal walk obeys at a cost edge (no pure reflection, no
  reflection along a piece unless along is cheaper), relaxed: they only
  tighten the lower bound, so dropping one loosens it, which only the
  tightness assertions could see, and those hold at the spacings tested
  for the scenes whose walks never touch a cost edge twice on one side.
"""
import pathlib, subprocess, sys

ROOT = pathlib.Path(__file__).resolve().parents[1]
W = "crates/algorithms/planar/route/src/weighted.rs"
TESTS = ["-p", "axiolid-route", "--test", "weighted", "--lib"]

MUTANTS = [
    ('hull held by corners alone', W, '            if meets_open_convex(p, q, hull, n)? {', '            if false && meets_open_convex(p, q, hull, n)? {', TESTS),
    ('leaving takes the dearer side always', W, '        (Sign::Positive, Sign::Positive) => left,\n        (Sign::Negative, Sign::Negative) => right,\n        _ => left.min(right),', '        _ => left.max(right),', TESTS),
    ('transmission forbidden, not reflection', W, '            if tag.side == edge.arrive {', '            if tag.side != edge.arrive && tag.side <= ALONG {', TESTS),
    ('a hull across a notch held', W, '    Ok(!on_boundary && in_polygon(polygon, centre)?)', '    let _ = on_boundary;\n    Ok(true)', TESTS),
    ('a touch counted with an end on the line', W, '    if side(p, q, u)? == Sign::Zero || side(p, q, v)? == Sign::Zero {\n        return Ok(false);\n    }', '', TESTS),
    ('collinear hop at its full length', W, '            return Ok(Some(self.weights.segment(u, v, self.scale)?.0));', '            return Ok(Some(self.weights.segment(a1, b1, self.scale)?.0.max(self.weights.segment(a2, b2, self.scale)?.0)));', TESTS),
    ('a doubtful stretch at the greatest factor below', W, '                lower += stretch;', '                lower += self.greatest * stretch;', TESTS),
    ('a doubtful stretch at factor 1 above', W, '                upper += self.greatest * stretch;', '                upper += stretch;', TESTS),
    ('both sides of a wall free', W, '        Ok(if walled { (left, right) } else { (true, true) })', '        Ok((true, true))', TESTS),
    ('no points on any cost edge', W, '        if weights.free_sides(p, q)? != (true, true) {', '        if true {', TESTS),
    ('every interval end a vertex', W, '                a_vertex: k == 0,', '                a_vertex: true,', TESTS),
    ('cost edges do not block a cell hop', W, '        for piece in &self.weights.pieces {\n            if blocks(piece.p, piece.q)? {\n                return Ok(None);\n            }\n        }', '', TESTS),
    ('any end touching a blocker blocks', W, '!((skip_u || skip_v) && through_end(u, v, p, q)?)', '!through_end(u, v, p, q)?', TESTS),
    ('overlap takes the lesser factor', W, '                    if in_polygon(polygon, m)? {\n                        left = left.max(extra);\n                        right = right.max(extra);', '                    if in_polygon(polygon, m)? {\n                        left = if left == 0.0 { extra } else { left.min(extra) };\n                        right = if right == 0.0 { extra } else { right.min(extra) };', TESTS),
    ('farthest point ignores the factor', W, '        let k = steep[root];', '        let k = 1.0;', TESTS),
    ('factor below 1 accepted', W, '            return Err(MapError::InvalidFactor { index });', '            let _ = index;', TESTS),
    ('any spacing accepted', W, '    if !(spacing.is_finite() && spacing > 0.0) {', '    if false {', TESTS),
]

def run(target):
    return subprocess.run(
        ["cargo", "test", "-q", *target],
        cwd=ROOT, capture_output=True, text=True, timeout=1800,
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
