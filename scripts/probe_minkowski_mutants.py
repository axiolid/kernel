"""Mutation probe for 2D Minkowski sums by a polygon or region, convex or
not (#145), and disc morphology on a stated side (#163).

Each mutant weakens one step -- the translated region, the edge pieces,
which polygon brackets which side, the convex fast path, the rings' boxes
that prune the arrangement, the convex cut of a non-convex shape (hole
bridges, ears, the tiling certificate, the convex merge), the edge-pair
route taken without a cut, or the erosion's anchor -- and must turn a test
red.

Equivalent mutants, deliberately not listed:
- a complement box with no room about the region: any box holding the
  region serves, since a translate of the polygon leaving the box crosses
  its boundary, whose own edge pieces then cover the point (the room only
  keeps the box's edges off the region's);
- cutting the larger operand of a region sum instead of the smaller: the
  sum is the same, only its cost changes;
- anchoring every erosion on the region moved by `-k0`, also when the shape
  holds the origin: both anchors give the same set up to the vertex sums'
  rounding, the moved rings being the ones the complement's copy already
  uses (the region's own anchor is kept so that those results stay bit
  for bit as before);
- dropping the rounding margin of the disc polygons: it keeps output
  rounding and the presentation's dropping of short edges from carrying a
  point across the exact result, and the tests' samples stay clear of
  boundaries by more than that.
"""
import pathlib, subprocess, sys

ROOT = pathlib.Path(__file__).resolve().parents[1]
M = "crates/algorithms/planar/overlay/src/minkowski.rs"
E = "crates/algorithms/planar/overlay/src/exact_arc/edge.rs"
C = "crates/algorithms/planar/overlay/src/convex_parts.rs"
TESTS = ["-p", "axiolid-overlay", "--test", "minkowski"]
# The crate's unit tests too: the cut and the edge-pair route are checked
# there.
UNIT = ["-p", "axiolid-overlay", "--lib", "--test", "minkowski"]

MUTANTS = [
    ('translated region dropped', M, '        self.sets.iter().any(|(at, set)| set.holds(&flags[*at..]))\n            ||', '        false\n            ||', TESTS),
    ('edge pieces dropped', M, '            || flags[self.pieces..].iter().any(|&f| f)', '            || false', TESTS),
    ('bracketing polygons swapped', M, '        let inscribed = (side == BoundSide::Inner) != erode;', '        let inscribed = (side == BoundSide::Inner) == erode;', TESTS),
    ('non-convex polygon taken as convex', M, '    if !on_hull || hull.len() < 3 {', '    if hull.len() < 3 {', TESTS),
    ('ring boxes miss edges to the left', E, '            x0: a.x0.min(b.x0),', '            x0: a.x0,', TESTS),
    ('only the first convex piece summed', M, '                        for k in convex {', '                        for k in convex.iter().take(1) {', TESTS),
    ('erosion anchored on the region for any shape', M, '        let origin = parts.iter().any(Part::holds_origin);', '        let origin = true;', TESTS),
    ('ring orientation ignored', C, '    if (turn == Sign::Positive) != ccw {', '    if false {', UNIT),
    ('hole bridges cross edges', C, '            if walls.iter().any(|&(p, q)| blocks(h, v, p, q)) {', '            if false {', UNIT),
    ('ears may hold vertices', C, '            && !live.iter().any(|&k| {', '            && !live.iter().take(0).any(|&k| {', UNIT),
    ('tiling certificate accepts anything', C, '    chain.values().all(|&c| c == 0)', '    !chain.is_empty()', UNIT),
    ('merge keeps a non-convex union', C, '    (convex_at(p.len() - 1) && convex_at(0)).then_some(out)', '    Some(out)', UNIT),
    ('edge-pair route misses the shape moved into the region', M, '                sets.push((out.len(), Set::Region(shape.clone())));', '', UNIT),
    ('edge-pair route without parallelograms', M, '                        for r in own {', '                        for r in own.iter().take(0) {', UNIT),
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
