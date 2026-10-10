"""Mutation probe for 2D Minkowski sums by a polygon or region, convex or
not (#145), and disc morphology on a stated side (#163).

Each mutant weakens one step -- the translated region, the edge pieces,
which polygon brackets which side, the convex fast path, the rings' boxes
that prune the arrangement, the convex cut of a non-convex shape (hole
bridges, ears, the tiling certificate, the convex merge), the edge-pair
route taken without a cut, the erosion's anchor, or the sparse rule that
reads a face from the rings holding it (#292: a polygon's outer ring
without its holes, the span of each translated set, the pieces counted
from their first ring, the rings holding each side of a piece) -- and
must turn a test red.

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
  boundaries by more than that;
- reading the polygon after an outer ring as one of its holes
  (`r - offset > end`): the polygons of a set are disjoint, so no side of
  a piece lies in two of their outer rings;
- a complement read without its frame: the rule asks a set only when one
  of its rings holds the point, and the frame holds every ring of its
  copy.
"""
import os, pathlib, signal, subprocess, sys

ROOT = pathlib.Path(__file__).resolve().parents[1]
M = "crates/algorithms/planar/overlay/src/minkowski.rs"
E = "crates/algorithms/planar/overlay/src/exact_arc/edge.rs"
C = "crates/algorithms/planar/overlay/src/convex_parts.rs"
A = "crates/algorithms/planar/overlay/src/arrangement.rs"
TESTS = ["-p", "axiolid-overlay", "--test", "minkowski"]
# The crate's unit tests too: the cut and the edge-pair route are checked
# there.
UNIT = ["-p", "axiolid-overlay", "--lib", "--test", "minkowski"]
# The sparse rule (#292): the unit test compares it with the flags on every
# piece, and the arrangement's public flags are read through it too.
SPARSE = ["-p", "axiolid-overlay", "--lib", "--test", "minkowski", "--test", "arrangement"]

MUTANTS = [
    ('translated region dropped', M, '            if set.holds(&rings[i..next], offset + first) {', '            if false {', TESTS),
    ('edge pieces dropped', M, '        if rings.last().is_some_and(|&r| r - offset >= self.pieces) {', '        if false {', TESTS),
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
    ('holes of a polygon ignored', M, '            ring == self.starts[polygon] && rings.get(i + 1).is_none_or(|&r| r - offset >= end)', '            ring == self.starts[polygon]', SPARSE),
    ('a piece counted from the one after the first', M, 'is_some_and(|&r| r - offset >= self.pieces)', 'is_some_and(|&r| r - offset > self.pieces)', SPARSE),
    ('only the first translated set asked', M, '            i = next;', '            break;', SPARSE),
    ('a set spans one ring too few', M, '            let end = offset + first + set.span();', '            let end = offset + first + set.span() - 1;', SPARSE),
    ('the frame not counted in a complement', M, '            Self::Outside(shape) => 1 + shape.len,', '            Self::Outside(shape) => shape.len,', SPARSE),
    ('a carrier holds both sides', A, '            if on_left == left {', '            if true {', SPARSE),
    ('a carrier read on the wrong side', A, '            Some(&(_, on_left)) => on_left == left,', '            Some(&(_, on_left)) => on_left,', SPARSE),
    ('rings inside a piece dropped', A, '                inside: std::mem::take(&mut edge.inside),', '                inside: Vec::new(),', SPARSE),
    ('a carrier listed out of order', A, '                let at = out.partition_point(|&r| r < ring);', '                let at = out.len();', SPARSE),
    ('flags left set between pieces', A, '                flags[ring] = false;', '                let _ = ring;', SPARSE),
]

# A mutant may make a test loop forever (the convex merge keeping a
# non-convex union did). Past this many seconds it counts as killed, by
# timeout, and its whole process group is stopped.
TIMEOUT = 600


def run(target):
    """The test run's exit code, or `None` when it timed out."""
    process = subprocess.Popen(
        ["cargo", "test", "-q", *target],
        cwd=ROOT, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL,
        start_new_session=True,
    )
    try:
        return process.wait(timeout=TIMEOUT)
    except subprocess.TimeoutExpired:
        os.killpg(process.pid, signal.SIGKILL)
        process.wait()
        return None

survivors = []
for name, rel, old, new, target in MUTANTS:
    path = ROOT / rel
    original = path.read_text()
    assert original.count(old) == 1, f"anchor for '{name}' not unique/found"
    path.write_text(original.replace(old, new))
    try:
        code = run(target)
    finally:
        path.write_text(original)
    if code is None:
        status = "killed (timeout)"
    else:
        status = "killed" if code != 0 else "SURVIVED"
    print(f"{status:8} {name}", flush=True)
    if code == 0:
        survivors.append(name)
print(f"{len(MUTANTS) - len(survivors)}/{len(MUTANTS)} killed")
sys.exit(1 if survivors else 0)
