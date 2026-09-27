"""Mutation probe for 2D Minkowski sums by a convex polygon (#145) and disc
morphology on a stated side (#163).

Each mutant weakens one step -- the translated region, the edge pieces,
which polygon brackets which side, the convexity refusal, or the rings'
boxes that prune the arrangement -- and must turn a test red.

Equivalent mutants, deliberately not listed:
- a complement box with no room about the region: any box holding the
  region serves, since a translate of the polygon leaving the box crosses
  its boundary, whose own edge pieces then cover the point (the room only
  keeps the box's edges off the region's);
- dropping the rounding margin of the disc polygons: it keeps output
  rounding and the presentation's dropping of short edges from carrying a
  point across the exact result, and the tests' samples stay clear of
  boundaries by more than that.
"""
import pathlib, subprocess, sys

ROOT = pathlib.Path(__file__).resolve().parents[1]
M = "crates/algorithms/planar/overlay/src/minkowski.rs"
E = "crates/algorithms/planar/overlay/src/exact_arc/edge.rs"
TESTS = ["-p", "axiolid-overlay", "--test", "minkowski"]

MUTANTS = [
    ('translated region dropped', M, '            |flags| member(&shape, &flags[..n]) || flags[n..].iter().any(|&f| f),', '            |flags| flags[n..].iter().any(|&f| f),', TESTS),
    ('edge pieces dropped', M, '            |flags| member(&shape, &flags[..n]) || flags[n..].iter().any(|&f| f),', '            |flags| member(&shape, &flags[..n]),', TESTS),
    ('bracketing polygons swapped', M, '        let inscribed = (side == BoundSide::Inner) != erode;', '        let inscribed = (side == BoundSide::Inner) == erode;', TESTS),
    ('non-convex polygon accepted', M, '    if !on_hull || hull.len() < 3 {', '    if hull.len() < 3 {', TESTS),
    ('ring boxes miss edges to the left', E, '            x0: a.x0.min(b.x0),', '            x0: a.x0,', TESTS),
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
