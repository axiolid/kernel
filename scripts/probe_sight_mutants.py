"""Mutation probe for certified line of sight (#185).

Each mutant weakens one side of the argument -- lets a sub-cone count as
covered without being strictly inside a piece's cone or without a
separating plane, lets a witness ray through a blocker, or accepts a
non-convex solid -- and must turn a test red.

Not listed, deliberately:

- a non-strict cone test. Under the convention that a ray touching a
  blocker's edge is blocked (the witness uses it too), a sub-cone on the
  boundary of a piece's cone is hidden as well; the test is kept strict
  so that a graze reads undecided, the more cautious of two sound
  answers.
- accepting a witness ray that only touches the target triangle. The
  sample points are interior points of the triangle, so the ray crosses
  its interior unless the triangle is seen edge-on, and then the plane
  parameter is undefined and the sample is refused anyway.
"""
import pathlib, subprocess, sys

ROOT = pathlib.Path(__file__).resolve().parents[1]
S = "crates/algorithms/query/inspect/src/sight.rs"
TESTS = ["-p", "axiolid-inspect", "--test", "sight", "--lib"]

MUTANTS = [
    ('no separating plane needed', S, '    piece.faces.iter().any(|face| {', '    true || piece.faces.iter().any(|face| {', TESTS),
    ('witness ignores blockers', S, '                if walls.iter().all(|wall| !blocks(eye, through, wall, &near)) {', '                if true || walls.iter().all(|wall| !blocks(eye, through, wall, &near)) {', TESTS),
    ('grazing a blocker edge does not block', S, '    if has_pos && has_neg {\n        return false;\n    }', '    if has_pos && has_neg || s.contains(&Sign::Zero) {\n        return false;\n    }', TESTS),
    ('blockers behind the target block', S, '    param.sign() != Sign::Negative && param.compare(near) != Sign::Positive', '    param.sign() != Sign::Negative', TESTS),
    ('convex solids unchecked for convexity', S, '                s if s != side => return None,', '                s if s != side => {}', TESTS),
    ('solid pieces dropped', S, '    if let Some(solid) = convex_solid(eye, blocker, own, mesh) {', '    if let Some(solid) = convex_solid(eye, blocker, own, mesh).filter(|_| false) {', TESTS),
    ('walls not merged', S, '            if let Some(quad) = convex_quad(a, b) {', '            if let Some(quad) = convex_quad(a, b).filter(|_| false) {', TESTS),
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
