"""Mutation probe for the exact planar-faced boolean (#199, #200).

Each mutant reintroduces a way `boolean_polyhedra_exact` drifted, left a
hole, or misclassified -- a split point rounded twice, a thin fragment
classified at a rounded centroid, a coplanar pair decided the wrong way --
and must turn a test red.

Equivalent mutants, deliberately not listed: dropping the shared-coordinate
case of `centroid_of` (thin axis-aligned fragments then take the exact probe
path: slower, same answer), a shorter or longer seed walk in
`nearest_ratio` (the bisection it falls back to returns the same double).
"""
import pathlib, subprocess, sys

ROOT = pathlib.Path(__file__).resolve().parents[1]
P = "crates/algorithms/construction/construct/src/polyhedron.rs"
X = "crates/algorithms/construction/construct/src/polyhedron_exact.rs"
TESTS = [
    ["-p", "axiolid-construct", "--lib", "polyhedron"],
    ["-p", "axiolid-construct", "--test", "grid_subtraction"],
    ["-p", "axiolid-construct", "--test", "contact_matrix"],
    ["-p", "axiolid-construct", "--test", "boolean_polyhedra"],
]

MUTANTS = [
    ("split points from the twice-rounded f64 formula", P,
     "    exact::plane_crossing(plane, a, b, a + (b - a) * t)",
     "    Some(a + (b - a) * t)"),
    ("a tie met walking down leaves the even double", X,
     "                Sign::Zero if !even(r) => r = down,",
     "                Sign::Zero => r = down,"),
    ("a bisected tie rounds up whatever the parity", X,
     "        Sign::Zero if !even(floor) => next,",
     "        Sign::Zero => next,"),
    ("bisection keeps the wrong half", X,
     "        if versus(unkey(mid))? == Sign::Negative {",
     "        if versus(unkey(mid))? != Sign::Negative {"),
    ("every fragment classified at its f64 centroid", P,
     "        if side_of_face(fragment, centroid) == Some(Sign::Zero)\n            && ring_position(fragment, centroid) == Some(RingPosition::Inside)\n        {",
     "        if true {"),
    ("a centroid on its own edge accepted as interior", P,
     "            && ring_position(fragment, centroid) == Some(RingPosition::Inside)",
     "            && ring_position(fragment, centroid) != Some(RingPosition::Outside)"),
    ("a centroid off the fragment's plane accepted", P,
     "        if side_of_face(fragment, centroid) == Some(Sign::Zero)\n",
     "        if side_of_face(fragment, centroid).is_some()\n"),
    ("exact ring parity counted on the wrong side", X,
     "            if (upward && s == Sign::Negative) || (!upward && s == Sign::Positive) {",
     "            if (upward && s == Sign::Negative) || (!upward && s == Sign::Negative) {"),
    ("an exact point on an edge reported inside", X,
     "            return Some(RingPosition::OnEdge);",
     "            return Some(RingPosition::Inside);"),
    ("an exact interior point accepted on an edge", X,
     "        if ring_position(ring, &point)? == RingPosition::Inside {",
     "        if ring_position(ring, &point)? != RingPosition::Outside {"),
    ("coplanar normals compared the wrong way", X,
     "    Some(sign(&dot(&normal(first)?, &normal(second)?)) == Sign::Positive)",
     "    Some(sign(&dot(&normal(first)?, &normal(second)?)) != Sign::Positive)"),
]

def run(target):
    return subprocess.run(
        ["cargo", "test", "-q", *target],
        cwd=ROOT, capture_output=True, text=True, timeout=1800,
    ).returncode

survivors = []
for name, rel, old, new in MUTANTS:
    path = ROOT / rel
    original = path.read_text()
    assert original.count(old) == 1, f"anchor for '{name}' not unique/found"
    path.write_text(original.replace(old, new))
    try:
        try:
            code = max(run(t) for t in TESTS)
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
