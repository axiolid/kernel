"""Mutation probe for the bounded free faces of an obstacle set (#490, ADR 0109,
`crates/algorithms/planar/overlay/src/free_faces/` and the `faces` and
distance predicates of `crates/algorithms/planar/exact-arrangement`).

Each mutant weakens one decision -- which pairs bridge a gap and which a
barrier shadows, how a face is read as obstacle or free, which faces are
closed gaps, what the rings keep and list, the order of the output, how a
point is located, and the numbering of faces and distances underneath -- and
must turn a test red.

Equivalent mutants, deliberately not listed:
- `dot <= 0` against `dot < 0` in `nearest_on_segment`: the foot of a
  perpendicular that lands on the end is the end itself;
- leaving out the dedup of vertex pairs found from both sides: two identical
  bridge segments are one piece, and equal rounded ends merge into one record;
- the vertex a closed-gap scan starts at (`count / 2`): the scan stops at the
  first vertex that is not near, and the verdict does not depend on the order;
- the box inflation by more than the gap: boxes only skip work;
- starting rings at `canonical_start`: `snap_round` numbers its vertices by `point_order`, so a
  traced ring already starts at its smallest point, at its lowest rotation where it visits that
  point twice; the call stays so that the contract does not rest on the numbering;
- skipping the pieces a vertex is on in `bridges::find`: such a pair is at distance
  zero, and the piece itself, which covers the vertex, shadows it.
"""
import probe_lib

OVERLAY = "crates/algorithms/planar/overlay/src/free_faces"
B = OVERLAY + "/bridges.rs"
C = OVERLAY + "/classify.rs"
R = OVERLAY + "/rings.rs"
L = OVERLAY + "/locate.rs"
M = OVERLAY + ".rs"
ARR = "crates/algorithms/planar/exact-arrangement/src"
F = ARR + "/faces.rs"
D = ARR + "/distance.rs"

OVERLAY_TESTS = ["-p", "axiolid-overlay", "--test", "free_faces"]
FACES_TESTS = ["-p", "axiolid-exact-arrangement", "--test", "faces"]
DISTANCE_TESTS = ["-p", "axiolid-exact-arrangement", "--test", "distance"]
BRIDGE_UNITS = ["-p", "axiolid-overlay", "--lib", "free_faces::"]
TESTS = [OVERLAY_TESTS, BRIDGE_UNITS]

MUTANTS = [
    ("a distance equal to the gap is bridged", B,
     "if closer != Ordering::Less {",
     "if closer == Ordering::Greater {", TESTS),
    ("nothing shadows a bridge", B,
     "if meets_open_segment(points, (here, there), (a, b))? {",
     "if false && meets_open_segment(points, (here, there), (a, b))? {", TESTS),
    ("the open-segment test never sees a collinear overlap", B,
     "    Ok(o1 == Sign::Zero\n        && o2 == Sign::Zero",
     "    Ok(false && o1 == Sign::Zero\n        && o2 == Sign::Zero", [BRIDGE_UNITS]),
    ("candidate boxes are not grown by the gap", B,
     "BoxTree::new(boxes(gap))",
     "BoxTree::new(boxes(0.0))", TESTS),
    ("only the first carrier of a piece is a source", B,
     "sorted(piece.segments.iter().map(|&segment| sources[segment]))",
     "sorted(piece.segments.iter().take(1).map(|&segment| sources[segment]))", TESTS),
    ("an end of a piece is never the nearest point", D,
     "        if ask(Question::Beyond)? != Sign::Negative {",
     "        if false && ask(Question::Beyond)? != Sign::Negative {", [DISTANCE_TESTS] + TESTS),
    ("the start of a piece is never the nearest point", D,
     "        if ask(Question::Behind)? != Sign::Positive {",
     "        if false && ask(Question::Behind)? != Sign::Positive {", [DISTANCE_TESTS] + TESTS),
    ("the foot's weight is dropped", D,
     "            a[2].mul(&scale),",
     "            a[2].clone(),", [DISTANCE_TESTS] + TESTS),
    ("a bridge inside an obstacle is kept", C,
     "bridge_alive.push(!obstacle[faces.left(piece).0]);",
     "bridge_alive.push(true);", TESTS),
    ("the holes of a polygon are ignored", C,
     "            && self.holes_in[polygon] == 0",
     "            && true", TESTS),
    ("no free face is a closed gap", C,
     "let gap = gap > 0.0 && fits_in_strip(",
     "let gap = false && fits_in_strip(", TESTS),
    ("a piece with one face on both sides bounds it", C,
     "        if left == right {\n            continue;\n        }",
     "        if false && left == right {\n            continue;\n        }", TESTS),
    ("the kept-vertex rule ignores the sources", C,
     "keep[vertex] = labels[first] != labels[second]\n                || !straight(",
     "keep[vertex] = false\n                || !straight(", TESTS),
    ("faces are not sorted", R,
     "        list.sort_by(|a, b| sequence_order(&a.outer.points, &b.outer.points));",
     "        let _ = list;", TESTS),
    ("an obstacle ring on one line is accepted", M,
     "            .any(|ring| exact_overlay::on_one_line(&ring.points))",
     "            .any(|ring| false && exact_overlay::on_one_line(&ring.points))", TESTS),
    ("locate skips the boundary pass", L,
     "                if let Some(edge) = on_edge(&found.points, point) {",
     "                if let Some(edge) = on_edge(&found.points, point).filter(|_| false) {", TESTS),
    ("a point on the edge of a closed gap is in the obstacle", L,
     "                .any(|ring| on_edge(&ring.points, point).is_some())\n                || inside(gap, point)",
     "                .any(|ring| false && on_edge(&ring.points, point).is_some())\n                || inside(gap, point)", TESTS),
    ("closed gaps are never located", L,
     "self.closed_gaps.iter().position(in_gap)",
     "self.closed_gaps.iter().position(|_| false)", TESTS),
    ("faces: the wide wedge is the one turning counter-clockwise", F,
     "points.orient(here, far(spoke), far(after))? == Sign::Negative",
     "points.orient(here, far(spoke), far(after))? == Sign::Positive", [FACES_TESTS] + TESTS),
    ("faces: a component goes to the outermost face that holds it", F,
     "                        .is_gt(),",
     "                        .is_lt(),", [FACES_TESTS] + TESTS),
    ("faces: components are not nested", F,
     "    if roots.len() > 1 {",
     "    if false && roots.len() > 1 {", [FACES_TESTS] + TESTS),
]

probe_lib.main(MUTANTS, timeout=900)
