"""Mutation probe for splitting constraints at the vertices on them (#360)
and for the domain triangulation (#353) of `axiolid-triangulate`.

The domain (#353) is probed in the rule that decides it: a ring's
orientation (kept as given, read at the first vertex, not the lowest, or the
role of a hole ignored), the pieces that carry winding (a breakline that
does, a split ring edge looked up whole), the flood fill (the sign of the
change across a piece, the seed's hull edge, a winding outside 0 or 1
accepted, an empty domain returned), the ring checks (a repeat read on
input indices or skipped, the fold-back check skipped), the outside points
(a merged point listed, none reported), the constraints outside the domain
kept, and refinement that rebuilds the whole hull. Each must turn
`tests/domain.rs` or the unit tests in `src/domain.rs` red.

Not probed: the debug-only assertion that a revisited triangle agrees with
its first winding, and that a hull edge agrees with the fill. The winding
difference between two triangles is the signed number of ring boundaries
between them, whatever path joins them, so release behaviour cannot depend
on the assertion.

A constraint that passes exactly through input vertices survives as the
chain of pieces between them. Each mutant here breaks one part of that:
the split itself, the closed bounding box that finds the vertices on a
segment, the order of the pieces along the constraint, and the rule that a
merged duplicate is no vertex. Each must turn `tests/cdt.rs` or
`tests/duplicates.rs` red.

Runs on the shared runner (`scripts/probe_lib.py`): each mutant in its own
process group under a timeout that kills the group, the file restored
whatever happens, and a missing or ambiguous anchor aborts the probe. Run
it in a worktree, never in a checkout others are building in.
"""
import probe_lib

TIMEOUT = 900
R = "crates/algorithms/planar/triangulate/src/recover.rs"
D = "crates/algorithms/planar/triangulate/src/domain.rs"
F = "crates/algorithms/planar/triangulate/src/refine.rs"

CDT = ["-p", "axiolid-triangulate", "--test", "cdt"]
DUP = ["-p", "axiolid-triangulate", "--test", "duplicates"]
LIB = ["-p", "axiolid-triangulate", "--lib"]
DOM = ["-p", "axiolid-triangulate", "--test", "domain"]


MUTANTS = [
    # Splitting constraints at the vertices on them (#360).
    ("a constraint through a vertex is not split", R,
     "        chain.extend(vertices_on(tri, &is_vertex, c.a, c.b));",
     "        chain.extend(Vec::<u32>::new());", [CDT]),
    ("the bounding-box filter is strict", R,
     "                && lo_x <= q.x\n",
     "                && lo_x < q.x\n", [CDT]),
    ("pieces are not ordered along the constraint", R,
     "    if (pa.x, pa.y) > (pb.x, pb.y) {",
     "    if false {", [CDT]),
    ("a merged duplicate counts as a vertex on the constraint", R,
     "            is_vertex[v as usize]\n                && v != a",
     "            v != a", [DUP]),
    # The domain (#353).
    ("a ring keeps the orientation it was given", D,
     "        let sign = if is_ccw(points, &seq) == outer { 1 } else { -1 };",
     "        let sign = if outer { 1 } else { -1 };", [DOM]),
    ("orientation is read at the first vertex, not the lowest", D,
     "    let lowest = (0..m)\n        .min_by(|&i, &j| lex(at(i), at(j)))\n        .unwrap_or_default();",
     "    let lowest = 0;", [DOM]),
    ("holes count like outer rings", D,
     "        let outer = ring.kind == RingKind::Outer;",
     "        let outer = true;", [DOM]),
    ("a breakline carries winding", D,
     "            let w = winding[t] - net;",
     "            let w = winding[t] - if net == 0 && tri.is_constrained(u, v) { 1 } else { net };",
     [DOM]),
    ("a split ring edge is looked up whole", D,
     "                .chain(key)\n",
     "                .chain(Constraint::new(0, 0))\n", [DOM]),
    ("the winding change across a piece has the wrong sign", D,
     "            let w = winding[t] - net;",
     "            let w = winding[t] + net;", [DOM]),
    ("the seed ignores its hull edge", D,
     "    winding[seed.0] = net;",
     "    winding[seed.0] = 0;", [DOM]),
    ("a winding outside 0 or 1 is accepted", D,
     "    if (0..=1).contains(&winding) {",
     "    if true {", [DOM]),
    ("an empty domain is returned", D,
     "    if triangles.is_empty() {",
     "    if false {", [DOM]),
    ("the repeated-vertex check reads input indices", D,
     "        for (i, &k) in seq.iter().enumerate() {",
     "        for (i, &k) in ring.vertices.iter().enumerate() {", [LIB]),
    ("the repeated-vertex check is skipped", D,
     "            if stamp[k as usize] == r {",
     "            if false {", [LIB]),
    ("the fold-back check is skipped", D,
     "            if folds {",
     "            if false {", [LIB]),
    ("a merged point is listed as outside", D,
     "        .filter(|&v| !used[v as usize] && !merged[v as usize])",
     "        .filter(|&v| !used[v as usize])", [DOM]),
    ("outside points are not reported", D,
     "    tri.outside = (0..tri.points.len() as u32)",
     "    tri.outside = Vec::new();\n    let _unused: Vec<u32> = (0..tri.points.len() as u32)", [DOM]),
    ("constraints outside the domain are kept", D,
     "        .retain(|c| edges.binary_search(&(c.a, c.b)).is_ok());",
     "        .retain(|c| edges.binary_search(&(c.a, c.b)).is_ok() || true);", [DOM]),
    ("refine rebuilds the whole hull", F,
     "rebuild(&points, &constraints, tri.domain.as_deref())",
     "rebuild(&points, &constraints, None)", [DOM]),
]


probe_lib.main(MUTANTS, timeout=TIMEOUT, reject_uncompiled=True)
