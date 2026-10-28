"""Mutation probe for splitting constraints at the vertices on them (#360)
and for the domain triangulation (#353) of `axiolid-triangulate`.

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

CDT = ["-p", "axiolid-triangulate", "--test", "cdt"]
DUP = ["-p", "axiolid-triangulate", "--test", "duplicates"]


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
]


probe_lib.main(MUTANTS, timeout=TIMEOUT, reject_uncompiled=True)
