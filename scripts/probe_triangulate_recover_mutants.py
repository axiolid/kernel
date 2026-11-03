"""Mutation probe for constraint recovery through the vertex-to-triangle
index and the segment walk (#355).

`axiolid-triangulate` no longer scans every triangle per missing
constraint. It finds an edge through an index that names one triangle per
vertex (`flip` keeps it current), rotates around the vertex through the
twins (`outgoing`), and finds the edges a constraint crosses and the
vertices on it by walking the segment with the exact orientation
predicate (`walk`). Each mutant misses one corner of that: a stale index
entry after a flip, a fan that stops at the boundary or never leaves one
side, an edge on the hull that is not found, the wrong halfedge of an edge
(which changes the flip, so the output), a vertex on the segment that is
not recorded, a walk that does not stop at the end, a wedge test that
drops the rays, a crossing edge listed once or out of order, or another
constraint offered for flipping. Each must turn the crate's unit tests
(`src/recover.rs`), which hold the walk and the index to a scan of every
triangle and the output to the hash of the scanning implementation, or
`tests/cdt.rs` red.

Not probed: which triangle the index names for a vertex (any triangle of
its fan serves, so a different choice is an equivalent mutant); the refusal of a walk that meets a vertex on the segment
(`find_crossing` returns `None`): the pieces are split at those vertices
first, so no input reaches it; and the walk budgets, which only bound a
corrupted mesh.

Runs on the shared runner (`scripts/probe_lib.py`): each mutant in its own
process group under a timeout that kills the group, the file restored
whatever happens, and a missing or ambiguous anchor aborts the probe. Run
it in a worktree, never in a checkout others are building in.
"""
import probe_lib

TIMEOUT = 900
R = "crates/algorithms/planar/triangulate/src/recover.rs"
B = "crates/algorithms/planar/triangulate/src/build.rs"

UNIT = ["-p", "axiolid-triangulate", "--lib"]
CDT = ["-p", "axiolid-triangulate", "--test", "cdt"]
BOTH = [UNIT, CDT]


def drop(line):
    """Delete one whole line, kept unique by its indentation."""
    return line, ""


MUTANTS = [
    # The index.
    ("a flip leaves a's triangle stale", R,
     *drop("    incident[a as usize] = ot as u32;\n"), [UNIT]),
    ("a flip leaves the other three corners stale", R,
     "\n        incident[v as usize] = t as u32;\n", "\n        let _ = v;\n", [UNIT]),
    # The fan around a vertex.
    ("a fan is walked one way only", R,
     "    if !closed {\n        h = start;",
     "    if false {\n        h = start;", BOTH),
    ("a fan rotates by the wrong halfedge", R,
     "        let twin = tri.halfedges[prev_halfedge(h)];\n        if twin == NO_HALFEDGE {\n            break;\n        }\n        if twin as usize == start {",
     "        let twin = tri.halfedges[h];\n        if twin == NO_HALFEDGE {\n            break;\n        }\n        if twin as usize == start {",
     [UNIT]),
    ("a fan that never ends is walked for ever, first way", R,
     "        h = twin as usize;\n        if fan.len() > bound {",
     "        h = twin as usize;\n        if false && fan.len() > bound {", [UNIT]),
    ("a fan that never ends is walked for ever, second way", R,
     "            fan.push(h);\n            if fan.len() > bound {",
     "            fan.push(h);\n            if false && fan.len() > bound {", [UNIT]),
    # Self-loops and non-finite points.
    ("a constraint from a vertex to itself is looked up as an edge", R,
     "    for c in tri.constraints.clone() {\n        if c.a == c.b || has_edge(",
     "    for c in tri.constraints.clone() {\n        if has_edge(", [UNIT]),
    ("a non-finite point is not refused", B,
     "        .position(|p| !(p.x.is_finite() && p.y.is_finite()))",
     "        .position(|_| false)", [UNIT]),
    ("only a non-finite x is refused", B,
     "        .position(|p| !(p.x.is_finite() && p.y.is_finite()))",
     "        .position(|p| !p.x.is_finite())", [UNIT]),
    # The edge lookup.
    ("a boundary edge into the vertex is not found", R,
     "        if tri.triangles[into] == v {",
     "        if tri.triangles[into] == u32::MAX {", BOTH),
    ("an edge out of the vertex is not found", R,
     "        if tri.triangles[next_halfedge(h)] == v {",
     "        if tri.triangles[next_halfedge(h)] == u32::MAX {", BOTH),
    ("an edge answers with the higher halfedge", R,
     "best = Some(best.map_or(h, |b| b.min(h)))",
     "best = Some(best.map_or(h, |b| b.max(h)))", [UNIT]),
    # The walk.
    ("a vertex on the segment is not recorded after a crossing", R,
     *drop("                    out.through.push(apex);\n"), [UNIT]),
    ("a vertex on the segment is not recorded at the start", R,
     *drop("                out.through.push(v);\n"), [UNIT]),
    ("the walk does not stop at the end after a crossing", R,
     "            if apex == b {\n                cur = b;\n                break;\n            }\n",
     "", [UNIT]),
    ("the wedge test drops the rays", R,
     "(sp >= 0 && sq <= 0).then_some(",
     "(sp > 0 && sq < 0).then_some(", BOTH),
    ("the wedge test takes the opposite direction too", R,
     "(sp >= 0 && sq <= 0).then_some(",
     "(sp >= 0 || sq <= 0).then_some(", BOTH),
    ("a point to the right of the line counts as on it", R,
     "    } else if turns_left(p(b), p(a), p(c)) {\n        -1\n",
     "    } else if false {\n        -1\n", BOTH),
    ("the walk leaves by the wrong edge when the apex is left", R,
     "                    e = next_halfedge(twin);\n",
     "                    e = prev_halfedge(twin);\n", BOTH),
    # The crossing list.
    ("a crossing edge is listed once", R,
     *drop("        halfedges.push(tri.halfedges[e] as usize);\n"), [UNIT]),
    ("the crossing edges are not in the order of the triangle array", R,
     *drop("    halfedges.sort_unstable();\n"), [UNIT]),
    ("another constraint is offered for flipping", R,
     "        .filter(|&h| !tri.is_constrained(tri.triangles[h], tri.triangles[next_halfedge(h)]))\n",
     "", [UNIT]),
]


probe_lib.main(MUTANTS, timeout=TIMEOUT, reject_uncompiled=True)
