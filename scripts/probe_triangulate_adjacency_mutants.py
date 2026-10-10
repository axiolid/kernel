"""Mutation probe for the locally maintained triangulation adjacency (#268).

`axiolid-triangulate` keeps its twin pointers current by re-linking only the
halfedges a split or flip rewrites and their outside neighbours. Each
mutant drops or misdirects one of those links -- in the one-into-three
split, the two-into-four edge split, the Lawson flip after insertion and
the flip of constraint recovery and restoration -- or weakens the guard
that holds the local result to a full recomputation: the debug check
itself, its comparison, its skip on a repeated directed edge, its size
limit, the check after each kind of update, and the claim that no full
recomputation maintains the mesh. The grid hint that starts point
location is probed too: ignored, trusted for a point on an edge, or
admitting the edge as inside. Each must turn the crate's unit tests
(`src/adjacency.rs`, `src/hint.rs`) or `tests/cdt.rs` red.

Not probed: the order of independent `link` calls within one update (any
order of disjoint writes gives the same array) and the size limit's value
outside the unit tests (it only bounds debug-build cost downstream). A
stale per-vertex triangle for the hint would only lengthen a walk, so a
debug assertion holds it exact, and dropping its update after a flip is
probed.

Runs on the shared runner (`scripts/probe_lib.py`): each mutant in its own
process group under a timeout that kills the group, the file restored
whatever happens, and a missing or ambiguous anchor aborts the probe. Run
it in a worktree, never in a checkout others are building in.
"""
import probe_lib

TIMEOUT = 900
A = "crates/algorithms/planar/triangulate/src/adjacency.rs"
B = "crates/algorithms/planar/triangulate/src/build.rs"
R = "crates/algorithms/planar/triangulate/src/recover.rs"

UNIT = ["-p", "axiolid-triangulate", "--lib"]
CDT = ["-p", "axiolid-triangulate", "--test", "cdt"]
BOTH = [UNIT, CDT]


def drop(line):
    """Delete one whole line, kept unique by its indentation and comment."""
    return line, ""


H = "crates/algorithms/planar/triangulate/src/hint.rs"


MUTANTS = [
    # One into three.
    ("split: the old neighbour across b -> c is not re-pointed", B,
     *drop("        link(h, 3 * t1, across_bc);\n"), BOTH),
    ("split: the old neighbour across c -> a is not re-pointed", B,
     *drop("        link(h, 3 * t2, across_ca);\n"), BOTH),
    ("split: b -> v is paired with the wrong spoke", B,
     "        link(h, 3 * t + 1, (3 * t1 + 2) as u32); // b -> v | v -> b",
     "        link(h, 3 * t + 1, (3 * t2 + 2) as u32); // b -> v | v -> b", BOTH),
    ("split: c -> v is left unpaired", B,
     *drop("        link(h, 3 * t1 + 1, (3 * t2 + 2) as u32); // c -> v | v -> c\n"), BOTH),
    ("split: a -> v is left unpaired", B,
     *drop("        link(h, 3 * t2 + 1, (3 * t + 2) as u32); // a -> v | v -> a\n"), BOTH),
    # Two into four.
    ("edge split: x -> u keeps a stale twin", B,
     *drop("        link(h, 3 * t, across_xu);\n"), BOTH),
    ("edge split: y -> w keeps a stale twin", B,
     *drop("        link(h, 3 * ot, across_yw);\n"), BOTH),
    ("edge split: w -> x is not re-pointed", B,
     *drop("        link(h, 3 * t1, across_wx);\n"), BOTH),
    ("edge split: u -> y is not re-pointed", B,
     *drop("        link(h, 3 * t2, across_uy);\n"), BOTH),
    ("edge split: u -> v is paired with the wrong spoke", B,
     "        link(h, 3 * t + 1, (3 * t2 + 2) as u32); // u -> v | v -> u",
     "        link(h, 3 * t + 1, (3 * t1 + 2) as u32); // u -> v | v -> u", BOTH),
    ("edge split: w -> v is left unpaired", B,
     *drop("        link(h, 3 * ot + 1, (3 * t1 + 2) as u32); // w -> v | v -> w\n"), BOTH),
    # The flip after insertion.
    ("flip: p0 -> other keeps the old diagonal as its twin", B,
     *drop("            link(h, 3 * t + ti, across_p0_other);\n"), BOTH),
    ("flip: p1 -> apex keeps the old diagonal as its twin", B,
     *drop("            link(h, 3 * ot + oi, across_p1_apex);\n"), BOTH),
    ("flip: the new diagonal is left unpaired", B,
     *drop("            link(h, 3 * t + (ti + 1) % 3, (3 * ot + (oi + 1) % 3) as u32);\n"), BOTH),
    ("flip: p1 -> apex takes the wrong outside twin", B,
     "            let across_p1_apex = self.halfedges[3 * t + (ti + 1) % 3];",
     "            let across_p1_apex = self.halfedges[3 * ot + (oi + 1) % 3];", BOTH),
    # Recovery's flip.
    ("recovery flip: the new diagonal is left unpaired", R,
     *drop("    link(h, 3 * t, (3 * ot) as u32); // apex -> other | other -> apex\n"), BOTH),
    ("recovery flip: other -> b is not re-pointed", R,
     *drop("    link(h, 3 * t + 1, across_other_b);\n"), BOTH),
    ("recovery flip: b -> apex is not re-pointed", R,
     *drop("    link(h, 3 * t + 2, across_b_apex);\n"), BOTH),
    ("recovery flip: apex -> a is not re-pointed", R,
     *drop("    link(h, 3 * ot + 1, across_apex_a);\n"), BOTH),
    ("recovery flip: a -> other is not re-pointed", R,
     *drop("    link(h, 3 * ot + 2, across_a_other);\n"), BOTH),
    # The guard.
    ("the debug check never runs", A,
     "    if !cfg!(debug_assertions) || too_large {",
     "    if true || too_large {", [UNIT]),
    ("the comparison reports every array equal", A,
     "        Some(e) => Comparison::Differs(e),",
     "        Some(_) => Comparison::Equal,", [UNIT]),
    ("every triangle list counts as ambiguous", A,
     "    if !manifold {\n        return Comparison::Ambiguous;",
     "    if manifold || !manifold {\n        return Comparison::Ambiguous;", [UNIT]),
    ("a difference is skipped, not reported", A,
     '        Comparison::Differs(e) => panic!(',
     '        Comparison::Differs(e) => println!(', [UNIT]),
    ("the unit tests check only small meshes", A,
     "if cfg!(test) { None } else { Some(1024) }",
     "if cfg!(test) { Some(64) } else { Some(1024) }", [UNIT]),
    ("no check after a split", B,
     "        self.corner[v as usize] = t as u32;\n"
     "        note_update();\n        debug_check(&self.triangles, &self.halfedges);\n",
     "        self.corner[v as usize] = t as u32;\n"
     "        note_update();\n", [UNIT]),
    ("no check after an edge split", B,
     "            self.corner[vertex as usize] = tri as u32;\n        }\n"
     "        note_update();\n        debug_check(&self.triangles, &self.halfedges);\n",
     "            self.corner[vertex as usize] = tri as u32;\n        }\n"
     "        note_update();\n", [UNIT]),
    ("no check after a flip", B,
     "            note_update();\n            debug_check(&self.triangles, &self.halfedges);\n",
     "            note_update();\n", [UNIT]),
    ("no check after a recovery flip", R,
     "    note_update();\n    debug_check(&tri.triangles, &tri.halfedges);\n",
     "    note_update();\n", [UNIT]),
    ("a split recomputes adjacency in full again", B,
     "        link(h, 3 * t2 + 1, (3 * t + 2) as u32); // a -> v | v -> a\n",
     "        link(h, 3 * t2 + 1, (3 * t + 2) as u32); // a -> v | v -> a\n"
     "        *h = crate::adjacency::twins(&self.triangles);\n", [UNIT]),
    ("a recovery flip recomputes adjacency in full again", R,
     "    link(h, 3 * ot + 2, across_a_other);\n",
     "    link(h, 3 * ot + 2, across_a_other);\n"
     "    *h = crate::adjacency::twins(&tri.triangles);\n", [UNIT]),
    # Where the walk starts.
    ("the walk ignores the hint", B,
     "            t\n        });",
     "            last\n        });", [UNIT]),
    ("a hinted walk that ends on an edge is trusted", B,
     "                if self.strictly_inside(t, p) {",
     "                if true {", [UNIT]),
    ("a point on an edge counts as strictly inside", B,
     "decided(orient2d(u, w, p)) == Sign::Positive",
     "decided(orient2d(u, w, p)) != Sign::Negative", [UNIT]),
    ("the hint grid records no vertex", H,
     "        self.cells[cell] = v;",
     "        self.cells[cell] = v.max(EMPTY);", [UNIT]),
    ("a split leaves c's triangle stale", B,
     *drop("        self.corner[c as usize] = t1 as u32;\n"), [UNIT]),
    ("a flip leaves p1's triangle stale", B,
     *drop("            self.corner[p1 as usize] = ot as u32;\n"), [UNIT]),
]


probe_lib.main(MUTANTS, timeout=TIMEOUT, reject_uncompiled=True)
