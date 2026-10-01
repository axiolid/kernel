"""Mutation probe for the halfedge surface mesh (#140, ledger row D1).

Each mutant weakens one decision in
`crates/representations/discrete/mesh/src/halfedge.rs` or its `build`,
`edit` submodules -- a non-manifold refusal, the link condition, a
degenerate-collapse guard, a stored-halfedge repair after an edit, the
rotation direction -- and `crates/representations/discrete/mesh/tests/halfedge.rs`
must turn a test red.

Deliberately not listed: mutants of `validate` (`halfedge/check.rs`). It is
the oracle the tests use, and every mesh the public API can produce already
satisfies it, so weakening one of its rules changes no outcome. The random
edit tests do not rely on it alone: after every step they also convert to a
`TriMesh` and rebuild, which refuses non-manifold results independently.
"""
import pathlib, subprocess, sys

ROOT = pathlib.Path(__file__).resolve().parents[1]
M = "crates/representations/discrete/mesh/src/halfedge.rs"
B = "crates/representations/discrete/mesh/src/halfedge/build.rs"
E = "crates/representations/discrete/mesh/src/halfedge/edit.rs"
T = ["-p", "axiolid-mesh", "--test", "halfedge"]

MUTANTS = [
    (
        "an edge shared by three faces accepted",
        B,
        "        if group.len() > 2 {",
        "        if group.len() > 3 {",
    ),
    (
        "inconsistent winding accepted",
        B,
        "        if group.len() == 2 && group[1].from == first.from {",
        "        if false {",
    ),
    (
        "a vertex on two holes (bowtie) accepted",
        B,
        "            if boundary_out[source] != NONE {",
        "            if false {",
    ),
    (
        "a vertex joining two closed fans accepted",
        B,
        "        if orbit != count {",
        "        if false {",
    ),
    (
        "rotation turns clockwise",
        M,
        "        self.opposite(self.prev(h))",
        "        self.next(self.opposite(h))",
    ),
    (
        "triangle export starts at the wrong corner",
        M,
        "            for corner in [self.source(h), self.target(h), self.target(self.next(h))] {",
        "            for corner in [self.target(h), self.target(self.next(h)), self.source(h)] {",
    ),
    (
        "flip with coinciding opposite corners not refused",
        E,
        "        if c == d {\n            return Err(HalfedgeEditError::WouldDegenerate { edge: e });",
        "        if false {\n            return Err(HalfedgeEditError::WouldDegenerate { edge: e });",
    ),
    (
        "flip onto an existing edge not refused",
        E,
        "        if self.find_halfedge(VertexId(c), VertexId(d)).is_some() {",
        "        if false {",
    ),
    (
        "flip leaves a stale stored halfedge on the source",
        E,
        "        if self.vertex_out[a as usize] == h {\n            self.vertex_out[a as usize] = on;\n        }",
        "",
    ),
    (
        "flip leaves a stale stored halfedge on the target",
        E,
        "        if self.vertex_out[b as usize] == o {\n            self.vertex_out[b as usize] = hn;\n        }",
        "",
    ),
    (
        "edge split ignores the hole when storing the new vertex's halfedge",
        E,
        "        self.vertex_out[m as usize] = if fo == NONE { o } else { g };",
        "        self.vertex_out[m as usize] = g;",
    ),
    (
        "edge split leaves a stale stored halfedge on the target",
        E,
        "        if self.vertex_out[b as usize] == o {\n            self.vertex_out[b as usize] = g ^ 1;\n        }",
        "",
    ),
    (
        "face split leaves the centre without a halfedge",
        E,
        "        self.vertex_out[center as usize] = spokes[0];",
        "",
    ),
    (
        "collapse next to a polygon not refused",
        E,
        "            self.require_triangle(f)?;",
        "",
    ),
    (
        "collapse of a pillow not refused",
        E,
        "        if c.is_some() && c == d {",
        "        if false {",
    ),
    (
        "link condition ignores common neighbours",
        E,
        "Some(n) != c && Some(n) != d)",
        "Some(n) != c && Some(n) != d && n == a)",
    ),
    (
        "link condition forgets one opposite corner",
        E,
        "Some(n) != c && Some(n) != d)",
        "Some(n) != c)",
    ),
    (
        "link condition ignores the boundary as a common neighbour",
        E,
        "        if self.is_boundary_vertex(a) && self.is_boundary_vertex(b) && !self.is_boundary_edge(edge)",
        "        if false",
    ),
    (
        "collapse of a tetrahedron not refused",
        E,
        "                    .all(|&v| !self.is_boundary_vertex(v) && self.degree(v) == 3)",
        "                    .all(|_| false)",
    ),
    (
        "collapse of a lone triangle not refused",
        E,
        "                if self.loop_halfedges(HalfedgeId(hole)).count() == 3 {",
        "                if self.loop_halfedges(HalfedgeId(hole)).count() == 2 {",
    ),
    (
        "collapse leaves the kept vertex storing the removed halfedge",
        E,
        "        if self.vertex_out[b.index()] == o {",
        "        if false {",
    ),
    (
        "collapse skips the boundary repair on the kept vertex",
        E,
        "            self.adjust_outgoing(v.0);",
        "",
    ),
    (
        "merged edge leaves the outer face storing the removed halfedge",
        E,
        "        if outer != NONE && self.face_halfedge[outer as usize] == o1 {",
        "        if false {",
    ),
    (
        "merged edge leaves the far corner storing the removed halfedge",
        E,
        "        if self.vertex_out[v0 as usize] == h1 {",
        "        if false {",
    ),
    (
        "merged edge leaves the kept vertex storing the removed halfedge",
        E,
        "        if self.vertex_out[v1 as usize] == o1 {",
        "        if false {",
    ),
    (
        "diagonal back along a side not refused",
        E,
        "        if a == b || self.next(a) == b || self.next(b) == a {",
        "        if a == b || self.next(a) == b {",
    ),
    (
        "filling from a face halfedge not refused",
        E,
        "        if !self.is_boundary_halfedge(h) {",
        "        if false {",
    ),
]


def run(target):
    return subprocess.run(
        ["cargo", "test", "-q", *target],
        cwd=ROOT, capture_output=True, text=True, timeout=1200,
    ).returncode


survivors = []
for name, rel, old, new in MUTANTS:
    path = ROOT / rel
    original = path.read_text()
    assert original.count(old) == 1, f"anchor for '{name}' not unique/found"
    path.write_text(original.replace(old, new))
    try:
        try:
            code = run(T)
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
