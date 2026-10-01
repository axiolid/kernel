//! Full structural validation, independent of how the mesh was produced.

use core::fmt;

use super::{HalfedgeId, HalfedgeMesh, VertexId, NONE};

/// The first structural invariant a [`HalfedgeMesh`] was found to break.
#[non_exhaustive]
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct HalfedgeInvariantError {
    /// Which rule broke, in words.
    pub rule: &'static str,
    /// Kind of element the rule was checked on.
    pub element: &'static str,
    /// Index of that element.
    pub index: u32,
}

impl fmt::Display for HalfedgeInvariantError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "{} {}: {}", self.element, self.index, self.rule)
    }
}

impl std::error::Error for HalfedgeInvariantError {}

fn fail(rule: &'static str, element: &'static str, index: usize) -> HalfedgeInvariantError {
    HalfedgeInvariantError {
        rule,
        element,
        index: index as u32,
    }
}

pub(super) fn validate(mesh: &HalfedgeMesh) -> Result<(), HalfedgeInvariantError> {
    let live_count = |flags: &[bool]| flags.iter().filter(|&&alive| alive).count();
    if live_count(&mesh.vertex_live) != mesh.live_vertices
        || live_count(&mesh.edge_live) != mesh.live_edges
        || live_count(&mesh.face_live) != mesh.live_faces
        || mesh.vertex_live.len() != mesh.positions.len()
        || mesh.vertex_out.len() != mesh.positions.len()
        || mesh.links.len() != mesh.edge_live.len() * 2
        || mesh.face_halfedge.len() != mesh.face_live.len()
    {
        return Err(fail("element counts disagree with storage", "mesh", 0));
    }
    let halfedge_live = |h: u32| h != NONE && mesh.edge_live.get((h / 2) as usize) == Some(&true);
    let vertex_live = |v: u32| v != NONE && mesh.vertex_live.get(v as usize) == Some(&true);
    let face_live = |f: u32| f != NONE && mesh.face_live.get(f as usize) == Some(&true);

    let mut out_count = vec![0u32; mesh.positions.len()];
    let mut boundary_out = vec![0u32; mesh.positions.len()];
    for e in mesh.edges() {
        let [a, b] = [e.0 * 2, e.0 * 2 + 1];
        if mesh.links[a as usize].face == NONE && mesh.links[b as usize].face == NONE {
            return Err(fail("both sides are boundary", "edge", e.index()));
        }
        for h in [a, b] {
            let link = mesh.links[h as usize];
            if !halfedge_live(link.next) || !halfedge_live(link.prev) {
                return Err(fail("links to a removed halfedge", "halfedge", h as usize));
            }
            if mesh.links[link.next as usize].prev != h || mesh.links[link.prev as usize].next != h
            {
                return Err(fail(
                    "next and prev are not inverse",
                    "halfedge",
                    h as usize,
                ));
            }
            if !vertex_live(link.target) {
                return Err(fail("targets a removed vertex", "halfedge", h as usize));
            }
            if link.face != NONE && !face_live(link.face) {
                return Err(fail("lies on a removed face", "halfedge", h as usize));
            }
            if mesh.links[link.next as usize].face != link.face {
                return Err(fail("next lies on another face", "halfedge", h as usize));
            }
            let source = mesh.links[(h ^ 1) as usize].target;
            if mesh.links[link.prev as usize].target != source {
                return Err(fail(
                    "prev does not end at the source",
                    "halfedge",
                    h as usize,
                ));
            }
            if source == link.target {
                return Err(fail("joins a vertex to itself", "halfedge", h as usize));
            }
            out_count[source as usize] += 1;
            if link.face == NONE {
                boundary_out[source as usize] += 1;
            }
        }
    }

    for f in mesh.faces() {
        let h = mesh.face_halfedge[f.index()];
        if !halfedge_live(h) || mesh.links[h as usize].face != f.0 {
            return Err(fail(
                "stored halfedge is not on the face",
                "face",
                f.index(),
            ));
        }
        let mut corners: Vec<u32> = Vec::new();
        for g in mesh.loop_halfedges(HalfedgeId(h)) {
            corners.push(mesh.links[g.index()].target);
            if corners.len() > mesh.links.len() {
                return Err(fail("loop does not close", "face", f.index()));
            }
        }
        if corners.len() < 3 {
            return Err(fail("fewer than three sides", "face", f.index()));
        }
        corners.sort_unstable();
        if corners.windows(2).any(|pair| pair[0] == pair[1]) {
            return Err(fail("visits a vertex twice", "face", f.index()));
        }
    }

    let mut neighbours: Vec<VertexId> = Vec::new();
    for v in mesh.vertices() {
        let h = mesh.vertex_out[v.index()];
        if h == NONE {
            if out_count[v.index()] != 0 {
                return Err(fail(
                    "has edges but no stored halfedge",
                    "vertex",
                    v.index(),
                ));
            }
            continue;
        }
        if !halfedge_live(h) || mesh.links[(h ^ 1) as usize].target != v.0 {
            return Err(fail(
                "stored halfedge does not leave it",
                "vertex",
                v.index(),
            ));
        }
        if boundary_out[v.index()] > 1 {
            return Err(fail("lies on more than one hole", "vertex", v.index()));
        }
        if boundary_out[v.index()] == 1 && mesh.links[h as usize].face != NONE {
            return Err(fail(
                "stored halfedge is not its boundary halfedge",
                "vertex",
                v.index(),
            ));
        }
        neighbours.clear();
        for g in mesh.outgoing_halfedges(v) {
            neighbours.push(mesh.target(g));
            if neighbours.len() > out_count[v.index()] as usize {
                break;
            }
        }
        if neighbours.len() != out_count[v.index()] as usize {
            return Err(fail(
                "rotation does not reach every outgoing halfedge",
                "vertex",
                v.index(),
            ));
        }
        neighbours.sort_unstable();
        if neighbours.windows(2).any(|pair| pair[0] == pair[1]) {
            return Err(fail(
                "two edges join the same neighbour",
                "vertex",
                v.index(),
            ));
        }
    }
    Ok(())
}
