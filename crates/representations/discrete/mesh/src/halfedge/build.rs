//! Building a halfedge mesh from an indexed face list, refusing what a
//! halfedge structure cannot represent.

use core::fmt;

use axiolid_core::Point3;

use super::{HalfedgeMesh, Link, VertexId, NONE};

/// Why a face list cannot become a [`HalfedgeMesh`].
///
/// Each refusal names the configuration and the elements involved, so a
/// caller can repair or report the input instead of guessing.
#[non_exhaustive]
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum HalfedgeBuildError {
    /// The triangle index buffer is not divisible by three.
    IncompleteTriangle {
        /// Indices found.
        index_count: usize,
    },
    /// A face has fewer than three corners.
    FaceTooSmall {
        /// The face.
        face: usize,
        /// Corners found.
        corners: usize,
    },
    /// A corner names a position that does not exist.
    IndexOutOfRange {
        /// The face.
        face: usize,
        /// The offending index.
        index: u32,
        /// Positions available.
        position_count: usize,
    },
    /// A face visits the same vertex twice.
    DegenerateFace {
        /// The face.
        face: usize,
        /// The repeated vertex.
        vertex: u32,
    },
    /// More than two faces share an edge.
    NonManifoldEdge {
        /// Lower endpoint.
        a: u32,
        /// Higher endpoint.
        b: u32,
        /// Faces using the edge.
        faces: usize,
    },
    /// Two faces traverse a shared edge in the same direction: their windings
    /// disagree (or the same face appears twice).
    InconsistentOrientation {
        /// Source of the doubly-used direction.
        from: u32,
        /// Target of the doubly-used direction.
        to: u32,
    },
    /// The faces around a vertex form more than one fan: two cones touching
    /// at a point, or a boundary vertex where two holes meet.
    NonManifoldVertex {
        /// The vertex.
        vertex: u32,
    },
    /// The mesh needs more elements than 32-bit ids can address.
    TooLarge,
}

impl fmt::Display for HalfedgeBuildError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::IncompleteTriangle { index_count } => {
                write!(f, "index count {index_count} is not divisible by three")
            }
            Self::FaceTooSmall { face, corners } => {
                write!(f, "face {face} has {corners} corners, fewer than three")
            }
            Self::IndexOutOfRange {
                face,
                index,
                position_count,
            } => write!(
                f,
                "face {face} names position {index}, but there are {position_count}"
            ),
            Self::DegenerateFace { face, vertex } => {
                write!(f, "face {face} visits vertex {vertex} twice")
            }
            Self::NonManifoldEdge { a, b, faces } => {
                write!(f, "edge {a}-{b} is shared by {faces} faces")
            }
            Self::InconsistentOrientation { from, to } => write!(
                f,
                "two faces traverse edge {from}->{to} in the same direction"
            ),
            Self::NonManifoldVertex { vertex } => {
                write!(f, "vertex {vertex} joins more than one fan of faces")
            }
            Self::TooLarge => write!(f, "mesh exceeds 32-bit element ids"),
        }
    }
}

impl std::error::Error for HalfedgeBuildError {}

/// One face corner's directed edge, sortable by its undirected key.
#[derive(Clone, Copy)]
struct Use {
    lower: u32,
    upper: u32,
    from: u32,
    to: u32,
    corner: u32,
}

pub(super) fn build<F: AsRef<[u32]>>(
    positions: Vec<Point3>,
    faces: &[F],
) -> Result<HalfedgeMesh, HalfedgeBuildError> {
    let position_count = positions.len();
    let corner_count: usize = faces.iter().map(|face| face.as_ref().len()).sum();
    // Every corner contributes at most one edge, so `2 * corners` halfedges
    // bounds the build; ids must stay below the `NONE` sentinel.
    if position_count >= NONE as usize
        || faces.len() >= NONE as usize
        || corner_count.saturating_mul(2) >= NONE as usize
    {
        return Err(HalfedgeBuildError::TooLarge);
    }

    let mut face_start = Vec::with_capacity(faces.len() + 1);
    let mut uses = Vec::with_capacity(corner_count);
    let mut sorted = Vec::new();
    for (face_index, face) in faces.iter().enumerate() {
        let face = face.as_ref();
        if face.len() < 3 {
            return Err(HalfedgeBuildError::FaceTooSmall {
                face: face_index,
                corners: face.len(),
            });
        }
        if let Some(&index) = face.iter().find(|&&i| i as usize >= position_count) {
            return Err(HalfedgeBuildError::IndexOutOfRange {
                face: face_index,
                index,
                position_count,
            });
        }
        sorted.clear();
        sorted.extend_from_slice(face);
        sorted.sort_unstable();
        if let Some(pair) = sorted.windows(2).find(|pair| pair[0] == pair[1]) {
            return Err(HalfedgeBuildError::DegenerateFace {
                face: face_index,
                vertex: pair[0],
            });
        }
        face_start.push(uses.len() as u32);
        for (i, &from) in face.iter().enumerate() {
            let to = face[(i + 1) % face.len()];
            uses.push(Use {
                lower: from.min(to),
                upper: from.max(to),
                from,
                to,
                corner: uses.len() as u32,
            });
        }
    }
    face_start.push(uses.len() as u32);

    // Group corners by undirected edge; within a group, by corner so the
    // first face to use an edge decides its halfedge `2e`.
    let mut order = uses.clone();
    order.sort_unstable_by_key(|u| (u.lower, u.upper, u.corner));

    let mut corner_halfedge = vec![NONE; corner_count];
    let mut links: Vec<Link> = Vec::with_capacity(corner_count * 2);
    let mut group_start = 0;
    while group_start < order.len() {
        let first = order[group_start];
        let mut group_end = group_start + 1;
        while group_end < order.len()
            && (order[group_end].lower, order[group_end].upper) == (first.lower, first.upper)
        {
            group_end += 1;
        }
        let group = &order[group_start..group_end];
        if group.len() > 2 {
            return Err(HalfedgeBuildError::NonManifoldEdge {
                a: first.lower,
                b: first.upper,
                faces: group.len(),
            });
        }
        if group.len() == 2 && group[1].from == first.from {
            return Err(HalfedgeBuildError::InconsistentOrientation {
                from: first.from,
                to: first.to,
            });
        }
        let h = links.len() as u32;
        // Boundary until a face claims the side.
        links.push(Link {
            target: first.to,
            ..super::DEAD_LINK
        });
        links.push(Link {
            target: first.from,
            ..super::DEAD_LINK
        });
        corner_halfedge[first.corner as usize] = h;
        if let Some(second) = group.get(1) {
            corner_halfedge[second.corner as usize] = h + 1;
        }
        group_start = group_end;
    }

    // Face loops.
    let mut face_halfedge = Vec::with_capacity(faces.len());
    for face in 0..faces.len() {
        let start = face_start[face] as usize;
        let end = face_start[face + 1] as usize;
        let len = end - start;
        for corner in start..end {
            let h = corner_halfedge[corner] as usize;
            let next = corner_halfedge[start + (corner - start + 1) % len];
            let prev = corner_halfedge[start + (corner - start + len - 1) % len];
            links[h].next = next;
            links[h].prev = prev;
            links[h].face = face as u32;
        }
        face_halfedge.push(corner_halfedge[start]);
    }

    // Boundary loops: each vertex may leave by at most one boundary halfedge.
    let mut vertex_out = vec![NONE; position_count];
    let mut boundary_out = vec![NONE; position_count];
    let mut out_count = vec![0u32; position_count];
    for h in 0..links.len() {
        let source = links[h ^ 1].target as usize;
        out_count[source] += 1;
        if links[h].face == NONE {
            if boundary_out[source] != NONE {
                return Err(HalfedgeBuildError::NonManifoldVertex {
                    vertex: source as u32,
                });
            }
            boundary_out[source] = h as u32;
        }
        if vertex_out[source] == NONE {
            vertex_out[source] = h as u32;
        }
    }
    for h in 0..links.len() {
        if links[h].face == NONE {
            // A vertex has as many incoming as outgoing boundary halfedges,
            // so the target of a boundary halfedge always leaves by one.
            let next = boundary_out[links[h].target as usize];
            links[h].next = next;
            links[next as usize].prev = h as u32;
        }
    }
    for v in 0..position_count {
        if boundary_out[v] != NONE {
            vertex_out[v] = boundary_out[v];
        }
    }

    let edge_count = links.len() / 2;
    let mesh = HalfedgeMesh {
        positions,
        vertex_live: vec![true; position_count],
        vertex_out,
        edge_live: vec![true; edge_count],
        links,
        face_live: vec![true; face_halfedge.len()],
        live_faces: face_halfedge.len(),
        face_halfedge,
        live_vertices: position_count,
        live_edges: edge_count,
    };

    // A manifold vertex's outgoing halfedges form one rotation orbit.
    for (v, &count) in out_count.iter().enumerate() {
        if count == 0 {
            continue;
        }
        let orbit = mesh.outgoing_halfedges(VertexId(v as u32)).count() as u32;
        if orbit != count {
            return Err(HalfedgeBuildError::NonManifoldVertex { vertex: v as u32 });
        }
    }
    debug_assert!(mesh.validate().is_ok(), "{:?}", mesh.validate());
    Ok(mesh)
}
