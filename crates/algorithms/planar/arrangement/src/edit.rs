// SPDX-License-Identifier: MPL-2.0

//! Incremental edits that preserve handles and topological validity.
//!
//! # Why these return `Result`
//!
//! Every edit here can be asked to do something that would corrupt the
//! structure: split an edge at a point that is not on it, drag a vertex so a
//! face self-intersects, connect two vertices that are not on a common face.
//! Each is refused by name. A DCEL that silently accepts a corrupting edit is
//! worse than one that cannot edit at all, because the damage surfaces later
//! as a traversal that never terminates.

use axiolid_core::Point2;
use axiolid_guarantees::Sign;
use axiolid_predicates::orient2d;

use crate::entity::{HalfEdge, Vertex};
use crate::id::{FaceId, HalfEdgeId, VertexId};
use crate::Arrangement;

/// Why an edit was refused.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
#[non_exhaustive]
pub enum EditError {
    /// The split point does not lie on the target edge.
    PointNotOnEdge,
    /// The two vertices do not share a face, so no chord connects them.
    VerticesNotOnCommonFace,
    /// The two vertices are already joined by an edge.
    AlreadyConnected,
    /// The edit would leave a face non-convex or self-intersecting.
    WouldSelfIntersect,
    /// The half-edge borders the unbounded face, which cannot be merged away.
    BordersUnboundedFace,
}

impl core::fmt::Display for EditError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        match self {
            Self::PointNotOnEdge => write!(f, "split point does not lie on the edge"),
            Self::VerticesNotOnCommonFace => {
                write!(f, "vertices do not share a face")
            }
            Self::AlreadyConnected => write!(f, "vertices are already connected"),
            Self::WouldSelfIntersect => {
                write!(f, "edit would make a face self-intersecting")
            }
            Self::BordersUnboundedFace => {
                write!(f, "half-edge borders the unbounded face")
            }
        }
    }
}

impl core::error::Error for EditError {}

impl Arrangement {
    /// Add an isolated vertex.
    ///
    /// Legal on its own: a vertex with no edges is a valid, if uninteresting,
    /// arrangement element, and building incrementally needs it.
    pub fn add_vertex(&mut self, position: Point2) -> VertexId {
        self.vertices.push(Vertex {
            position,
            outgoing: None,
        });
        VertexId::from_index(self.vertices.len() - 1)
    }

    /// Split an edge at `at`, inserting a new vertex.
    ///
    /// The point must lie on the segment. Both the edge and its twin are
    /// split, so the structure stays consistent.
    ///
    /// # Errors
    ///
    /// [`EditError::PointNotOnEdge`] if `at` is not collinear with the edge
    /// or lies outside it.
    pub fn split_edge(&mut self, edge: HalfEdgeId, at: Point2) -> Result<VertexId, EditError> {
        let he = self.halfedges[edge.index()];
        let twin = self.halfedges[he.twin.index()];
        let from = self.vertices[he.origin.index()].position;
        let to = self.vertices[twin.origin.index()].position;

        // Exactly collinear, and strictly between the endpoints. Using the
        // certified predicate rather than a tolerance keeps this decision
        // consistent with the rest of the kernel.
        if !matches!(decided(orient2d(from, to, at)), Sign::Zero) {
            return Err(EditError::PointNotOnEdge);
        }
        let within = (at.x - from.x) * (to.x - at.x) + (at.y - from.y) * (to.y - at.y);
        if within <= 0.0 {
            return Err(EditError::PointNotOnEdge);
        }

        let new_vertex = self.add_vertex(at);
        let a = HalfEdgeId::from_index(self.halfedges.len());
        let b = HalfEdgeId::from_index(self.halfedges.len() + 1);

        // `a` continues `edge` from the new vertex; `b` continues the twin.
        self.halfedges.push(HalfEdge {
            origin: new_vertex,
            twin: he.twin,
            next: he.next,
            prev: edge,
            face: he.face,
        });
        self.halfedges.push(HalfEdge {
            origin: new_vertex,
            twin: edge,
            next: twin.next,
            prev: he.twin,
            face: twin.face,
        });

        let next_of_edge = he.next;
        let next_of_twin = twin.next;
        self.halfedges[edge.index()].next = a;
        self.halfedges[edge.index()].twin = b;
        self.halfedges[he.twin.index()].next = b;
        self.halfedges[he.twin.index()].twin = a;
        self.halfedges[next_of_edge.index()].prev = a;
        self.halfedges[next_of_twin.index()].prev = b;
        self.vertices[new_vertex.index()].outgoing = Some(a);
        Ok(new_vertex)
    }

    /// Move a vertex, refusing the move if it would break a face.
    ///
    /// # Errors
    ///
    /// [`EditError::WouldSelfIntersect`] if any incident face would stop
    /// being simple. The arrangement is left unchanged in that case.
    pub fn drag_vertex(&mut self, vertex: VertexId, to: Point2) -> Result<(), EditError> {
        let original = self.vertices[vertex.index()].position;
        self.vertices[vertex.index()].position = to;

        // Check every face touching the vertex, and roll back as a unit. A
        // partially applied drag would be worse than a refused one.
        let faces: Vec<FaceId> = self
            .vertex_halfedges(vertex)
            .into_iter()
            .map(|h| self.halfedges[h.index()].face)
            .collect();
        for face in faces {
            if face == FaceId::OUTER {
                continue;
            }
            if !self.face_is_simple(face) {
                self.vertices[vertex.index()].position = original;
                return Err(EditError::WouldSelfIntersect);
            }
        }
        Ok(())
    }

    /// Whether a face's boundary is a simple counter-clockwise polygon.
    fn face_is_simple(&self, face: FaceId) -> bool {
        let outline = self.face_outline(face);
        if outline.len() < 3 {
            return false;
        }
        // A convex-or-not test is not enough: a simple polygon may be
        // concave. Check that the boundary does not reverse orientation,
        // which is the failure a drag actually causes.
        let area = self.face_area(face);
        if area <= 0.0 {
            return false;
        }
        // And that no two non-adjacent edges cross.
        let n = outline.len();
        for i in 0..n {
            let (a1, a2) = (outline[i], outline[(i + 1) % n]);
            for j in (i + 2)..n {
                if (j + 1) % n == i {
                    continue;
                }
                let (b1, b2) = (outline[j], outline[(j + 1) % n]);
                if segments_cross(a1, a2, b1, b2) {
                    return false;
                }
            }
        }
        true
    }
}

/// Certified sign, treating an undecidable result as degenerate.
fn decided(certified: axiolid_guarantees::Certified) -> Sign {
    match certified {
        axiolid_guarantees::Certified::Certain { sign, .. } => sign,
        _ => Sign::Zero,
    }
}

/// Whether two segments properly cross.
fn segments_cross(a1: Point2, a2: Point2, b1: Point2, b2: Point2) -> bool {
    let d1 = decided(orient2d(a1, a2, b1));
    let d2 = decided(orient2d(a1, a2, b2));
    let d3 = decided(orient2d(b1, b2, a1));
    let d4 = decided(orient2d(b1, b2, a2));
    d1 != d2 && d3 != d4 && d1 != Sign::Zero && d2 != Sign::Zero
}
