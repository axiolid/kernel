// SPDX-License-Identifier: MPL-2.0

//! Traversal and measurement over the arrangement.
//!
//! Everything here is derived on demand. Nothing is cached, so no edit can
//! leave a stale answer behind.

use axiolid_core::Point2;

use crate::id::{FaceId, HalfEdgeId, VertexId};
use crate::Arrangement;

impl Arrangement {
    /// Half-edges around a face's boundary, counter-clockwise.
    ///
    /// Returns empty for a face with no boundary (the unbounded face of an
    /// empty arrangement).
    #[must_use]
    pub fn face_halfedges(&self, face: FaceId) -> Vec<HalfEdgeId> {
        let Some(start) = self.faces[face.index()].boundary else {
            return Vec::new();
        };
        let mut out = vec![start];
        let mut current = self.halfedges[start.index()].next;
        // Bounded by the arena: a corrupted `next` cycle cannot hang the
        // caller, it just yields a short walk that `validate` will flag.
        while current != start && out.len() <= self.halfedges.len() {
            out.push(current);
            current = self.halfedges[current.index()].next;
        }
        out
    }

    /// Corner positions of a face, counter-clockwise.
    #[must_use]
    pub fn face_outline(&self, face: FaceId) -> Vec<Point2> {
        self.face_halfedges(face)
            .into_iter()
            .map(|h| self.vertices[self.halfedges[h.index()].origin.index()].position)
            .collect()
    }

    /// Signed area of a face, positive when its boundary runs
    /// counter-clockwise.
    ///
    /// The unbounded face has no meaningful area; it reports 0.0 rather than
    /// a negative number that a caller might sum into a total.
    #[must_use]
    pub fn face_area(&self, face: FaceId) -> f64 {
        if face == FaceId::OUTER {
            return 0.0;
        }
        let outline = self.face_outline(face);
        if outline.len() < 3 {
            return 0.0;
        }
        // Shoelace about the first vertex rather than the world origin: same
        // reasoning as the mass-properties fix, and free here.
        let base = outline[0];
        let mut twice = 0.0;
        for window in outline[1..].windows(2) {
            let a = window[0] - base;
            let b = window[1] - base;
            twice += a.x * b.y - a.y * b.x;
        }
        twice / 2.0
    }

    /// The face on the other side of a half-edge.
    ///
    /// Always defined: the unbounded face is a real face, so an edge on the
    /// outer boundary reports it rather than `None`.
    #[must_use]
    pub fn neighbour_across(&self, edge: HalfEdgeId) -> FaceId {
        let twin = self.halfedges[edge.index()].twin;
        self.halfedges[twin.index()].face
    }

    /// Every bounded face, in arena order.
    ///
    /// Deliberately excludes the unbounded face: a caller asking for "the
    /// regions" almost never means the infinite one, and including it is the
    /// kind of default that produces a wrong total on the first use.
    pub fn bounded_faces(&self) -> impl Iterator<Item = FaceId> + '_ {
        (1..self.faces.len()).map(FaceId::from_index)
    }

    /// Half-edges leaving a vertex, counter-clockwise around it.
    #[must_use]
    pub fn vertex_halfedges(&self, vertex: VertexId) -> Vec<HalfEdgeId> {
        let Some(start) = self.vertices[vertex.index()].outgoing else {
            return Vec::new();
        };
        let mut out = vec![start];
        // Around a vertex: take the twin (now pointing in), then its next
        // (pointing out again, one step around).
        let mut current = self.halfedges[self.halfedges[start.index()].twin.index()].next;
        while current != start && out.len() <= self.halfedges.len() {
            out.push(current);
            current = self.halfedges[self.halfedges[current.index()].twin.index()].next;
        }
        out
    }

    /// Number of edges meeting at a vertex.
    #[must_use]
    pub fn degree(&self, vertex: VertexId) -> usize {
        self.vertex_halfedges(vertex).len()
    }

    /// Vertex a half-edge leaves from.
    #[must_use]
    pub fn halfedge_origin(&self, edge: HalfEdgeId) -> VertexId {
        self.halfedges[edge.index()].origin
    }

    /// Face lying to the left of a half-edge.
    #[must_use]
    pub fn halfedge_face(&self, edge: HalfEdgeId) -> FaceId {
        self.halfedges[edge.index()].face
    }
}
