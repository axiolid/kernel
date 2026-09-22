// SPDX-License-Identifier: MPL-2.0
#![forbid(unsafe_code)]
#![warn(missing_docs)]

//! Editable planar subdivision with persistent half-edge topology.
//!
//! # What this is for
//!
//! [`axiolid_overlay`] answers "what is the union of these polygons" in one
//! shot: polygons in, polygons out, no structure retained. That is the right
//! shape for a query, and the wrong shape for editing. A caller who moves one
//! vertex has to rebuild everything and then re-derive which output polygon
//! corresponds to which input -- identity is lost on every call.
//!
//! This crate keeps the subdivision itself. Vertices, half-edges and faces
//! have stable handles that survive edits, so "this face" means the same face
//! before and after a vertex moves, and an edit touches only the affected
//! neighbourhood instead of rebuilding the plane.
//!
//! # Deliberately neutral
//!
//! A planar arrangement is a general structure: it does not know about rooms,
//! walls, storeys, or net floor area. Those are domain concepts and belong to
//! the consumer that has the domain. This crate exposes faces, their
//! boundaries, their areas, and their adjacencies; deciding that a particular
//! face is a room is the caller's judgement, made with information this crate
//! does not have.
//!
//! # Structure
//!
//! Standard doubly-connected edge list. Each edge is two opposite half-edges;
//! each half-edge knows its origin vertex, its twin, and the next half-edge
//! around its face. A face is identified by any half-edge on its boundary.
//! Walking `next` traverses a face's boundary; walking `twin`/`next`
//! traverses the edges around a vertex.
//!
//! The unbounded outer region is a real face ([`Arrangement::outer_face`]),
//! not a `None`. Making it explicit removes a special case from every
//! traversal: "the face across this edge" always has an answer.

use axiolid_core::Point2;

mod build;
mod edit;
mod entity;
mod id;
mod query;
mod validate;

pub use build::BuildError;
pub use edit::EditError;
pub use entity::{Face, HalfEdge, Vertex};
pub use id::{FaceId, HalfEdgeId, VertexId};
pub use validate::ArrangementHealth;

/// A planar subdivision as a doubly-connected edge list.
///
/// Handles stay valid across edits unless the element they name is removed,
/// which is what makes incremental editing possible at all.
#[derive(Debug, Clone)]
pub struct Arrangement {
    pub(crate) vertices: Vec<Vertex>,
    pub(crate) halfedges: Vec<HalfEdge>,
    pub(crate) faces: Vec<Face>,
}

impl Arrangement {
    /// An empty plane: one unbounded face, no vertices or edges.
    #[must_use]
    pub fn new() -> Self {
        Self {
            vertices: Vec::new(),
            halfedges: Vec::new(),
            // The unbounded face exists from the start, so `outer_face` is
            // always a valid handle and callers never special-case an empty
            // arrangement.
            faces: vec![Face { boundary: None }],
        }
    }

    /// The unbounded region surrounding every bounded face.
    #[must_use]
    pub const fn outer_face(&self) -> FaceId {
        FaceId::OUTER
    }

    /// Number of vertices, including any left isolated by edits.
    #[must_use]
    pub fn vertex_count(&self) -> usize {
        self.vertices.len()
    }

    /// Number of faces, including the unbounded one.
    #[must_use]
    pub fn face_count(&self) -> usize {
        self.faces.len()
    }

    /// Number of half-edges; always twice the number of edges.
    #[must_use]
    pub fn halfedge_count(&self) -> usize {
        self.halfedges.len()
    }

    /// Position of a vertex.
    #[must_use]
    pub fn position(&self, vertex: VertexId) -> Point2 {
        self.vertices[vertex.index()].position
    }
}

impl Default for Arrangement {
    fn default() -> Self {
        Self::new()
    }
}
