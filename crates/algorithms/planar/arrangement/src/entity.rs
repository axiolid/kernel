// SPDX-License-Identifier: MPL-2.0

//! Arena records for the arrangement.
//!
//! Records hold only topology and position. Anything derived -- area,
//! outline, adjacency -- is computed in `query.rs` rather than cached here,
//! because a cached derivative is a second source of truth that edits have to
//! remember to invalidate.

use axiolid_core::Point2;

use crate::id::{FaceId, HalfEdgeId, VertexId};

impl FaceId {
    /// The unbounded face, which every arrangement has from creation.
    pub const OUTER: Self = Self::outer();
}

/// A point in the plane, with one outgoing half-edge as an entry into the
/// edges around it.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Vertex {
    /// Where the vertex is.
    pub position: Point2,
    /// Any half-edge whose origin is this vertex, or `None` if isolated.
    ///
    /// Isolated vertices are legal: an edit can remove the last edge touching
    /// a vertex without the vertex itself becoming invalid.
    pub outgoing: Option<HalfEdgeId>,
}

/// One direction of an edge.
///
/// The twin runs the other way along the same geometric segment. Storing both
/// directions is what makes "the face on the other side" a constant-time
/// question rather than a search.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct HalfEdge {
    /// Vertex this half-edge leaves from.
    pub origin: VertexId,
    /// The opposite half-edge along the same segment.
    pub twin: HalfEdgeId,
    /// Next half-edge counter-clockwise around `face`.
    pub next: HalfEdgeId,
    /// Previous half-edge around `face`.
    pub prev: HalfEdgeId,
    /// Face lying to the left of this half-edge.
    pub face: FaceId,
}

/// A region of the plane bounded by a cycle of half-edges.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Face {
    /// Any half-edge on this face's outer boundary.
    ///
    /// `None` only for the unbounded face of an arrangement with no edges.
    pub boundary: Option<HalfEdgeId>,
}
