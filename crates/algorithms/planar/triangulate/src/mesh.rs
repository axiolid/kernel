// SPDX-License-Identifier: MPL-2.0

//! The triangulation structure and the constrained Delaunay build.
//!
//! # Representation
//!
//! Triangles are stored in a flat array with a parallel neighbour array, the
//! standard "triangle + halfedge" layout: halfedge `3t + i` belongs to
//! triangle `t`, runs from vertex `i` to vertex `i + 1 mod 3`, and its twin
//! is `halfedge[3t + i]`. A half-edge structure with explicit records would
//! carry a pointer per edge and buy nothing here -- the triangulation is
//! rebuilt rather than edited, so compactness and cache locality win.
//!
//! Vertices 0..n are the caller's points. Three extra vertices form a super
//! triangle large enough to contain them all; they are removed at the end,
//! together with every triangle that touches them.

use axiolid_core::Point2;

use crate::Constraint;

/// Sentinel for "no neighbour across this halfedge".
pub(crate) const NO_HALFEDGE: u32 = u32::MAX;

/// Why a triangulation could not be produced.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
#[non_exhaustive]
pub enum TriangulationError {
    /// Fewer than three points were supplied, so no triangle exists.
    TooFewPoints,
    /// Every input point is collinear, so the triangulation has zero area.
    ///
    /// Reported rather than returning an empty triangulation: a caller that
    /// asked to triangulate a degenerate polygon has a bug upstream, and an
    /// empty result would hide it.
    AllPointsCollinear,
    /// A vertex could not be placed in any triangle during insertion.
    ///
    /// Indicates a corrupted adjacency rather than bad input. Reported so a
    /// missing vertex surfaces as an error instead of as a silent hole in an
    /// otherwise valid-looking mesh.
    VertexUnplaceable {
        /// Index of the vertex that could not be placed.
        index: u32,
    },
    /// A constraint referenced a vertex index that does not exist.
    ConstraintOutOfRange {
        /// The offending index.
        index: u32,
    },
    /// A constraint could not be recovered after insertion.
    ///
    /// Only reachable when two constraints cross: an edge cannot survive if
    /// another required edge passes through it. The crossing point would have
    /// to be inserted as a vertex, which changes the caller's input, so this
    /// is reported instead of silently repaired.
    CrossingConstraints {
        /// First endpoint of the constraint that could not be recovered.
        a: u32,
        /// Second endpoint.
        b: u32,
    },
}

impl core::fmt::Display for TriangulationError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        match self {
            Self::TooFewPoints => write!(f, "a triangulation needs at least three points"),
            Self::AllPointsCollinear => write!(f, "every input point is collinear"),
            Self::VertexUnplaceable { index } => {
                write!(f, "vertex {index} could not be placed in any triangle")
            }
            Self::ConstraintOutOfRange { index } => {
                write!(f, "constraint references out-of-range vertex {index}")
            }
            Self::CrossingConstraints { a, b } => {
                write!(f, "constraint ({a}, {b}) crosses another constraint")
            }
        }
    }
}

impl core::error::Error for TriangulationError {}

/// A planar triangulation over a point set, honouring constraint edges.
#[derive(Debug, Clone)]
pub struct Triangulation {
    pub(crate) points: Vec<Point2>,
    /// Three vertex indices per triangle, counter-clockwise.
    pub(crate) triangles: Vec<u32>,
    /// Twin halfedge per halfedge, or [`NO_HALFEDGE`].
    pub(crate) halfedges: Vec<u32>,
    /// Constraint edges, sorted, as vertex index pairs.
    pub(crate) constraints: Vec<Constraint>,
}

impl Triangulation {
    /// The triangle list, three counter-clockwise vertex indices each.
    #[must_use]
    pub fn triangles(&self) -> &[u32] {
        &self.triangles
    }

    /// The vertex positions, including any Steiner points appended by
    /// refinement.
    #[must_use]
    pub fn points(&self) -> &[Point2] {
        &self.points
    }

    /// Number of triangles.
    #[must_use]
    pub fn triangle_count(&self) -> usize {
        self.triangles.len() / 3
    }

    /// The constraint edges this triangulation was built to honour.
    #[must_use]
    pub fn constraints(&self) -> &[Constraint] {
        &self.constraints
    }

    /// Whether `edge` is a constraint.
    pub(crate) fn is_constrained(&self, a: u32, b: u32) -> bool {
        self.constraints
            .binary_search(&Constraint::new(a, b))
            .is_ok()
    }
}
