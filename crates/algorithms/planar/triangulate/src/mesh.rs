// SPDX-License-Identifier: MPL-2.0

//! The triangulation structure and the constrained Delaunay build.
//!
//! # Representation
//!
//! Triangles are stored in a flat array with a parallel neighbour array, the
//! standard "triangle + halfedge" layout: halfedge `3t + i` belongs to
//! triangle `t`, runs from vertex `i` to vertex `i + 1 mod 3`, and its twin
//! is `halfedge[3t + i]`. A half-edge structure with explicit records would
//! carry a pointer per edge and buy nothing here: compactness and cache
//! locality win. Insertion and flips keep the twin array current locally;
//! the `adjacency` module holds the rules and the debug check against a
//! full recomputation.
//!
//! Vertices 0..n are the caller's points. During the build one more
//! vertex, at infinity and without coordinates, closes the convex hull
//! (see the `build` module); it is removed at the end, together with
//! every triangle that touches it, and the rest covers the convex hull.

use axiolid_core::Point2;

use crate::domain::Ring;
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
    /// Two constraint pieces cross; a constraint is split at the vertices
    /// on it first.
    ///
    /// Only reachable when two constraints cross at a point that is no
    /// input vertex: an edge cannot survive if another required edge passes
    /// through it. The crossing point would have to be inserted as a
    /// vertex, which changes the caller's input, so this is reported
    /// instead of silently repaired. A constraint that runs through input
    /// vertices, or ends on another constraint, is not a crossing.
    CrossingConstraints {
        /// First endpoint of the piece that could not be recovered.
        a: u32,
        /// Second endpoint.
        b: u32,
    },
    /// A constraint joins two input points at exactly the same position.
    ///
    /// Coincident points are merged into one vertex, so the constraint
    /// would be an edge from that vertex to itself. Reported rather than
    /// dropped: a caller who asked for an edge between two distinct
    /// indices has two points where they meant one, or one where they
    /// meant two.
    CoincidentConstraint {
        /// First endpoint, as an input index.
        a: u32,
        /// Second endpoint, as an input index.
        b: u32,
    },
    /// A ring referenced a vertex index that does not exist.
    RingOutOfRange {
        /// Position of the ring in the ring list.
        ring: u32,
        /// The offending index.
        index: u32,
    },
    /// A ring has fewer than three vertices, so it encloses nothing.
    TooFewRingVertices {
        /// Position of the ring in the ring list.
        ring: u32,
    },
    /// A ring visits the same vertex twice.
    ///
    /// Two indices at one position count as one vertex, and so does
    /// repeating the first index at the end to close the ring: a ring is
    /// closed implicitly.
    RepeatedRingVertex {
        /// Position of the ring in the ring list.
        ring: u32,
        /// The later of the two occurrences, as an input index.
        index: u32,
    },
    /// A ring doubles back on itself at a vertex.
    ///
    /// Both neighbours of the vertex lie on one side of it on a common
    /// line, so the ring has a zero-width spike (or no area at all) there.
    FoldedRing {
        /// Position of the ring in the ring list.
        ring: u32,
        /// The vertex where the ring folds, as an input index.
        index: u32,
    },
    /// The rings nest in a way that is neither inside nor outside: a region
    /// has a winding other than 0 or 1.
    ///
    /// A hole outside every outer ring, a hole inside a hole, a hole listed
    /// twice and an outer ring inside another all end here. Rings (or a
    /// ring and itself) that cross between vertices are refused earlier, as
    /// [`CrossingConstraints`](Self::CrossingConstraints).
    MisnestedRing {
        /// The lowest-numbered ring that carries the boundary crossed into
        /// the region.
        ring: u32,
        /// The winding of that region.
        winding: i32,
    },
    /// The rings enclose no triangle, so the domain is empty.
    EmptyDomain,
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
            Self::CoincidentConstraint { a, b } => {
                write!(f, "constraint ({a}, {b}) joins two coincident points")
            }
            Self::RingOutOfRange { ring, index } => {
                write!(f, "ring {ring} references out-of-range vertex {index}")
            }
            Self::TooFewRingVertices { ring } => {
                write!(f, "ring {ring} has fewer than three vertices")
            }
            Self::RepeatedRingVertex { ring, index } => {
                write!(f, "ring {ring} visits vertex {index} twice")
            }
            Self::FoldedRing { ring, index } => {
                write!(f, "ring {ring} doubles back on itself at vertex {index}")
            }
            Self::MisnestedRing { ring, winding } => {
                write!(
                    f,
                    "ring {ring} borders a region of winding {winding}, not 0 or 1"
                )
            }
            Self::EmptyDomain => write!(f, "the rings enclose no triangle"),
        }
    }
}

impl core::error::Error for TriangulationError {}

/// An input point merged into an exactly coincident one.
///
/// Two input points at the same position are one vertex of the
/// triangulation: the lowest input index at that position is kept, and the
/// others are reported here. No triangle uses a merged index.
#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Hash)]
#[non_exhaustive]
pub struct MergedPoint {
    /// The input index that was merged.
    pub input: u32,
    /// The input index that stands for it in the triangles: the lowest
    /// index at the same position.
    pub kept: u32,
}

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
    /// Coincident input points merged into a kept one, by input index.
    pub(crate) merged: Vec<MergedPoint>,
    /// Kept vertices that no triangle uses, ascending; set by the domain
    /// filter.
    pub(crate) outside: Vec<u32>,
    /// The rings the triangles were filtered by, so refinement can filter
    /// again; `None` for a hull triangulation.
    pub(crate) domain: Option<Vec<Ring>>,
}

impl Triangulation {
    /// The triangle list, three counter-clockwise vertex indices each.
    ///
    /// Indices are into [`points`](Self::points); an input point merged
    /// into a coincident one (see [`merged_points`](Self::merged_points))
    /// appears under the kept index.
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

    /// The constraint edges this triangulation was built to honour, between
    /// kept vertices: an endpoint that was merged into a coincident point
    /// is replaced by the kept index. A constraint through input vertices
    /// is listed as the pieces between them, each an edge of the output.
    /// For a domain triangulation, pieces outside the domain are dropped.
    #[must_use]
    pub fn constraints(&self) -> &[Constraint] {
        &self.constraints
    }

    /// The input points merged into an exactly coincident point, in
    /// ascending input order; empty when every input point is distinct.
    ///
    /// Every input index not listed is its own vertex. With this list the
    /// merge is lossless: a caller can map each of its indices to the
    /// vertex that carries it.
    #[must_use]
    pub fn merged_points(&self) -> &[MergedPoint] {
        &self.merged
    }

    /// The kept vertices that no triangle uses, ascending: the points of a
    /// domain triangulation that lie outside the domain (or inside a hole)
    /// and were dropped, and ring vertices that bound no domain triangle
    /// (an outline corner covered by a hole). Empty for a hull triangulation, where every kept
    /// point is a corner. A merged point is not listed; see
    /// [`merged_points`](Self::merged_points).
    #[must_use]
    pub fn outside_points(&self) -> &[u32] {
        &self.outside
    }

    /// Whether `edge` is a constraint.
    pub(crate) fn is_constrained(&self, a: u32, b: u32) -> bool {
        self.constraints
            .binary_search(&Constraint::new(a, b))
            .is_ok()
    }
}
