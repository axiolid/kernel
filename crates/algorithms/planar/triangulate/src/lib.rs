// SPDX-License-Identifier: MPL-2.0
#![forbid(unsafe_code)]
#![warn(missing_docs)]

//! Constrained Delaunay triangulation with bounded quality refinement.
//!
//! # Why this exists
//!
//! The mesh boolean's retriangulation path was ear clipping. Ear clipping
//! terminates and it is fast, but it offers no angle guarantee: it fans from
//! whichever ear is convex, so a narrow opening near a far boundary corner
//! produces slivers whose aspect ratio is bounded only by the input's own
//! geometry. A sliver is not a cosmetic problem downstream -- normals of a
//! near-degenerate triangle are numerically meaningless, and every consumer
//! that shades, offsets, or measures from those normals inherits the error.
//!
//! This crate replaces "some valid triangulation" with a triangulation whose
//! quality is *stated and checked*:
//!
//! - every constraint edge survives as a union of output edges,
//! - the result is Delaunay away from the constraints (empty-circumcircle,
//!   decided by the certified `incircle` predicate rather than by a
//!   floating-point circumcircle test),
//! - with [`Quality`] refinement, interior angles meet a caller-chosen
//!   minimum, or the call reports that it could not get there.
//!
//! # The guarantee is conditional, and says so
//!
//! Ruppert refinement does not terminate for every input. Two boundary
//! segments meeting at a small angle cannot be fixed by inserting interior
//! points: splitting one segment to fix the angle creates a shorter segment
//! that is itself too close to its neighbour, and the process diverges. This
//! implementation therefore carries an explicit Steiner budget and reports
//! [`RefineOutcome::Capped`] when it stops early. The triangulation is still
//! valid and still constrained-Delaunay when capped -- only the angle bound
//! is unmet. Silently returning a worse mesh than requested would make the
//! quality parameter a lie.

use axiolid_core::Point2;
use axiolid_guarantees::{Certified, Sign};
use axiolid_predicates::{incircle, orient2d};

mod build;
mod mesh;
mod recover;
mod refine;

pub use build::triangulate;
pub use mesh::{Triangulation, TriangulationError};
pub use refine::{refine, triangulate_refined, Quality, RefineOutcome};

/// A constraint edge, as indices into the input point slice.
///
/// Held as a pair rather than as two points so the caller's vertex identity
/// survives the triangulation: a consumer that knows "edge 3 was my window
/// head" can still find it in the output.
#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord)]
pub struct Constraint {
    /// Index of the first endpoint.
    pub a: u32,
    /// Index of the second endpoint.
    pub b: u32,
}

impl Constraint {
    /// Construct a constraint between two vertex indices.
    ///
    /// Endpoints are stored in ascending order so that `(a, b)` and `(b, a)`
    /// are the same constraint. An undirected edge compared directionally is
    /// a classic source of "constraint lost" bugs that only appear when the
    /// caller happens to list an edge the other way round.
    #[must_use]
    pub fn new(a: u32, b: u32) -> Self {
        if a <= b {
            Self { a, b }
        } else {
            Self { a: b, b: a }
        }
    }
}

/// The proven sign of a predicate, or `Sign::Zero` when the configuration is
/// exactly degenerate.
///
/// The certified predicates escalate to exact arithmetic internally, so
/// `Uncertain` cannot survive a top-level call. Treating it as `Zero` rather
/// than unwrapping keeps this total: a degenerate answer is a real geometric
/// outcome here (three collinear points, four cocircular ones), and the
/// callers below all branch on strict positivity.
pub(crate) fn decided(certified: Certified) -> Sign {
    match certified {
        Certified::Certain { sign, .. } => sign,
        // `Certified` is `#[non_exhaustive]`, so a wildcard is required. Any
        // future variant means "not proven", which is the same conservative
        // answer as `Uncertain`.
        _ => Sign::Zero,
    }
}

/// Orientation of three points, decided exactly.
///
/// Wraps the certified predicate so the sign convention is stated once here
/// rather than re-derived at each call site.
pub(crate) fn turns_left(a: Point2, b: Point2, c: Point2) -> bool {
    decided(orient2d(a, b, c)) == Sign::Positive
}

/// Whether `a`, `b`, `c` are exactly collinear.
pub(crate) fn collinear(a: Point2, b: Point2, c: Point2) -> bool {
    decided(orient2d(a, b, c)) == Sign::Zero
}

/// Whether `d` lies strictly inside the circumcircle of `a`, `b`, `c`.
///
/// `a`, `b`, `c` must be counter-clockwise; the caller guarantees that by
/// construction. Decided by the certified `incircle` predicate: a
/// floating-point circumcircle test flips sign on nearly-cocircular input,
/// and a flip decision made on a wrong sign can cycle forever.
pub(crate) fn in_circumcircle(a: Point2, b: Point2, c: Point2, d: Point2) -> bool {
    decided(incircle(a, b, c, d)) == Sign::Positive
}
