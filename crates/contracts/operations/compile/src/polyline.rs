//! Curves compiled to polylines, with a certified chord bound (#311).
//!
//! A curve node has no area, so it has no triangle mesh: compiling it with
//! [`crate::MeshCompiler::compile_mesh`] is refused. It is a separate
//! operation instead of a second output kind of `compile_mesh`, so the mesh
//! path's signature and every consumer of it stay as they are.

use axiolid_contracts::{Backend, ExecutionOptions, GeomResult};
use axiolid_core::{Point3, Scalar};
use axiolid_model::{GeometryGraph, NodeId};

/// How far the exact curve may lie from its polyline, reported the way a
/// compiled mesh's deviation is: the chord budget asked for, and a
/// certified upper bound on the distance from every point of the exact
/// curve to the nearest point of the polyline, both in world units.
///
/// The bound is one-sided, exact curve to polyline, and it is always
/// present: a compiler that cannot bound a curve refuses it rather than
/// returning a polyline with no bound.
#[non_exhaustive]
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct ChordDeviation {
    /// The chord budget the compilation was asked for:
    /// `ExecutionOptions::chord_error`, else the linear tolerance.
    pub requested: Scalar,
    /// An upper bound on the distance from every point of the exact curve
    /// to the polyline.
    pub bound: Scalar,
}

impl ChordDeviation {
    /// A bound against the budget it was asked for.
    #[must_use]
    pub const fn new(requested: Scalar, bound: Scalar) -> Self {
        Self { requested, bound }
    }

    /// Whether the polyline is certified within the requested budget.
    ///
    /// `false` means "not certified", not "missed": the bound is an upper
    /// bound, and the true deviation may still be below the budget.
    #[must_use]
    pub fn meets_requested(&self) -> bool {
        self.bound <= self.requested
    }
}

/// A curve compiled to a 3D polyline.
#[non_exhaustive]
#[derive(Debug, Clone, PartialEq)]
pub struct PolylineOutcome {
    /// The vertices in the curve's direction of travel. A closed curve
    /// ends on the vertex it starts on.
    pub points: Vec<Point3>,
    /// Whether the curve is closed: the last vertex is the first.
    pub closed: bool,
    /// How far the exact curve may lie from the polyline.
    pub deviation: ChordDeviation,
}

impl PolylineOutcome {
    /// A polyline and its deviation. `closed` is read from the points: at
    /// least three of them, the last equal to the first.
    #[must_use]
    pub fn new(points: Vec<Point3>, deviation: ChordDeviation) -> Self {
        let closed = points.len() > 2 && points.first() == points.last();
        Self {
            points,
            closed,
            deviation,
        }
    }
}

/// Backend that compiles a curve node of a graph to a 3D polyline.
///
/// Implementations must return [`axiolid_contracts::GeomError::Unsupported`]
/// or [`axiolid_contracts::GeomError::UnsupportedInput`] for a node they
/// cannot compile or cannot bound, naming it; they never return a polyline
/// whose [`ChordDeviation::bound`] is not an upper bound.
pub trait PolylineCompiler: Backend {
    /// Compile one root to a polyline within the options' chord budget.
    fn compile_polyline(
        &self,
        graph: &GeometryGraph,
        root: NodeId,
        options: &ExecutionOptions,
    ) -> GeomResult<PolylineOutcome>;
}
