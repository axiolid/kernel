//! Piecewise-linear curves.

use axiolid_core::{Point2, Point3};

/// Piecewise-linear curve preserving source vertex order.
///
/// Its native parameter is one unit per segment, not arc length and not a
/// normalised `(0, 1)`: an open n-point polyline spans `(0, n - 1)` and a
/// closed one `(0, n)`. A trim or profile-segment domain of `(0, 1)` on a
/// multi-segment polyline therefore selects only its first edge; the scalar
/// flattener in `axiolid-evaluate` refuses that rather than silently dropping
/// the remaining vertices.
#[derive(Debug, Clone, Default, PartialEq)]
pub struct Polyline<P> {
    /// Ordered control points.
    pub points: Vec<P>,
    /// Whether the final point connects back to the first.
    pub closed: bool,
}

/// Two-dimensional polyline.
pub type Polyline2 = Polyline<Point2>;
/// Three-dimensional polyline.
pub type Polyline3 = Polyline<Point3>;
