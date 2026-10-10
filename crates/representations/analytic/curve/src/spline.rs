//! B-spline and rational B-spline curve data.

use axiolid_core::{Point2, Point3, Scalar};

/// How the knot vector was specified by the source representation.
#[non_exhaustive]
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
#[cfg_attr(
    feature = "serde",
    derive(serde::Serialize, serde::Deserialize),
    serde(deny_unknown_fields)
)]
pub enum KnotSpec {
    /// Uniform spacing, not necessarily clamped.
    Uniform,
    /// Quasi-uniform spacing with clamped ends.
    QuasiUniform,
    /// Piecewise Bezier knot multiplicities.
    PiecewiseBezier,
    /// Explicit knots and multiplicities.
    Unspecified,
}

/// Exact B-spline data. Rational curves carry one weight per control point.
#[derive(Debug, Clone, PartialEq)]
#[cfg_attr(
    feature = "serde",
    derive(serde::Serialize, serde::Deserialize),
    serde(deny_unknown_fields)
)]
pub struct BSplineCurve<P> {
    /// Polynomial degree.
    pub degree: u16,
    /// Ordered control points.
    pub control_points: Vec<P>,
    /// Distinct knot values.
    pub knots: Vec<Scalar>,
    /// Multiplicity for each distinct knot.
    pub multiplicities: Vec<u32>,
    /// Optional rational weights, one per control point.
    pub weights: Option<Vec<Scalar>>,
    /// Whether the source declares the curve closed.
    pub closed: bool,
    /// Whether the source declares self intersection.
    pub self_intersect: Option<bool>,
    /// Source knot convention.
    pub knot_spec: KnotSpec,
}

/// Two-dimensional B-spline curve.
pub type BSplineCurve2 = BSplineCurve<Point2>;
/// Three-dimensional B-spline curve.
pub type BSplineCurve3 = BSplineCurve<Point3>;
