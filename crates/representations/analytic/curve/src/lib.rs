#![forbid(unsafe_code)]

//! Exact, format-neutral curve representations and evaluation contracts.
//!
//! Composite, trimmed, offset, and surface-bound curves are graph relations in
//! `axiolid-model`; keeping them there avoids a curve/surface dependency cycle.

pub mod conic;
pub mod evaluate;
pub mod intrinsic;
pub mod linear;
pub mod spline;

pub use conic::{Circle2, Circle3, Ellipse2, Ellipse3};
pub use evaluate::CurveEvaluator;
pub use intrinsic::{CurvatureLaw, Intrinsic2};
pub use linear::{Line, Line2, Line3, Polyline, Polyline2, Polyline3};
pub use spline::{BSplineCurve, BSplineCurve2, BSplineCurve3, KnotSpec};

/// The focused linear vocabulary, re-exported for consumers that want to name
/// its origin explicitly. A line-only consumer should depend on
/// `axiolid-linear` directly instead of paying for this aggregate.
pub mod linear_vocabulary {
    pub use axiolid_linear::*;
}

/// Atomic two-dimensional curve values.
#[non_exhaustive]
#[derive(Debug, Clone, PartialEq)]
pub enum Curve2 {
    /// Infinite line.
    Line(Line2),
    /// Circle.
    Circle(Circle2),
    /// Ellipse.
    Ellipse(Ellipse2),
    /// Piecewise linear curve.
    Polyline(Polyline2),
    /// Polynomial or rational B-spline.
    BSpline(BSplineCurve2),
    /// Curve given by its natural equation: curvature as a function of arc
    /// length, anchored to a start frame. Carries clothoid and other
    /// transition spirals exactly, which no parametric variant can.
    Intrinsic(Intrinsic2),
}

/// Atomic three-dimensional curve values.
#[non_exhaustive]
#[derive(Debug, Clone, PartialEq)]
pub enum Curve3 {
    /// Infinite line.
    Line(Line3),
    /// Circle in a plane.
    Circle(Circle3),
    /// Ellipse in a plane.
    Ellipse(Ellipse3),
    /// Piecewise linear curve.
    Polyline(Polyline3),
    /// Polynomial or rational B-spline.
    BSpline(BSplineCurve3),
}
