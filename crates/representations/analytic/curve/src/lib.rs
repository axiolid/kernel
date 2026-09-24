#![forbid(unsafe_code)]

//! Exact, format-neutral curve representations and evaluation contracts.
//!
//! Composite, trimmed, offset, and surface-bound curves are graph relations in
//! `axiolid-model`; keeping them there avoids a curve/surface dependency cycle.

pub mod conic;
pub mod elevation;
pub mod evaluate;
mod intrinsic;
mod intrinsic3;
pub mod linear;
pub mod sinusoid;
pub mod spline;

pub use conic::{Circle2, Circle3, Ellipse2, Ellipse3};
pub use elevation::{Elevated3, ElevationLaw};
pub use evaluate::CurveEvaluator;
pub use intrinsic::{CurvatureLaw, Harmonic, Intrinsic2};
pub use intrinsic3::Intrinsic3;
pub use linear::{Line, Line2, Line3, Polyline, Polyline2, Polyline3};
pub use sinusoid::Sinusoid2;
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
    /// The graph `v = mean + a cos(t) + b sin(t)`: the exact pcurve of a
    /// plane's cut across a cylinder, read in the cylinder's (angle, height)
    /// parameters. See [`Sinusoid2`] (ADR 0071).
    Sinusoid(Sinusoid2),
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
    /// Curve given by its natural equations: curvature AND torsion as
    /// functions of arc length, anchored to a start frame. Carries helices
    /// and general space spirals as exact values.
    Intrinsic(Intrinsic3),
    /// A planar curve paired with an elevation law: the exact composition an
    /// alignment centreline is authored as. See [`Elevated3`].
    Elevated(Elevated3),
}
