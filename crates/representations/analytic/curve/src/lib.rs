#![forbid(unsafe_code)]

//! Exact, format-neutral curve representations and evaluation contracts.
//!
//! Composite, trimmed, offset, and surface-bound curves are graph relations in
//! `axiolid-model`; keeping them there avoids a curve/surface dependency cycle.

pub mod banked;
pub mod chain;
pub mod conic;
pub mod elevation;
pub mod evaluate;
pub mod implicit;
mod intrinsic;
mod intrinsic3;
pub mod linear;
pub mod pair_section;
pub mod quadric_section;
pub mod seam;
pub mod sinusoid;
pub mod spline;
pub mod spline_surface;
pub mod torus_section;

pub use banked::{
    bank_angle, BankConvention, BankError, Banked3, CantForm, CantLaw, CantPiece, CantValue,
};
pub use chain::{Chain2, ChainPiece2};
pub use conic::{Circle2, Circle3, Ellipse2, Ellipse3};
pub use elevation::{Elevated3, ElevationLaw};
pub use evaluate::CurveEvaluator;
pub use implicit::{
    Axis, Basis, Carrier, Field2, ImplicitCell, ImplicitCurve2, ImplicitSection3, Jet2,
    LiftedCurve2, PatchField2, SeriesField2, SurfaceJet,
};
pub use intrinsic::{CurvatureLaw, Harmonic, Intrinsic2};
pub use intrinsic3::Intrinsic3;
pub use linear::{Line, Line2, Line3, Polyline, Polyline2, Polyline3};
pub use pair_section::{PairNode, PairSection3};
pub use quadric_section::{Branch, QuadraticGraph2, RuledCarrier, RuledSection3, Trig2};
pub use seam::SeamSide;
pub use sinusoid::Sinusoid2;
pub use spline::{BSplineCurve, BSplineCurve2, BSplineCurve3, KnotSpec};
pub use spline_surface::BSplineSurface;
pub use torus_section::{AngleGraph2, TorusCarrier, TorusSection3};

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
    /// One root of `a(t) v^2 + b(t) v + c(t) = 0` as a graph over `t`: the
    /// exact pcurve of a quadric's cut across a cylinder or cone. See
    /// [`QuadraticGraph2`] (ADR 0076).
    QuadraticGraph(QuadraticGraph2),
    /// A solution of `a(t) cos u + b(t) sin u = c(t)` as `u` over `t`: the
    /// exact pcurve of a plane's or sphere's cut across a torus. See
    /// [`AngleGraph2`] (ADR 0076).
    AngleGraph(AngleGraph2),
    /// A stretch of a field's zero set, the field being another analytic
    /// surface's equation read in this surface's parameters: the exact
    /// pcurve of any section between planes, quadrics and tori. See
    /// [`ImplicitCurve2`] (ADR 0077).
    Implicit(ImplicitCurve2),
    /// A space curve read in an analytic surface's parameters, sharing the
    /// curve's parameter: the pcurve of a B-spline's section on the analytic
    /// face it meets. See [`LiftedCurve2`] (ADR 0077).
    Lifted(LiftedCurve2),
    /// Pieces placed rigidly end to end and parameterised by cumulative arc
    /// length: curvature-law pieces and parametric curves read by arc
    /// length, such as a cubic parabola between a line and an arc. Usable
    /// as the plan of an [`Elevated3`]. See [`Chain2`].
    Chain(Chain2),
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
    /// A quadric's cut across a cylinder or cone: the carrier evaluated along
    /// a [`QuadraticGraph2`]. See [`RuledSection3`] (ADR 0076).
    RuledSection(RuledSection3),
    /// A plane's or sphere's cut across a torus: the torus evaluated along
    /// an [`AngleGraph2`]. See [`TorusSection3`] (ADR 0076).
    TorusSection(TorusSection3),
    /// An [`ImplicitCurve2`] on its analytic carrier: the section of a
    /// torus by a cylinder, cone or torus, and any other analytic pair with
    /// no closed form. See [`ImplicitSection3`] (ADR 0077).
    ImplicitSection(ImplicitSection3),
    /// Where two B-spline surfaces meet, carried by nodes on both and
    /// defined between them by the surfaces themselves. See
    /// [`PairSection3`] (ADR 0077).
    PairSection(PairSection3),
    /// An elevated centreline carrying a cant law: the section rolls about
    /// the tangent by a bank angle, under a named convention. See
    /// [`Banked3`] (ADR 0081).
    Banked(Banked3),
}
