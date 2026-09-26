//! Scalar reference implementation of curve evaluation (ADR 0012).
//!
//! # What this closes
//!
//! `axiolid-curve` declares `Curve2`/`Curve3` and a `CurveEvaluator` trait. Until
//! now nothing in the workspace implemented that trait, so every declared curve
//! family was inert data. `axiolid-mesh-compile` worked around this with its own
//! private circle flattener and refused ellipses and B-splines outright.
//!
//! # Design
//!
//! Evaluation is analytic per family, never a generic subdivision fallback:
//!
//! - `Line`     -- `origin + t * direction`, exact.
//! - `Circle`   -- `origin + r*(cos t * x + sin t * y)`, `t` in radians.
//! - `Ellipse`  -- same with independent semi-axes. Note `t` is the
//!   *parametric* angle, not the polar angle; they differ except on axis.
//! - `Polyline` -- `t` in `[0, n)`, integer part selects the segment. Chosen
//!   over arc-length parameterization because it is exact and stable under
//!   degenerate (zero-length) segments, which imported data contains.
//! - `BSpline`  -- de Boor. Rational curves evaluate in homogeneous space and
//!   project, which is the only way to get correct rational derivatives.
//!
//! Derivatives are closed-form. A finite-difference derivative would make the
//! curvature oracle in `tests/curve.rs` self-referential: it would be checking
//! a difference quotient against a difference quotient.
//!
//! # Frames are used as given
//!
//! Imported frames may be non-orthonormal. Evaluation applies the frame axes as
//! written rather than orthonormalizing, so a caller sees the geometry its
//! source actually declared. Validation is a separate concern (`axiolid-heal`).

use axiolid_contracts::{GeomError, GeomResult};
use axiolid_core::{Frame2, Frame3, Interval, Point2, Point3, Scalar, Tolerance, Vec2, Vec3};
use axiolid_curve::{
    BSplineCurve, BSplineCurve2, BSplineCurve3, Circle2, Circle3, Curve2, Curve3, CurveEvaluator,
    Ellipse2, Ellipse3, Line2, Line3, Polyline2, Polyline3,
};

use crate::nurbs::SplineAxis;

/// Portable scalar curve evaluator.
///
/// Stateless: every method is a pure function of its arguments, so one instance
/// is freely shareable across threads.
#[derive(Debug, Clone, Copy, Default)]
pub struct ScalarCurve;

impl ScalarCurve {
    /// Construct the evaluator.
    #[must_use]
    pub const fn new() -> Self {
        Self
    }
}

/// Position and the first two parameter derivatives of a curve.
///
/// Keeping the derivatives with the point prevents callers from accidentally
/// mixing results evaluated at different parameters. All derivatives use the
/// curve's native parameter, not arc length.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct CurveJet<P, D> {
    /// Position at the requested parameter.
    pub point: P,
    /// First derivative with respect to the native parameter.
    pub first: D,
    /// Second derivative with respect to the native parameter.
    pub second: D,
}

// --- parameter domains ------------------------------------------------------

/// Domain of a 2D curve.
#[must_use]
pub fn domain2(curve: &Curve2) -> Interval {
    match curve {
        // A line is infinite; the unit interval is the conventional finite
        // window. Bounded use always arrives via `ProfileSegment::domain`.
        Curve2::Line(_) => Interval::UNIT,
        Curve2::Circle(_) | Curve2::Ellipse(_) => full_turn(),
        Curve2::Polyline(p) => polyline_domain(p.points.len(), p.closed),
        Curve2::BSpline(b) => spline_domain(b),
        // Parameterised by ARC LENGTH, not by a unit parameter: the domain is
        // the span the law is declared over. A non-finite or non-positive
        // length claims no domain rather than a guessed one.
        Curve2::Intrinsic(i) if i.length.is_finite() && i.length > 0.0 => Interval {
            start: 0.0,
            end: i.length,
        },
        // Periodic in the angle it is parameterised by, like a circle.
        Curve2::Sinusoid(_) => full_turn(),
        // Unknown family: no domain is knowable, so claim none.
        _ => Interval {
            start: 0.0,
            end: 0.0,
        },
    }
}

/// Domain of a 3D curve.
#[must_use]
pub fn domain3(curve: &Curve3) -> Interval {
    match curve {
        Curve3::Line(_) => Interval::UNIT,
        Curve3::Circle(_) | Curve3::Ellipse(_) => full_turn(),
        Curve3::Polyline(p) => polyline_domain(p.points.len(), p.closed),
        Curve3::BSpline(b) => spline_domain(b),
        // Parameterised by ARC LENGTH, not by a unit parameter: the domain is
        // the span the laws are declared over. A non-finite or non-positive
        // length claims no domain rather than a guessed one.
        Curve3::Intrinsic(i) if i.length.is_finite() && i.length > 0.0 => Interval {
            start: 0.0,
            end: i.length,
        },
        _ => Interval {
            start: 0.0,
            end: 0.0,
        },
    }
}

fn full_turn() -> Interval {
    Interval {
        start: 0.0,
        end: core::f64::consts::TAU,
    }
}

/// Polyline parameter runs `[0, segment_count]`.
fn polyline_domain(count: usize, closed: bool) -> Interval {
    let segments = if closed {
        count
    } else {
        count.saturating_sub(1)
    };
    Interval {
        start: 0.0,
        end: segments as Scalar,
    }
}

/// Domain of a validated B-spline axis. Invalid imported data reports an empty
/// domain through the infallible evaluator trait and is rejected by evaluation.
fn spline_domain<P>(b: &BSplineCurve<P>) -> Interval {
    SplineAxis::new(
        &b.knots,
        &b.multiplicities,
        b.degree,
        b.control_points.len(),
        "B-spline curve",
    )
    .map_or(
        Interval {
            start: 0.0,
            end: 0.0,
        },
        |axis| {
            let (start, end) = axis.domain();
            Interval { start, end }
        },
    )
}

// --- 2D evaluation ----------------------------------------------------------

/// Position on a 2D curve.
pub fn evaluate2(curve: &Curve2, t: Scalar) -> GeomResult<Point2> {
    finite(t)?;
    let value = match curve {
        Curve2::Line(l) => Ok(line_point(l.origin, l.direction, t)),
        Curve2::Circle(c) => Ok(conic_point2(&c.frame, c.radius, c.radius, t)),
        Curve2::Ellipse(e) => Ok(conic_point2(&e.frame, e.semi_axis_x, e.semi_axis_y, t)),
        Curve2::Polyline(p) => polyline_point(&p.points, p.closed, t),
        Curve2::BSpline(b) => de_boor(b, t, |p| [p.x, p.y], |c| Point2::new(c[0], c[1])),
        // Position has no elementary closed form for a general curvature law
        // (a clothoid needs Fresnel integrals), so it is quadrature over the
        // exact heading rather than a parametric formula.
        Curve2::Intrinsic(i) => crate::arc_length::intrinsic_point(i, t),
        // The parameter is the first coordinate; closed form, no sampling.
        Curve2::Sinusoid(w) => Ok(Point2::new(t, w.height(t))),
        // Likewise, one root of a quadratic in the height (ADR 0076).
        Curve2::QuadraticGraph(g) => g
            .height(t)
            .map(|v| Point2::new(t, v))
            .ok_or_else(|| outside_graph(t)),
        // The angle is the first coordinate, the parameter the second.
        Curve2::AngleGraph(g) => g
            .angle(t)
            .map(|u| Point2::new(u, t))
            .ok_or_else(|| outside_graph(t)),
        // `Curve*` is #[non_exhaustive]. An unknown family is refused by name
        // rather than approximated by whichever arm happens to be nearest.
        _ => Err(GeomError::Unsupported {
            backend: axiolid_contracts::BackendId::new("axiolid-reference"),
            operation: axiolid_contracts::Operation::CurveEvaluation,
        }),
    }?;
    finite2(value, "curve point")
}

/// First derivative of a 2D curve.
pub fn derivative2(curve: &Curve2, t: Scalar) -> GeomResult<Vec2> {
    finite(t)?;
    let value = match curve {
        Curve2::Line(l) => Ok(l.direction),
        Curve2::Circle(c) => Ok(conic_tangent2(&c.frame, c.radius, c.radius, t)),
        Curve2::Ellipse(e) => Ok(conic_tangent2(&e.frame, e.semi_axis_x, e.semi_axis_y, t)),
        Curve2::Polyline(p) => polyline_tangent(&p.points, p.closed, t),
        Curve2::BSpline(b) => de_boor_derivative(b, t, |p| [p.x, p.y], |c| Vec2::new(c[0], c[1])),
        // Arc-length parameterised, so the derivative is the UNIT tangent, and
        // the heading it is built from is exact -- only position needs
        // quadrature, never the tangent.
        Curve2::Intrinsic(i) => crate::arc_length::intrinsic_tangent(i, t),
        Curve2::Sinusoid(w) => {
            let (sin, cos) = t.sin_cos();
            Ok(Vec2::new(1.0, -w.cosine * sin + w.sine * cos))
        }
        Curve2::QuadraticGraph(g) => g
            .slope(t)
            .map(|slope| Vec2::new(1.0, slope))
            .ok_or_else(|| outside_graph(t)),
        Curve2::AngleGraph(g) => g
            .slope(t)
            .map(|slope| Vec2::new(slope, 1.0))
            .ok_or_else(|| outside_graph(t)),
        // `Curve*` is #[non_exhaustive]. An unknown family is refused by name
        // rather than approximated by whichever arm happens to be nearest.
        _ => Err(GeomError::Unsupported {
            backend: axiolid_contracts::BackendId::new("axiolid-reference"),
            operation: axiolid_contracts::Operation::CurveEvaluation,
        }),
    }?;
    finite2(value, "curve derivative")
}

/// Second derivative of a 2D curve with respect to its native parameter.
pub fn second_derivative2(curve: &Curve2, t: Scalar) -> GeomResult<Vec2> {
    finite(t)?;
    let value = match curve {
        Curve2::Line(_) | Curve2::Polyline(_) => Ok(Vec2::ZERO),
        Curve2::Circle(c) => Ok(conic_second2(&c.frame, c.radius, c.radius, t)),
        Curve2::Ellipse(e) => Ok(conic_second2(&e.frame, e.semi_axis_x, e.semi_axis_y, t)),
        Curve2::BSpline(b) => {
            de_boor_second_derivative(b, t, |p| [p.x, p.y], |c| Vec2::new(c[0], c[1]))
        }
        Curve2::Sinusoid(w) => {
            let (sin, cos) = t.sin_cos();
            Ok(Vec2::new(0.0, -w.cosine * cos - w.sine * sin))
        }
        Curve2::QuadraticGraph(g) => g
            .bend(t)
            .map(|bend| Vec2::new(0.0, bend))
            .ok_or_else(|| outside_graph(t)),
        Curve2::AngleGraph(g) => g
            .bend(t)
            .map(|bend| Vec2::new(bend, 0.0))
            .ok_or_else(|| outside_graph(t)),
        _ => Err(GeomError::Unsupported {
            backend: axiolid_contracts::BackendId::new("axiolid-reference"),
            operation: axiolid_contracts::Operation::CurveEvaluation,
        }),
    }?;
    finite2(value, "curve second derivative")
}

/// Second-order differential jet of a 2D curve.
pub fn bspline_jet2(curve: &BSplineCurve2, t: Scalar) -> GeomResult<CurveJet<Point2, Vec2>> {
    Ok(CurveJet {
        point: de_boor(curve, t, |p| [p.x, p.y], |c| Point2::new(c[0], c[1]))?,
        first: de_boor_derivative(curve, t, |p| [p.x, p.y], |c| Vec2::new(c[0], c[1]))?,
        second: de_boor_second_derivative(curve, t, |p| [p.x, p.y], |c| Vec2::new(c[0], c[1]))?,
    })
}

/// Second-order differential jet of a 2D curve.
pub fn jet2(curve: &Curve2, t: Scalar) -> GeomResult<CurveJet<Point2, Vec2>> {
    Ok(CurveJet {
        point: evaluate2(curve, t)?,
        first: derivative2(curve, t)?,
        second: second_derivative2(curve, t)?,
    })
}

/// A quadratic-graph parameter where its root does not exist, diverges, or
/// has a vertical tangent: outside the spans the curve was built for.
fn outside_graph(t: Scalar) -> GeomError {
    GeomError::Degenerate(format!(
        "quadratic section graph has no regular point at t = {t}"
    ))
}

// --- 3D evaluation ----------------------------------------------------------

/// Position on a 3D curve.
pub fn evaluate3(curve: &Curve3, t: Scalar) -> GeomResult<Point3> {
    finite(t)?;
    let value = match curve {
        Curve3::Line(l) => Ok(line_point(l.origin, l.direction, t)),
        Curve3::Circle(c) => Ok(conic_point3(&c.frame, c.radius, c.radius, t)),
        Curve3::Ellipse(e) => Ok(conic_point3(&e.frame, e.semi_axis_x, e.semi_axis_y, t)),
        Curve3::Polyline(p) => polyline_point(&p.points, p.closed, t),
        Curve3::BSpline(b) => de_boor(b, t, |p| [p.x, p.y, p.z], |c| Point3::new(c[0], c[1], c[2])),
        // Natural equations: `t` is ARC LENGTH, and the point comes from the
        // Frenet integrator (ADR 0061). Dispatching here is what lets the
        // graph's existing relation machinery -- trim, composite, sweep
        // directrix -- work on a torsion curve without special-casing it.
        Curve3::Intrinsic(i) => crate::frenet::frenet_point(i, t),
        // The carrier along its section graph (ADR 0076).
        Curve3::RuledSection(r) => r.point(t).ok_or_else(|| outside_graph(t)),
        Curve3::TorusSection(r) => r.point(t).ok_or_else(|| outside_graph(t)),
        // `Curve*` is #[non_exhaustive]. An unknown family is refused by name
        // rather than approximated by whichever arm happens to be nearest.
        _ => Err(GeomError::Unsupported {
            backend: axiolid_contracts::BackendId::new("axiolid-reference"),
            operation: axiolid_contracts::Operation::CurveEvaluation,
        }),
    }?;
    finite3(value, "curve point")
}

/// First derivative of a 3D curve.
pub fn derivative3(curve: &Curve3, t: Scalar) -> GeomResult<Vec3> {
    finite(t)?;
    let value = match curve {
        Curve3::Line(l) => Ok(l.direction),
        Curve3::Circle(c) => Ok(conic_tangent3(&c.frame, c.radius, c.radius, t)),
        Curve3::Ellipse(e) => Ok(conic_tangent3(&e.frame, e.semi_axis_x, e.semi_axis_y, t)),
        Curve3::Polyline(p) => polyline_tangent(&p.points, p.closed, t),
        Curve3::BSpline(b) => {
            de_boor_derivative(b, t, |p| [p.x, p.y, p.z], |c| Vec3::new(c[0], c[1], c[2]))
        }
        // Arc-length parameterised, so the derivative is the UNIT tangent.
        Curve3::Intrinsic(i) => crate::frenet::frenet_tangent(i, t),
        Curve3::RuledSection(r) => r.tangent(t).ok_or_else(|| outside_graph(t)),
        Curve3::TorusSection(r) => r.tangent(t).ok_or_else(|| outside_graph(t)),
        // `Curve*` is #[non_exhaustive]. An unknown family is refused by name
        // rather than approximated by whichever arm happens to be nearest.
        _ => Err(GeomError::Unsupported {
            backend: axiolid_contracts::BackendId::new("axiolid-reference"),
            operation: axiolid_contracts::Operation::CurveEvaluation,
        }),
    }?;
    finite3(value, "curve derivative")
}

/// Second derivative of a 3D curve with respect to its native parameter.
pub fn second_derivative3(curve: &Curve3, t: Scalar) -> GeomResult<Vec3> {
    finite(t)?;
    let value = match curve {
        Curve3::Line(_) | Curve3::Polyline(_) => Ok(Vec3::ZERO),
        Curve3::Circle(c) => Ok(conic_second3(&c.frame, c.radius, c.radius, t)),
        Curve3::Ellipse(e) => Ok(conic_second3(&e.frame, e.semi_axis_x, e.semi_axis_y, t)),
        Curve3::BSpline(b) => {
            de_boor_second_derivative(b, t, |p| [p.x, p.y, p.z], |c| Vec3::new(c[0], c[1], c[2]))
        }
        Curve3::RuledSection(r) => r.bend(t).ok_or_else(|| outside_graph(t)),
        Curve3::TorusSection(r) => r.bend(t).ok_or_else(|| outside_graph(t)),
        _ => Err(GeomError::Unsupported {
            backend: axiolid_contracts::BackendId::new("axiolid-reference"),
            operation: axiolid_contracts::Operation::CurveEvaluation,
        }),
    }?;
    finite3(value, "curve second derivative")
}

/// Second-order differential jet of a 3D curve.
pub fn bspline_jet3(curve: &BSplineCurve3, t: Scalar) -> GeomResult<CurveJet<Point3, Vec3>> {
    Ok(CurveJet {
        point: de_boor(
            curve,
            t,
            |p| [p.x, p.y, p.z],
            |c| Point3::new(c[0], c[1], c[2]),
        )?,
        first: de_boor_derivative(
            curve,
            t,
            |p| [p.x, p.y, p.z],
            |c| Vec3::new(c[0], c[1], c[2]),
        )?,
        second: de_boor_second_derivative(
            curve,
            t,
            |p| [p.x, p.y, p.z],
            |c| Vec3::new(c[0], c[1], c[2]),
        )?,
    })
}

/// Second-order differential jet of a 3D curve.
pub fn jet3(curve: &Curve3, t: Scalar) -> GeomResult<CurveJet<Point3, Vec3>> {
    Ok(CurveJet {
        point: evaluate3(curve, t)?,
        first: derivative3(curve, t)?,
        second: second_derivative3(curve, t)?,
    })
}

// --- family kernels ---------------------------------------------------------

fn finite(t: Scalar) -> GeomResult<()> {
    if t.is_finite() {
        Ok(())
    } else {
        Err(GeomError::InvalidInput(format!(
            "curve parameter must be finite, got {t}"
        )))
    }
}

fn finite2(value: Vec2, what: &str) -> GeomResult<Vec2> {
    if value.is_finite() {
        Ok(value)
    } else {
        Err(GeomError::Degenerate(format!("{what} is non-finite")))
    }
}

fn finite3(value: Vec3, what: &str) -> GeomResult<Vec3> {
    if value.is_finite() {
        Ok(value)
    } else {
        Err(GeomError::Degenerate(format!("{what} is non-finite")))
    }
}

fn line_point<P>(origin: P, direction: P, t: Scalar) -> P
where
    P: core::ops::Add<Output = P> + core::ops::Mul<Scalar, Output = P>,
{
    origin + direction * t
}

fn conic_point2(frame: &Frame2, rx: Scalar, ry: Scalar, t: Scalar) -> Point2 {
    frame.origin + frame.x * (rx * t.cos()) + frame.y * (ry * t.sin())
}

fn conic_tangent2(frame: &Frame2, rx: Scalar, ry: Scalar, t: Scalar) -> Vec2 {
    frame.x * (-rx * t.sin()) + frame.y * (ry * t.cos())
}

fn conic_second2(frame: &Frame2, rx: Scalar, ry: Scalar, t: Scalar) -> Vec2 {
    frame.x * (-rx * t.cos()) + frame.y * (-ry * t.sin())
}

fn conic_point3(frame: &Frame3, rx: Scalar, ry: Scalar, t: Scalar) -> Point3 {
    frame.origin + frame.x * (rx * t.cos()) + frame.y * (ry * t.sin())
}

fn conic_tangent3(frame: &Frame3, rx: Scalar, ry: Scalar, t: Scalar) -> Vec3 {
    frame.x * (-rx * t.sin()) + frame.y * (ry * t.cos())
}

fn conic_second3(frame: &Frame3, rx: Scalar, ry: Scalar, t: Scalar) -> Vec3 {
    frame.x * (-rx * t.cos()) + frame.y * (-ry * t.sin())
}

/// Segment index and local fraction for a polyline parameter.
///
/// Returns `None` when the polyline cannot be evaluated at all.
fn polyline_span(count: usize, closed: bool, t: Scalar) -> Option<(usize, usize, Scalar)> {
    let segments = if closed {
        count
    } else {
        count.saturating_sub(1)
    };
    if count < 2 || segments == 0 {
        return None;
    }
    // Clamp into range: the endpoint t == segments is the final vertex, which
    // would otherwise index one past the last segment.
    let clamped = t.clamp(0.0, segments as Scalar);
    let mut index = clamped.floor() as usize;
    if index >= segments {
        index = segments - 1;
    }
    let local = clamped - index as Scalar;
    let next = (index + 1) % count;
    Some((index, next, local))
}

fn polyline_point<P>(points: &[P], closed: bool, t: Scalar) -> GeomResult<P>
where
    P: Copy
        + core::ops::Add<Output = P>
        + core::ops::Sub<Output = P>
        + core::ops::Mul<Scalar, Output = P>,
{
    let (i, j, local) = polyline_span(points.len(), closed, t).ok_or_else(|| {
        GeomError::Degenerate(format!(
            "polyline with {} points has no evaluable segment",
            points.len()
        ))
    })?;
    Ok(points[i] + (points[j] - points[i]) * local)
}

fn polyline_tangent<P>(points: &[P], closed: bool, t: Scalar) -> GeomResult<P>
where
    P: Copy + core::ops::Sub<Output = P>,
{
    let (i, j, _) = polyline_span(points.len(), closed, t).ok_or_else(|| {
        GeomError::Degenerate(format!(
            "polyline with {} points has no evaluable segment",
            points.len()
        ))
    })?;
    // Derivative w.r.t. the unit-per-segment parameter is the full edge vector.
    Ok(points[j] - points[i])
}

// --- reusable de Boor core (shared with `crate::surface`) -------------------

/// Locate the knot span for `u` in a validated flat knot vector.
///
/// Extracted from [`spline_span`] so a tensor-product surface can reuse the
/// exact same span logic per axis. `n` is the control-point count, `d` the
/// degree; the caller has already checked `knots.len() == n + d + 1`.
pub(crate) fn span_in(knots: &[Scalar], n: usize, d: usize, u: Scalar) -> usize {
    let mut span = d;
    for (k, knot) in knots.iter().enumerate().take(n).skip(d) {
        if *knot <= u {
            span = k;
        } else {
            break;
        }
    }
    span
}

/// One de Boor recurrence over homogeneous coordinates.
///
/// `points` holds the `d+1` premultiplied control points influencing `span`,
/// `weights` their weights. Both are consumed in place. This is the numerical
/// heart shared by curve and surface evaluation: keeping one copy means a fix
/// to the recurrence cannot land in one and not the other.
pub(crate) fn de_boor_recurrence<const N: usize>(
    knots: &[Scalar],
    span: usize,
    d: usize,
    u: Scalar,
    points: &mut [[Scalar; N]],
    weights: &mut [Scalar],
) {
    for r in 1..=d {
        for j in (r..=d).rev() {
            let i = span - d + j;
            let denom = knots[i + d + 1 - r] - knots[i];
            let alpha = if denom.abs() > 0.0 {
                (u - knots[i]) / denom
            } else {
                0.0
            };
            for k in 0..N {
                points[j][k] = points[j - 1][k] * (1.0 - alpha) + points[j][k] * alpha;
            }
            weights[j] = weights[j - 1] * (1.0 - alpha) + weights[j] * alpha;
        }
    }
}

// --- de Boor ----------------------------------------------------------------

/// Shared setup: validated flat knots, degree, and the knot span for `t`.
fn spline_span<P>(b: &BSplineCurve<P>, t: Scalar) -> GeomResult<(Vec<Scalar>, usize, usize)> {
    let axis = SplineAxis::new(
        &b.knots,
        &b.multiplicities,
        b.degree,
        b.control_points.len(),
        "B-spline curve",
    )?;
    if let Some(weights) = &b.weights {
        if weights.len() != b.control_points.len() {
            return Err(GeomError::InvalidInput(format!(
                "B-spline has {} weights for {} control points",
                weights.len(),
                b.control_points.len()
            )));
        }
        if weights
            .iter()
            .any(|weight| !weight.is_finite() || *weight <= 0.0)
        {
            return Err(GeomError::InvalidInput(
                "B-spline weights must be finite and strictly positive".to_owned(),
            ));
        }
    }
    let t = axis.clamp(t);
    let span = span_in(&axis.knots, axis.count, axis.degree, t);
    Ok((axis.knots, span, axis.degree))
}

/// Convert and validate every control point before selecting a knot span.
/// Imported NaN/Inf coordinates must not be hidden in currently uninfluential
/// spans and surface later when the parameter changes.
fn finite_control_points<P, const N: usize, F>(
    control_points: &[P],
    to: &F,
) -> GeomResult<Vec<[Scalar; N]>>
where
    F: Fn(&P) -> [Scalar; N],
{
    let points: Vec<[Scalar; N]> = control_points.iter().map(to).collect();
    if points
        .iter()
        .flatten()
        .any(|coordinate| !coordinate.is_finite())
    {
        return Err(GeomError::InvalidInput(
            "B-spline control points must be finite".to_owned(),
        ));
    }
    Ok(points)
}

/// Position via de Boor's algorithm.
///
/// `to` and `from` convert between the point type and a fixed-size coordinate
/// array so 2D and 3D share one implementation. Rational curves are evaluated
/// in homogeneous coordinates `(w*x, w*y, [w*z], w)` and projected at the end.
fn de_boor<P, const N: usize, F, G, Q>(
    b: &BSplineCurve<P>,
    t: Scalar,
    to: F,
    from: G,
) -> GeomResult<Q>
where
    F: Fn(&P) -> [Scalar; N],
    G: Fn([Scalar; N]) -> Q,
{
    let (knots, span, d) = spline_span(b, t)?;
    let control_points = finite_control_points(&b.control_points, &to)?;
    let u = t.clamp(knots[d], knots[b.control_points.len()]);

    // Working set: the d+1 control points influencing this span, in homogeneous
    // form. The trailing slot holds the weight (1.0 for polynomial curves).
    let mut work: Vec<[Scalar; N]> = Vec::with_capacity(d + 1);
    let mut weights: Vec<Scalar> = Vec::with_capacity(d + 1);
    for j in 0..=d {
        let idx = span - d + j;
        let w = b.weights.as_ref().map_or(1.0, |ws| ws[idx]);
        let c = control_points[idx];
        // Premultiply by w: interpolating in homogeneous space is what makes
        // rational curves correct. Projecting first would be plain averaging.
        let homogeneous = core::array::from_fn(|k| c[k] * w);
        if homogeneous.iter().any(|value| !value.is_finite()) {
            return Err(GeomError::Degenerate(
                "B-spline homogeneous control point overflowed".to_owned(),
            ));
        }
        work.push(homogeneous);
        weights.push(w);
    }

    // A repeated knot makes an interval empty; the shared recurrence treats
    // that as alpha = 0, which is the correct limit.
    de_boor_recurrence(&knots, span, d, u, &mut work, &mut weights);

    let w = weights[d];
    if !w.is_finite() || w == 0.0 {
        return Err(GeomError::Degenerate(
            "B-spline weight collapsed to zero".to_owned(),
        ));
    }
    Ok(from(core::array::from_fn(|k| work[d][k] / w)))
}

/// First derivative via the hodograph, with the quotient rule for rationals.
///
/// The derivative of a degree-`d` B-spline is a degree-`(d-1)` B-spline over
/// the same knots minus their outermost entries, with control points
/// `d * (P[i+1] - P[i]) / (knots[i+d+1] - knots[i+1])`.
///
/// For a rational curve `C = A/w`, both `A` and `w` are differentiated in
/// homogeneous space and combined as `(A' - C * w') / w`.
fn de_boor_derivative<P, const N: usize, F, G, Q>(
    b: &BSplineCurve<P>,
    t: Scalar,
    to: F,
    from: G,
) -> GeomResult<Q>
where
    F: Fn(&P) -> [Scalar; N],
    G: Fn([Scalar; N]) -> Q,
{
    let (knots, _, d) = spline_span(b, t)?;
    let control_points = finite_control_points(&b.control_points, &to)?;
    let n = b.control_points.len();
    let u = t.clamp(knots[d], knots[n]);

    // Homogeneous control points, weight in a parallel array.
    let hom: Vec<[Scalar; N]> = (0..n)
        .map(|i| {
            let w = b.weights.as_ref().map_or(1.0, |ws| ws[i]);
            let c = control_points[i];
            core::array::from_fn(|k| c[k] * w)
        })
        .collect();
    if hom.iter().flatten().any(|value| !value.is_finite()) {
        return Err(GeomError::Degenerate(
            "B-spline homogeneous control point overflowed".to_owned(),
        ));
    }
    let hw: Vec<Scalar> = (0..n)
        .map(|i| b.weights.as_ref().map_or(1.0, |ws| ws[i]))
        .collect();

    // Hodograph control points.
    let mut dhom: Vec<[Scalar; N]> = Vec::with_capacity(n - 1);
    let mut dhw: Vec<Scalar> = Vec::with_capacity(n - 1);
    for i in 0..n - 1 {
        let denom = knots[i + d + 1] - knots[i + 1];
        let f = if denom.abs() > 0.0 {
            d as Scalar / denom
        } else {
            0.0
        };
        dhom.push(core::array::from_fn(|k| (hom[i + 1][k] - hom[i][k]) * f));
        dhw.push((hw[i + 1] - hw[i]) * f);
    }

    // Evaluate the hodograph at u with degree d-1 over the trimmed knots.
    let dknots = &knots[1..knots.len() - 1];
    let (da, dw) = eval_homogeneous(dknots, d - 1, &dhom, &dhw, u);
    // Evaluate the curve itself for the quotient rule.
    let (a, w) = eval_homogeneous(&knots, d, &hom, &hw, u);

    if !w.is_finite() || w == 0.0 {
        return Err(GeomError::Degenerate(
            "B-spline weight collapsed to zero".to_owned(),
        ));
    }
    // C = A/w  =>  C' = (A' - (A/w) * w') / w
    Ok(from(core::array::from_fn(|k| {
        (da[k] - (a[k] / w) * dw) / w
    })))
}

/// Second derivative via two homogeneous hodograph constructions.
///
/// For `C = A / w`, the rational recurrence is
/// `C'' = (A'' - 2 w' C' - w'' C) / w`.
fn de_boor_second_derivative<P, const N: usize, F, G, Q>(
    b: &BSplineCurve<P>,
    t: Scalar,
    to: F,
    from: G,
) -> GeomResult<Q>
where
    F: Fn(&P) -> [Scalar; N],
    G: Fn([Scalar; N]) -> Q,
{
    let (knots, _, degree) = spline_span(b, t)?;
    let control_points = finite_control_points(&b.control_points, &to)?;
    let count = b.control_points.len();
    let u = t.clamp(knots[degree], knots[count]);

    let points: Vec<[Scalar; N]> = (0..count)
        .map(|i| {
            let weight = b.weights.as_ref().map_or(1.0, |weights| weights[i]);
            core::array::from_fn(|axis| control_points[i][axis] * weight)
        })
        .collect();
    if points.iter().flatten().any(|value| !value.is_finite()) {
        return Err(GeomError::Degenerate(
            "B-spline homogeneous control point overflowed".to_owned(),
        ));
    }
    let weights: Vec<Scalar> = (0..count)
        .map(|i| b.weights.as_ref().map_or(1.0, |values| values[i]))
        .collect();

    let (point, weight) = eval_homogeneous(&knots, degree, &points, &weights, u);
    if !weight.is_finite() || weight == 0.0 {
        return Err(GeomError::Degenerate(
            "B-spline weight collapsed to zero".to_owned(),
        ));
    }

    let (first_points, first_weights) = derivative_controls(&points, &weights, &knots, degree);
    let first_knots = &knots[1..knots.len() - 1];
    let (first, first_weight) =
        eval_homogeneous(first_knots, degree - 1, &first_points, &first_weights, u);
    let position: [Scalar; N] = core::array::from_fn(|axis| point[axis] / weight);
    let first_projected: [Scalar; N] =
        core::array::from_fn(|axis| (first[axis] - position[axis] * first_weight) / weight);

    let (second, second_weight) = if degree == 1 {
        ([0.0; N], 0.0)
    } else {
        let (second_points, second_weights) =
            derivative_controls(&first_points, &first_weights, first_knots, degree - 1);
        let second_knots = &first_knots[1..first_knots.len() - 1];
        eval_homogeneous(second_knots, degree - 2, &second_points, &second_weights, u)
    };
    Ok(from(core::array::from_fn(|axis| {
        (second[axis] - 2.0 * first_weight * first_projected[axis] - second_weight * position[axis])
            / weight
    })))
}

/// Derivative control polygon for one homogeneous B-spline axis.
fn derivative_controls<const N: usize>(
    points: &[[Scalar; N]],
    weights: &[Scalar],
    knots: &[Scalar],
    degree: usize,
) -> (Vec<[Scalar; N]>, Vec<Scalar>) {
    let mut derivative_points = Vec::with_capacity(points.len() - 1);
    let mut derivative_weights = Vec::with_capacity(weights.len() - 1);
    for i in 0..points.len() - 1 {
        let denominator = knots[i + degree + 1] - knots[i + 1];
        let factor = if denominator.abs() > 0.0 {
            degree as Scalar / denominator
        } else {
            0.0
        };
        derivative_points.push(core::array::from_fn(|axis| {
            (points[i + 1][axis] - points[i][axis]) * factor
        }));
        derivative_weights.push((weights[i + 1] - weights[i]) * factor);
    }
    (derivative_points, derivative_weights)
}

/// de Boor over explicit homogeneous arrays; returns `(numerator, weight)`.
pub(crate) fn eval_homogeneous<const N: usize>(
    knots: &[Scalar],
    d: usize,
    hom: &[[Scalar; N]],
    hw: &[Scalar],
    u: Scalar,
) -> ([Scalar; N], Scalar) {
    let n = hom.len();
    if d == 0 {
        // Degree zero: piecewise constant, pick the containing span.
        let mut idx = 0;
        for (k, knot) in knots.iter().enumerate().take(n) {
            if *knot <= u {
                idx = k;
            }
        }
        return (hom[idx.min(n - 1)], hw[idx.min(n - 1)]);
    }
    let mut span = d;
    for (k, knot) in knots.iter().enumerate().take(n).skip(d) {
        if *knot <= u {
            span = k;
        } else {
            break;
        }
    }
    let mut work: Vec<[Scalar; N]> = (0..=d).map(|j| hom[span - d + j]).collect();
    let mut weights: Vec<Scalar> = (0..=d).map(|j| hw[span - d + j]).collect();
    for r in 1..=d {
        for j in (r..=d).rev() {
            let i = span - d + j;
            let denom = knots[i + d + 1 - r] - knots[i];
            let alpha = if denom.abs() > 0.0 {
                (u - knots[i]) / denom
            } else {
                0.0
            };
            for k in 0..N {
                work[j][k] = work[j - 1][k] * (1.0 - alpha) + work[j][k] * alpha;
            }
            weights[j] = weights[j - 1] * (1.0 - alpha) + weights[j] * alpha;
        }
    }
    (work[d], weights[d])
}

// --- adaptive flattening ----------------------------------------------------

/// Flatten a 2D curve over `domain` so the chord never deviates from the true
/// curve by more than `chord_tolerance`.
///
/// # Why bisection rather than a closed-form segment count
///
/// A count derived from radius and tolerance only works for circles. Bisecting
/// on measured sagitta works for every family, including rational splines whose
/// curvature varies along the span, and it degrades gracefully on the
/// degenerate inputs imported data actually contains.
///
/// The returned polyline includes both endpoints and is ordered along
/// increasing parameter. `max_depth` bounds the work: a caller gets a
/// deterministic result rather than an unbounded subdivision on a pathological
/// curve.
pub fn flatten2(
    curve: &Curve2,
    domain: Interval,
    chord_tolerance: Scalar,
    max_depth: u32,
) -> GeomResult<Vec<Point2>> {
    // A depth bound alone is not a resource bound: depth `d` permits `2^d`
    // segments. Cap the total point count too, so a curve that cannot meet
    // the tolerance fails fast instead of exhausting memory.
    const MAX_POINTS: usize = 1 << 16;
    if !(chord_tolerance.is_finite()
        && chord_tolerance.is_sign_positive()
        && chord_tolerance != 0.0)
    {
        return Err(GeomError::InvalidInput(format!(
            "chord tolerance must be positive and finite, got {chord_tolerance}"
        )));
    }
    // A line and a polyline are already exact between their breakpoints:
    // subdividing them adds vertices that carry no information.
    if let Curve2::Line(_) = curve {
        return Ok(vec![
            evaluate2(curve, domain.start)?,
            evaluate2(curve, domain.end)?,
        ]);
    }
    if let Curve2::Polyline(p) = curve {
        // A polyline's parameter is one unit per segment, so a caller passing
        // a normalized `(0, 1)` domain would silently collapse an n-vertex
        // ring to its first edge. That is data loss disguised as success, so
        // it is refused: a domain narrower than one segment can only be
        // intentional for a genuinely 1-segment polyline.
        let natural = polyline_domain(p.points.len(), p.closed);
        let requested = (domain.end - domain.start).abs();
        if natural.end > 1.0 && requested <= 1.0 {
            return Err(GeomError::InvalidInput(format!(
                "polyline domain {:?} spans {requested} of {} segments; a \
                 polyline parameter is one unit per segment, so this would \
                 discard {} vertices",
                domain,
                natural.end,
                p.points.len().saturating_sub(2)
            )));
        }
        return polyline_flatten(&p.points, p.closed, domain, |t| evaluate2(curve, t));
    }

    let mut out = vec![evaluate2(curve, domain.start)?];
    let eval = |t| evaluate2(curve, t);
    subdivide(
        &eval,
        domain.start,
        domain.end,
        chord_tolerance,
        max_depth.min(MAX_DEPTH_CEILING),
        MAX_POINTS,
        &mut out,
    )?;
    out.push(evaluate2(curve, domain.end)?);
    Ok(out)
}

/// Hard ceiling on recursion depth regardless of what a caller asks for.
///
/// 2^20 segments is already far past any usable tolerance; beyond this a
/// request is a bug, not a quality setting.
const MAX_DEPTH_CEILING: u32 = 20;

/// Emit interior points of `(a, b)` that are needed to meet the tolerance.
///
/// `budget` bounds total emitted points. Exceeding it is an error rather than
/// a truncation: silently returning a coarser polyline than the caller asked
/// for would violate the tolerance contract this function exists to honour.
/// The vector operations chord subdivision needs, in either dimension.
///
/// `Point2` and `Point3` are both glam vectors with the same surface, and
/// the subdivision below is genuinely dimension-independent: it measures a
/// sagitta and bisects a parameter interval, neither of which mentions a
/// coordinate count. This trait states that once instead of maintaining
/// two copies that can drift apart.
trait ChordPoint: Copy {
    fn sub(self, other: Self) -> Self;
    fn add_scaled(self, direction: Self, scale: Scalar) -> Self;
    fn dot(self, other: Self) -> Scalar;
    fn length(self) -> Scalar;
    fn length_squared(self) -> Scalar;
}

impl ChordPoint for Point2 {
    fn sub(self, other: Self) -> Self {
        self - other
    }
    fn add_scaled(self, direction: Self, scale: Scalar) -> Self {
        self + direction * scale
    }
    fn dot(self, other: Self) -> Scalar {
        Point2::dot(self, other)
    }
    fn length(self) -> Scalar {
        Point2::length(self)
    }
    fn length_squared(self) -> Scalar {
        Point2::length_squared(self)
    }
}

impl ChordPoint for Point3 {
    fn sub(self, other: Self) -> Self {
        self - other
    }
    fn add_scaled(self, direction: Self, scale: Scalar) -> Self {
        self + direction * scale
    }
    fn dot(self, other: Self) -> Scalar {
        Point3::dot(self, other)
    }
    fn length(self) -> Scalar {
        Point3::length(self)
    }
    fn length_squared(self) -> Scalar {
        Point3::length_squared(self)
    }
}

/// Perpendicular distance from `m` to the chord `a`-`b`.
fn sagitta<P: ChordPoint>(a: P, b: P, m: P) -> Scalar {
    let ab = b.sub(a);
    let len2 = ab.length_squared();
    if len2 <= 0.0 {
        // Degenerate chord: fall back to point distance so a closed curve
        // whose endpoints coincide still subdivides.
        return m.sub(a).length();
    }
    let t = (m.sub(a).dot(ab) / len2).clamp(0.0, 1.0);
    m.sub(a.add_scaled(ab, t)).length()
}

/// Emit the interior points of `(a, b)` needed to meet `tol`.
///
/// Shared by both dimensions; the caller supplies the evaluator. Depth
/// exhaustion is an error rather than a truncation, matching the 2D
/// contract: silently returning a coarser polyline than asked for would
/// break the tolerance guarantee the caller is relying on.
fn subdivide<P, F>(
    eval: &F,
    a: Scalar,
    b: Scalar,
    tol: Scalar,
    depth: u32,
    budget: usize,
    out: &mut Vec<P>,
) -> GeomResult<()>
where
    P: ChordPoint,
    F: Fn(Scalar) -> GeomResult<P>,
{
    if out.len() >= budget {
        return Err(GeomError::Degenerate(format!(
            "curve flattening exceeded {budget} points before meeting the \
             chord tolerance {tol}; the curve may be degenerate"
        )));
    }
    let mid = 0.5 * (a + b);
    let pa = eval(a)?;
    let pb = eval(b)?;
    // A parameter interval too small to bisect cannot be refined further:
    // `mid` equals `a` or `b` in floating point. Returning the chord anyway
    // would hand back an unverified approximation, so this fails closed --
    // unless the chord already meets the tolerance, in which case there was
    // nothing left to verify.
    if !(mid > a && mid < b) {
        if sagitta(pa, pb, pa) <= tol && (pb.sub(pa)).length() <= tol {
            return Ok(());
        }
        return Err(GeomError::Degenerate(format!(
            "curve parameter interval ({a}, {b}) is too small to bisect but \
             its chord still exceeds the tolerance {tol}"
        )));
    }
    let pm = eval(mid)?;
    if sagitta(pa, pb, pm) <= tol {
        // Within tolerance: the chord a->b stands, no interior point.
        return Ok(());
    }
    if depth == 0 {
        return Err(GeomError::BudgetExceeded {
            resource: "curve flattening depth",
        });
    }
    subdivide(eval, a, mid, tol, depth - 1, budget, out)?;
    out.push(pm);
    subdivide(eval, mid, b, tol, depth - 1, budget, out)?;
    Ok(())
}

/// Polylines flatten to their own breakpoints, restricted to `domain`.
fn polyline_flatten<P, F>(
    points: &[P],
    closed: bool,
    domain: Interval,
    eval: F,
) -> GeomResult<Vec<P>>
where
    P: Copy,
    F: Fn(Scalar) -> GeomResult<P>,
{
    let segments = if closed {
        points.len()
    } else {
        points.len().saturating_sub(1)
    };
    if segments == 0 {
        return Err(GeomError::Degenerate(
            "polyline has no evaluable segment".to_owned(),
        ));
    }
    let lo = domain.start.min(domain.end);
    let hi = domain.start.max(domain.end);
    let mut out = vec![eval(lo)?];
    // Interior breakpoints are the integer parameters strictly inside.
    let first = lo.floor() as i64 + 1;
    let last = hi.ceil() as i64 - 1;
    for k in first..=last {
        let t = k as Scalar;
        if t > lo && t < hi {
            out.push(eval(t)?);
        }
    }
    out.push(eval(hi)?);
    Ok(out)
}

// --- trait wiring -----------------------------------------------------------

impl CurveEvaluator<Curve2> for ScalarCurve {
    type Point = Point2;
    type Derivative = Vec2;
    type Error = GeomError;

    fn domain(&self, curve: &Curve2) -> Interval {
        domain2(curve)
    }

    fn evaluate(
        &self,
        curve: &Curve2,
        t: Scalar,
        _tolerance: Tolerance,
    ) -> Result<Self::Point, Self::Error> {
        evaluate2(curve, t)
    }

    fn derivative(
        &self,
        curve: &Curve2,
        t: Scalar,
        _tolerance: Tolerance,
    ) -> Result<Self::Derivative, Self::Error> {
        derivative2(curve, t)
    }
}

impl CurveEvaluator<Curve3> for ScalarCurve {
    type Point = Point3;
    type Derivative = Vec3;
    type Error = GeomError;

    fn domain(&self, curve: &Curve3) -> Interval {
        domain3(curve)
    }

    fn evaluate(
        &self,
        curve: &Curve3,
        t: Scalar,
        _tolerance: Tolerance,
    ) -> Result<Self::Point, Self::Error> {
        evaluate3(curve, t)
    }

    fn derivative(
        &self,
        curve: &Curve3,
        t: Scalar,
        _tolerance: Tolerance,
    ) -> Result<Self::Derivative, Self::Error> {
        derivative3(curve, t)
    }
}

// Silence unused-import warnings for types only named in signatures.
#[allow(unused)]
fn _type_anchors(_: Circle2, _: Circle3, _: Ellipse2, _: Ellipse3, _: Line2, _: Line3) {}
#[allow(unused)]
fn _poly_anchors(_: Polyline2, _: Polyline3) {}

/// Flatten a 3D curve to a polyline within `chord_tolerance`.
///
/// The 3D twin of [`flatten2`], sharing its subdivision, its resource
/// bounds and its polyline contract. Sampling is adaptive: a curve is
/// bisected only where the chord actually departs from it, so a gentle
/// arc costs few points and a tight one costs many, and neither is
/// decided by a fixed count chosen in advance.
pub fn flatten3(
    curve: &Curve3,
    domain: Interval,
    chord_tolerance: Scalar,
    max_depth: u32,
) -> GeomResult<Vec<Point3>> {
    // A depth bound alone is not a resource bound: depth `d` permits `2^d`
    // segments. Cap the total point count too, so a curve that cannot meet
    // the tolerance fails fast instead of exhausting memory.
    const MAX_POINTS: usize = 1 << 16;
    if !(chord_tolerance.is_finite()
        && chord_tolerance.is_sign_positive()
        && chord_tolerance != 0.0)
    {
        return Err(GeomError::InvalidInput(format!(
            "chord tolerance must be positive and finite, got {chord_tolerance}"
        )));
    }
    // A line is exact between its endpoints: subdividing adds vertices that
    // carry no information.
    if let Curve3::Line(_) = curve {
        return Ok(vec![
            evaluate3(curve, domain.start)?,
            evaluate3(curve, domain.end)?,
        ]);
    }
    if let Curve3::Polyline(p) = curve {
        // A polyline's parameter is one unit per segment, so a caller passing
        // a normalized `(0, 1)` domain would silently collapse an n-vertex
        // path to its first edge. That is data loss disguised as success.
        let natural = polyline_domain(p.points.len(), p.closed);
        let requested = (domain.end - domain.start).abs();
        if natural.end > 1.0 && requested <= 1.0 {
            return Err(GeomError::InvalidInput(format!(
                "polyline domain {:?} spans {requested} of {} segments; a \
                 polyline parameter is one unit per segment, so this would \
                 discard {} vertices",
                domain,
                natural.end,
                p.points.len().saturating_sub(2)
            )));
        }
        return polyline_flatten(&p.points, p.closed, domain, |t| evaluate3(curve, t));
    }

    let eval = |t| evaluate3(curve, t);
    let mut out = vec![eval(domain.start)?];
    subdivide(
        &eval,
        domain.start,
        domain.end,
        chord_tolerance,
        max_depth.min(MAX_DEPTH_CEILING),
        MAX_POINTS,
        &mut out,
    )?;
    out.push(eval(domain.end)?);
    Ok(out)
}

/// Refuse a family that has no closed-form inversion.
///
/// Named separately from `unsupported_family` so the message can say WHY:
/// the family is understood, its inversion simply is not algebraic.
fn no_closed_form_inversion() -> GeomError {
    GeomError::InvalidInput(
        "curve family has no closed form inversion; a point trim on this basis \
         would require iteration, which trim resolution does not perform"
            .to_owned(),
    )
}

/// The point is not on the curve, so no parameter names it.
fn point_not_on_curve(distance: Scalar, tolerance: Scalar) -> GeomError {
    GeomError::InvalidInput(format!(
        "point is {distance} from the curve, outside the {tolerance} tolerance; \
         refusing rather than projecting it onto the nearest parameter"
    ))
}

/// Parameter of `point` on `line`, or a refusal when it is off the line.
///
/// The direction need not be unit length, so the projection is normalised by
/// its squared length: that is what makes the returned value a parameter of
/// THIS line rather than an arc length.
fn invert_line(
    origin: impl Into<[Scalar; 3]>,
    direction: [Scalar; 3],
    point: [Scalar; 3],
    tolerance: Scalar,
) -> GeomResult<Scalar> {
    let origin = origin.into();
    let dd = direction.iter().map(|c| c * c).sum::<Scalar>();
    if !dd.is_finite() || dd <= Scalar::EPSILON {
        return Err(GeomError::InvalidInput(
            "line direction is degenerate, so no parameter names a point".to_owned(),
        ));
    }
    let offset = [
        point[0] - origin[0],
        point[1] - origin[1],
        point[2] - origin[2],
    ];
    let t = offset
        .iter()
        .zip(direction.iter())
        .map(|(o, d)| o * d)
        .sum::<Scalar>()
        / dd;
    // Verify rather than assume: the projection always yields a parameter, but
    // only a point actually ON the line is named by it.
    let residual = [
        offset[0] - direction[0] * t,
        offset[1] - direction[1] * t,
        offset[2] - direction[2] * t,
    ];
    let distance = residual.iter().map(|c| c * c).sum::<Scalar>().sqrt();
    if distance > tolerance {
        return Err(point_not_on_curve(distance, tolerance));
    }
    Ok(t)
}

/// Parametric angle of `point` about a conic frame, verified against the curve.
///
/// For a circle the parametric angle is the polar angle; for an ellipse it is
/// not, so the local coordinates are divided by their semi-axes BEFORE the
/// arctangent. Taking the polar angle directly would be wrong off-axis.
fn invert_conic(
    local_x: Scalar,
    local_y: Scalar,
    semi_x: Scalar,
    semi_y: Scalar,
) -> GeomResult<Scalar> {
    if !(semi_x.is_finite() && semi_y.is_finite()) || semi_x <= 0.0 || semi_y <= 0.0 {
        return Err(GeomError::InvalidInput(
            "conic semi-axes must be finite and positive to invert a point".to_owned(),
        ));
    }
    let angle = (local_y / semi_y).atan2(local_x / semi_x);
    if !angle.is_finite() {
        return Err(GeomError::InvalidInput(
            "conic inversion produced a non-finite angle".to_owned(),
        ));
    }
    // Report on the same domain evaluation uses, so invert then evaluate is a
    // round trip rather than an off-by-one-turn surprise.
    Ok(angle.rem_euclid(std::f64::consts::TAU))
}

/// Parameter naming `point` on a 2D curve, or a refusal.
///
/// Exact for the families whose inversion is algebraic. Anything else is
/// refused by name: introducing iteration here would put a tolerance and a
/// convergence failure mode into every consumer of a point trim, and the
/// certified iterative path belongs to a caller that can carry its evidence.
///
/// # Errors
///
/// Refuses a point further than `tolerance` from the curve rather than
/// projecting it, and refuses families with no closed-form inversion.
pub fn invert2(curve: &Curve2, point: Point2, tolerance: Tolerance) -> GeomResult<Scalar> {
    finite2(point, "inversion point")?;
    let linear = tolerance.linear();
    match curve {
        Curve2::Line(l) => invert_line(
            [l.origin.x, l.origin.y, 0.0],
            [l.direction.x, l.direction.y, 0.0],
            [point.x, point.y, 0.0],
            linear,
        ),
        Curve2::Circle(c) => {
            let t = invert_conic_in_frame2(&c.frame, point, c.radius, c.radius)?;
            verify2(curve, t, point, linear)
        }
        Curve2::Ellipse(e) => {
            let t = invert_conic_in_frame2(&e.frame, point, e.semi_axis_x, e.semi_axis_y)?;
            verify2(curve, t, point, linear)
        }
        // A graph over its parameter: the parameter of a point IS its first
        // coordinate, then the height is checked.
        Curve2::Sinusoid(_) => verify2(curve, point.x, point, linear),
        _ => Err(no_closed_form_inversion()),
    }
}

/// Project a point into a 2D conic frame and invert it there.
fn invert_conic_in_frame2(
    frame: &axiolid_core::Frame2,
    point: Point2,
    semi_x: Scalar,
    semi_y: Scalar,
) -> GeomResult<Scalar> {
    let offset = point - frame.origin;
    invert_conic(offset.dot(frame.x), offset.dot(frame.y), semi_x, semi_y)
}

/// Confirm the recovered parameter actually reproduces the point.
///
/// The algebra above assumes an orthonormal frame. Imported frames are not
/// always orthonormal, and this crate deliberately keeps dirty frames
/// representable, so the claim is checked against the real evaluator instead
/// of trusted.
fn verify2(curve: &Curve2, t: Scalar, point: Point2, tolerance: Scalar) -> GeomResult<Scalar> {
    let found = evaluate2(curve, t)?;
    let distance = (found - point).length();
    if distance > tolerance {
        return Err(point_not_on_curve(distance, tolerance));
    }
    Ok(t)
}

/// Parameter naming `point` on a 3D curve, or a refusal.
///
/// See [`invert2`] for the exactness policy.
///
/// # Errors
///
/// Refuses an off-curve point and any family without a closed-form inversion.
pub fn invert3(curve: &Curve3, point: Point3, tolerance: Tolerance) -> GeomResult<Scalar> {
    finite3(point, "inversion point")?;
    let linear = tolerance.linear();
    match curve {
        Curve3::Line(l) => invert_line(
            [l.origin.x, l.origin.y, l.origin.z],
            [l.direction.x, l.direction.y, l.direction.z],
            [point.x, point.y, point.z],
            linear,
        ),
        Curve3::Circle(c) => {
            let t = invert_conic_in_frame3(&c.frame, point, c.radius, c.radius)?;
            verify3(curve, t, point, linear)
        }
        Curve3::Ellipse(e) => {
            let t = invert_conic_in_frame3(&e.frame, point, e.semi_axis_x, e.semi_axis_y)?;
            verify3(curve, t, point, linear)
        }
        _ => Err(no_closed_form_inversion()),
    }
}

/// Project a point into a 3D conic frame and invert it there.
fn invert_conic_in_frame3(
    frame: &axiolid_core::Frame3,
    point: Point3,
    semi_x: Scalar,
    semi_y: Scalar,
) -> GeomResult<Scalar> {
    let offset = point - frame.origin;
    invert_conic(offset.dot(frame.x), offset.dot(frame.y), semi_x, semi_y)
}

/// Confirm the recovered parameter reproduces the point. See [`verify2`].
fn verify3(curve: &Curve3, t: Scalar, point: Point3, tolerance: Scalar) -> GeomResult<Scalar> {
    let found = evaluate3(curve, t)?;
    let distance = (found - point).length();
    if distance > tolerance {
        return Err(point_not_on_curve(distance, tolerance));
    }
    Ok(t)
}
