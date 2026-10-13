//! Curve-path piece bounds given as a native parameter or a point (#304,
//! ADR 0082 amendment).
//!
//! A [`PathBound`] names where a piece starts or ends without stating its
//! distance in the curve's station measure. A parameter is read as it is;
//! the distance to it is [`StationCurve`]'s own measure at that parameter,
//! the one [`super::StationPiece::between_parameters`] measures by, so a
//! bounded piece resolves to the piece between those parameters bitwise.
//! This module resolves a POINT to the one parameter naming it.
//!
//! # Points
//!
//! A point is given in the curve's own coordinates, before any placement,
//! and lifted to `z = 0` on a 2D curve. It must lie on the curve within the
//! joint tolerance, [`super::JOINT_TOLERANCE`]` * max(1, |p|)`, `|p|` the
//! point's largest coordinate magnitude, the tolerance within which two
//! pieces meet and within which a graph compiler accepts a trim point:
//!
//! - A line, a circle or an ellipse is inverted in closed form (the
//!   parameter `axiolid-evaluate`'s `invert2`/`invert3` reports, a closed
//!   conic's angle in `[0, 2 pi)`), then checked against the curve. The
//!   angle runs round, so a closed conic's seam is one parameter, not two.
//! - A polyline is projected onto each segment; a B-spline is sampled
//!   (64 samples per control point, and every knot), each local minimum
//!   of the distance refined by Newton on `|C(t) - p|^2` within its
//!   neighbouring samples. Every parameter that comes within tolerance is
//!   kept.
//! - Two parameters name the same place when the curve between them stays
//!   within tolerance of the point (every vertex between them on a
//!   polyline; every sample between them and their midpoint on a
//!   B-spline): one place is one parameter, the nearest of them.
//!
//! No parameter within tolerance is refused naming [`PATH_POINT_OFF_CURVE`]
//! (the point is never projected onto the curve); several places naming
//! [`PATH_POINT_AMBIGUOUS`]: a curve passing through the point more than
//! once, and the seam of a closed polyline or B-spline, whose two domain
//! ends both name it. Every other family (an intrinsic curve, a chain, an
//! elevated or banked curve, an offset) is refused naming
//! [`POINT_BOUND_UNSUPPORTED`]; a caller states its parameter instead.

use axiolid_contracts::{BackendId, GeomError, GeomResult, Operation};
use axiolid_core::{Point2, Point3, Scalar, Tolerance, Vec3};
use axiolid_curve::{Curve2, Curve3, PathBound};
pub use axiolid_curve_evaluate_contract::{PATH_POINT_AMBIGUOUS, PATH_POINT_OFF_CURVE};

use super::composite::{StationCurve, JOINT_TOLERANCE};
use crate::curve::{
    derivative2, derivative3, domain2, domain3, evaluate2, evaluate3, invert2, invert3,
    second_derivative2, second_derivative3,
};

/// The input a point bound on a curve this evaluator does not invert is
/// refused for.
pub const POINT_BOUND_UNSUPPORTED: &str =
    "a curve path piece bounded by a point on a curve this evaluator does not invert (an \
     intrinsic curve, a chain, an elevated or banked curve, an offset, or a family other than a \
     line, a circle, an ellipse, a polyline or a B-spline): state its parameter instead";

/// Samples per control point a B-spline is scanned at for a point.
const SAMPLES_PER_CONTROL_POINT: usize = 64;

/// Newton steps refining one candidate parameter.
const NEWTON_STEPS: usize = 60;

fn refused(input: &'static str) -> GeomError {
    GeomError::UnsupportedInput {
        backend: BackendId::new("axiolid-evaluate"),
        operation: Operation::CurveEvaluation,
        input,
    }
}

/// How far from its curve a point bound may lie: the joint tolerance at
/// the point.
#[must_use]
pub fn point_bound_tolerance(point: Point3) -> Scalar {
    JOINT_TOLERANCE * point.abs().max_element().max(1.0)
}

/// The native parameter `bound` names on `curve`: a parameter as it is
/// (finite), a point resolved by the
/// [module documentation](self#points)'s rule.
///
/// # Errors
///
/// A parameter that is not finite, a point that is not finite, off the
/// curve ([`PATH_POINT_OFF_CURVE`]), named by several parameters
/// ([`PATH_POINT_AMBIGUOUS`]) or on a curve this evaluator does not invert
/// ([`POINT_BOUND_UNSUPPORTED`]), and a bound of a kind it does not know,
/// by name.
pub(super) fn bound_parameter(curve: StationCurve<'_>, bound: PathBound) -> GeomResult<Scalar> {
    match bound {
        PathBound::Parameter(t) if t.is_finite() => Ok(t),
        PathBound::Parameter(t) => Err(GeomError::Degenerate(format!(
            "station: a curve path piece is bounded at parameter {t}, which is not finite"
        ))),
        PathBound::Point(point) => point_parameter(curve, point),
        _ => Err(refused(
            "a curve path piece bound of a kind this evaluator does not know",
        )),
    }
}

/// The one parameter of `curve` naming `point`; see the
/// [module documentation](self#points).
fn point_parameter(curve: StationCurve<'_>, point: Point3) -> GeomResult<Scalar> {
    if !(point.x.is_finite() && point.y.is_finite() && point.z.is_finite()) {
        return Err(GeomError::InvalidInput(format!(
            "station: a curve path piece is bounded at point {point:?}, which is not finite"
        )));
    }
    let traced = Traced::of(curve)?;
    let found = traced.parameters(point, point_bound_tolerance(point))?;
    match found.as_slice() {
        [] => Err(refused(PATH_POINT_OFF_CURVE)),
        [t] => Ok(*t),
        _ => Err(refused(PATH_POINT_AMBIGUOUS)),
    }
}

/// A curve whose points this module inverts, 2D lifted to `z = 0`.
#[derive(Clone, Copy)]
enum Traced<'c> {
    Two(&'c Curve2),
    Three(&'c Curve3),
}

impl<'c> Traced<'c> {
    fn of(curve: StationCurve<'c>) -> GeomResult<Self> {
        match curve {
            StationCurve::Two(
                curve @ (Curve2::Line(_)
                | Curve2::Circle(_)
                | Curve2::Ellipse(_)
                | Curve2::Polyline(_)
                | Curve2::BSpline(_)),
            ) => Ok(Self::Two(curve)),
            StationCurve::Three(
                curve @ (Curve3::Line(_)
                | Curve3::Circle(_)
                | Curve3::Ellipse(_)
                | Curve3::Polyline(_)
                | Curve3::BSpline(_)),
            ) => Ok(Self::Three(curve)),
            _ => Err(refused(POINT_BOUND_UNSUPPORTED)),
        }
    }

    fn point(self, t: Scalar) -> GeomResult<Point3> {
        match self {
            Self::Two(curve) => {
                let p = evaluate2(curve, t)?;
                Ok(Point3::new(p.x, p.y, 0.0))
            }
            Self::Three(curve) => evaluate3(curve, t),
        }
    }

    fn derivatives(self, t: Scalar) -> GeomResult<(Vec3, Vec3)> {
        match self {
            Self::Two(curve) => {
                let (d, dd) = (derivative2(curve, t)?, second_derivative2(curve, t)?);
                Ok((Vec3::new(d.x, d.y, 0.0), Vec3::new(dd.x, dd.y, 0.0)))
            }
            Self::Three(curve) => Ok((derivative3(curve, t)?, second_derivative3(curve, t)?)),
        }
    }

    fn distance(self, t: Scalar, point: Point3) -> GeomResult<Scalar> {
        Ok((self.point(t)? - point).length())
    }

    /// Every place of the curve within `tolerance` of `point`, one
    /// parameter each, in increasing order.
    fn parameters(self, point: Point3, tolerance: Scalar) -> GeomResult<Vec<Scalar>> {
        match self {
            Self::Two(Curve2::Polyline(polyline)) => {
                let points: Vec<Point3> = polyline
                    .points
                    .iter()
                    .map(|p| Point3::new(p.x, p.y, 0.0))
                    .collect();
                Ok(polyline_parameters(
                    &points,
                    polyline.closed,
                    point,
                    tolerance,
                ))
            }
            Self::Three(Curve3::Polyline(polyline)) => Ok(polyline_parameters(
                &polyline.points,
                polyline.closed,
                point,
                tolerance,
            )),
            Self::Two(Curve2::BSpline(_)) | Self::Three(Curve3::BSpline(_)) => {
                self.spline_parameters(point, tolerance)
            }
            _ => {
                // Closed form: the parameter the inversion reports, at any
                // distance, then checked here so that a point off the curve
                // is refused by name.
                let anywhere = Tolerance::new(Scalar::MAX, 0.0).map_err(|_| {
                    GeomError::InvalidInput("station: no inversion tolerance".into())
                })?;
                let t = match self {
                    Self::Two(curve) => invert2(curve, Point2::new(point.x, point.y), anywhere)?,
                    Self::Three(curve) => invert3(curve, point, anywhere)?,
                };
                Ok(if self.distance(t, point)? <= tolerance {
                    vec![t]
                } else {
                    Vec::new()
                })
            }
        }
    }

    /// A B-spline's places within `tolerance` of `point`: every local
    /// minimum of the sampled distance, refined between its neighbours.
    fn spline_parameters(self, point: Point3, tolerance: Scalar) -> GeomResult<Vec<Scalar>> {
        let (domain, count, knots) = match self {
            Self::Two(curve @ Curve2::BSpline(spline)) => (
                domain2(curve),
                spline.control_points.len(),
                spline.knots.as_slice(),
            ),
            Self::Three(curve @ Curve3::BSpline(spline)) => (
                domain3(curve),
                spline.control_points.len(),
                spline.knots.as_slice(),
            ),
            _ => return Err(refused(POINT_BOUND_UNSUPPORTED)),
        };
        let (lo, hi) = (domain.start.min(domain.end), domain.start.max(domain.end));
        if !(lo.is_finite() && hi.is_finite() && hi > lo) {
            return Err(GeomError::Degenerate(
                "station: a B-spline with no parameter domain bounds a curve path piece".into(),
            ));
        }
        let samples = count.max(2) * SAMPLES_PER_CONTROL_POINT;
        let mut at: Vec<Scalar> = (0..=samples)
            .map(|i| lo + (hi - lo) * i as Scalar / samples as Scalar)
            .chain(knots.iter().copied().filter(|&k| k > lo && k < hi))
            .collect();
        at.sort_by(Scalar::total_cmp);
        at.dedup();
        let off = at
            .iter()
            .map(|&t| self.distance(t, point))
            .collect::<GeomResult<Vec<_>>>()?;
        let last = at.len() - 1;
        let mut found = Vec::new();
        for i in 0..=last {
            let falls = i == 0 || off[i] <= off[i - 1];
            let rises = i == last || off[i] <= off[i + 1];
            if !(falls && rises) {
                continue;
            }
            let (a, b) = (at[i.saturating_sub(1)], at[(i + 1).min(last)]);
            let refined = self.refine(at[i], a, b, point)?;
            let near = self.distance(refined, point)?;
            let (t, distance) = if near <= off[i] {
                (refined, near)
            } else {
                (at[i], off[i])
            };
            if distance <= tolerance {
                found.push((t, distance));
            }
        }
        found.sort_by(|x, y| x.0.total_cmp(&y.0));
        places(found, |ta, tb| {
            let between = at
                .iter()
                .zip(&off)
                .filter(|(t, _)| **t > ta && **t < tb)
                .all(|(_, d)| *d <= tolerance);
            Ok(between && self.distance(0.5 * (ta + tb), point)? <= tolerance)
        })
    }

    /// Newton on `(C(t) - p) . C'(t)` from `t`, kept within `[a, b]`.
    fn refine(self, mut t: Scalar, a: Scalar, b: Scalar, point: Point3) -> GeomResult<Scalar> {
        for _ in 0..NEWTON_STEPS {
            let (d, dd) = self.derivatives(t)?;
            let r = self.point(t)? - point;
            let (g, h) = (r.dot(d), d.dot(d) + r.dot(dd));
            if !(h.is_finite() && h > 0.0) {
                break;
            }
            let next = (t - g / h).clamp(a, b);
            let moved = (next - t).abs();
            t = next;
            if moved <= 4.0 * Scalar::EPSILON * (1.0 + t.abs()) {
                break;
            }
        }
        Ok(t)
    }
}

/// A polyline's places within `tolerance` of `point`: each segment's
/// nearest parameter (its index plus the fraction along it), one per
/// place.
fn polyline_parameters(
    points: &[Point3],
    closed: bool,
    point: Point3,
    tolerance: Scalar,
) -> Vec<Scalar> {
    let n = points.len();
    let segments = if closed { n } else { n.saturating_sub(1) };
    let mut found = Vec::new();
    for i in 0..segments {
        let (a, b) = (points[i], points[(i + 1) % n]);
        let d = b - a;
        let dd = d.dot(d);
        let u = if dd > 0.0 {
            ((point - a).dot(d) / dd).clamp(0.0, 1.0)
        } else {
            0.0
        };
        let distance = (a + d * u - point).length();
        if distance <= tolerance {
            found.push((i as Scalar + u, distance));
        }
    }
    let vertex = |k: usize| points[k % n.max(1)];
    places(found, |ta, tb| {
        // The vertices strictly between: the segments joining them stay
        // within tolerance exactly when every vertex does.
        let first = ta.floor() as usize + 1;
        let last = tb.ceil() as usize;
        Ok((first..last)
            .filter(|&k| (k as Scalar) > ta && (k as Scalar) < tb)
            .all(|k| (vertex(k) - point).length() <= tolerance))
    })
    .unwrap_or_default()
}

/// The candidates `(parameter, distance)`, in increasing parameter, merged
/// into places: consecutive candidates are one place when `same` says the
/// curve stays within tolerance between them. Each place is its nearest
/// candidate (the first on a tie).
fn places(
    found: Vec<(Scalar, Scalar)>,
    same: impl Fn(Scalar, Scalar) -> GeomResult<bool>,
) -> GeomResult<Vec<Scalar>> {
    let mut out: Vec<(Scalar, Scalar)> = Vec::new();
    let mut previous: Option<Scalar> = None;
    for (t, distance) in found {
        let joins = match previous {
            Some(before) => same(before, t)?,
            None => false,
        };
        previous = Some(t);
        match out.last_mut() {
            Some(place) if joins => {
                if distance < place.1 {
                    *place = (t, distance);
                }
            }
            _ => out.push((t, distance)),
        }
    }
    Ok(out.into_iter().map(|(t, _)| t).collect())
}
