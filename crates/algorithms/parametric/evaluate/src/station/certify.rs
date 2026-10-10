//! Certified lengths of station pieces, and so certified joint and seam
//! distances (#294, ADR 0082 amendment).
//!
//! A composite basis places its joints at the running sum of its pieces'
//! lengths ([`super::CompositeBasis`]). Most of those lengths are read from
//! stored data -- a line's, a circle's, a polyline's, an intrinsic curve's,
//! a chain's, an elevated or banked curve's plan distance, a constant
//! offset of a line or of a circle -- and are exact, rounding aside. The
//! others are quadratures, estimates to [`super::OFFSET_TOLERANCE`] or the
//! arc-length tolerance: an ellipse, a B-spline, and every other offset.
//! This module bounds those from both sides, so that a joint's distance is
//! a [`DistanceInterval`] that CONTAINS the true sum of the lengths before
//! it, measured between the points the kernel reads as each piece's ends.
//!
//! # Chords and turning
//!
//! Where the measured curve `P(t)` is evaluated in closed form at its
//! native parameter -- a line, a circle, an ellipse or a B-spline, and any
//! offset of one, 2D or 3D, measured by arc length -- each panel `[a, b]`,
//! `h = b - a`, is bounded by its chord `c = |P(b) - P(a)|` and a bound `M`
//! on `|P''|` there:
//!
//! - `c` is a lower bound of the arc length (a chord is never longer);
//! - `|P'(v)| >= m = c / h - M h / 2` everywhere on the panel, since
//!   `|c - h P'(v)| <= M h^2 / 2`;
//! - so the tangent turns by at most `beta = h M / m` over the panel, the
//!   chord direction lies within `beta` of every tangent, and the arc
//!   length is at most `c / cos(beta)`.
//!
//! The width `c (1 / cos(beta) - 1)` falls as `h^2`; panels are halved
//! until each is within its share of the tolerance. `M` comes from the
//! certified derivative suprema of [`crate::bound`] (`|B''|`, `|B'''|` of
//! the base) and, for an offset, the frame's: a unit vector `w / |w|` has
//! `|u'| <= |w'| / r` and `|u''| <= 2 |w''| / r + 6 |w'|^2 / r^2` where `r`
//! bounds `|w|` from below (`|B'|` for the tangent, `|B'_xy|` for the
//! reference-up lateral, `|V x B'|` for a 3D offset's direction), and
//! `(a e)'' = 2 a' e' + a e''` for a displacement `a` linear in the
//! parameter. Each point carries a rounding error of [`POINT_ROUNDING`]
//! relative to its largest coordinate, which widens every chord.
//!
//! # Midpoints of a speed in closed form
//!
//! Where the base is parameterised by its own arc length in the plan and
//! its curvature is a law -- an intrinsic curve or a chain (2D), the plan
//! of an elevated curve (measured by plan length) -- the offset's speed is
//! a closed form in the base's curvature `k` and, in a section frame that
//! leans with the grade, the grade `p = z'`: for a displacement `a` along
//! the left normal, `b` along the plan tangent and `c` vertical, all
//! linear in the base distance `v`,
//!
//! `f = sqrt((1 + b' - a k)^2 + (a' + b k)^2 + c'^2)`
//!
//! (`c'` only where the offset is measured in 3D, beside a 2D base), with
//! `b = longitudinal * cos g - vertical * sin g`, `tan g = p`, in a
//! section frame beside an elevated curve and `b = longitudinal`
//! otherwise. Each panel is the midpoint rule, `h f(mid)`, within `h^3 /
//! 24 sup |f''|`; `sup |f''|` is assembled from suprema of `k`, `k'`,
//! `k''` (the law and its derivatives) and of `z'` to `z''''` (the
//! profile's) by the Leibniz and chain rules, with `f >= inf (1 + b' -
//! a k) > 0` (a cusp is refused).
//!
//! # Refusals
//!
//! By a typed [`GeomError::UnsupportedInput`] naming
//! [`UNCERTIFIED_LENGTH`]: an offset of a banked curve, of a 3D intrinsic
//! curve, of a chain's parametric piece or of an elevated curve with an
//! intrinsic profile in a leaning frame; an offset by distances (a law that
//! varies) along an ellipse or a B-spline, whose law runs in the base's
//! numerical arc length; a 3D offset direction that is not vertical beside
//! a plan-parameterised base; a base placed in a frame that tilts `+Z`
//! there; a skewed frame; and a tolerance the rounding of the points
//! already exceeds. A panel budget of [`MAX_CERTIFIED_PANELS`] is
//! [`GeomError::BudgetExceeded`].

use std::borrow::Cow;

use axiolid_contracts::{BackendId, GeomError, GeomResult, Operation};
use axiolid_core::{Frame2, Point2, Point3, Scalar, Vec2, Vec3};
use axiolid_curve::{
    ChainPiece2, CurvatureLaw, Curve2, Curve3, ElevationLaw, Line2, Line3, OffsetFrame, OffsetLaw,
};
pub use axiolid_curve_evaluate_contract::DistanceInterval;

use super::composite::{StationCurve, StationPiece};
use super::offset::{OffsetTable, StationOffset};
use super::{
    reference_up_frame, station_length2, station_length3, StationSeam, KEEPS_UP_TOLERANCE,
};
use crate::arc_parameter::{
    parameter_at_arc_length2, parameter_at_arc_length3, ARC_LENGTH_TOLERANCE,
};
use crate::bound::{
    continuity_breaks2, continuity_breaks3, curve_derivative_bounds2, curve_derivative_bounds3,
    CurveDerivativeBounds,
};
use crate::curve::{
    derivative2, derivative3, domain2, domain3, evaluate2, evaluate3, second_derivative2,
    second_derivative3,
};

/// The input a length no certified bound here reaches is refused for.
pub const UNCERTIFIED_LENGTH: &str =
    "a certified length of a station piece no bound here reaches: an offset of a banked curve, a \
     3D intrinsic curve, a chain's parametric piece or an elevated curve with an intrinsic \
     profile in a leaning frame; an offset by distances along an ellipse or a B-spline; a 3D \
     offset direction that is not vertical beside a plan-measured or intrinsic base; a base \
     placed in a frame that tilts +Z; a skewed frame; or a tolerance below the points' rounding";

/// The rounding error of a point evaluated in closed form, relative to the
/// largest coordinate magnitude of the curve's data, the point and the
/// displacement (at least 1): about 45 units in the last place, several
/// times what a de Boor or conic evaluation and a frame accumulate.
pub const POINT_ROUNDING: Scalar = 1e-14;

/// The rounding error of a speed evaluated in closed form, relative.
const VALUE_ROUNDING: Scalar = 1e-14;

/// Most panels one certified length may use.
pub const MAX_CERTIFIED_PANELS: usize = 1 << 16;

/// Suprema of `|cos(atan p)|` and its first three derivatives in `p`.
const COS_OF_ATAN: [Scalar; 4] = [1.0, 0.39, 1.0, 1.94];

/// Suprema of `|sin(atan p)|` and its first three derivatives in `p`.
const SIN_OF_ATAN: [Scalar; 4] = [1.0, 1.0, 0.86, 3.0];

fn uncertified() -> GeomError {
    GeomError::UnsupportedInput {
        backend: BackendId::new("axiolid-evaluate"),
        operation: Operation::CurveEvaluation,
        input: UNCERTIFIED_LENGTH,
    }
}

fn budget() -> GeomError {
    GeomError::BudgetExceeded {
        resource: "certified length panels",
    }
}

fn slack(distance: Scalar) -> Scalar {
    ARC_LENGTH_TOLERANCE * distance.abs().max(1.0)
}

/// A tolerance a certified length can be asked for: finite and positive.
///
/// # Errors
///
/// Any other, by name.
pub(crate) fn checked_tolerance(tolerance: Scalar) -> GeomResult<Scalar> {
    if tolerance.is_finite() && tolerance > 0.0 {
        Ok(tolerance)
    } else {
        Err(GeomError::InvalidInput(format!(
            "station: a certified distance needs a finite, positive tolerance, not {tolerance}"
        )))
    }
}

// --- sums ---------------------------------------------------------------------

/// A running sum of lengths: exact while every term is, widened for the
/// rounding of its additions once one is not.
#[derive(Debug, Clone, Copy)]
pub(crate) struct Running {
    lower: Scalar,
    upper: Scalar,
    terms: usize,
    exact: bool,
}

impl Running {
    pub(crate) const ZERO: Self = Self {
        lower: 0.0,
        upper: 0.0,
        terms: 0,
        exact: true,
    };

    /// The sum with `term` added.
    pub(crate) fn plus(self, term: DistanceInterval) -> Self {
        Self {
            lower: self.lower + term.lower,
            upper: self.upper + term.upper,
            terms: self.terms + 1,
            exact: self.exact && term.is_point(),
        }
    }

    /// The certified interval of the sum.
    pub(crate) fn interval(self) -> DistanceInterval {
        if self.exact {
            return DistanceInterval::new(self.lower, self.upper);
        }
        let rounding = 2.0 * (self.terms as Scalar + 1.0) * Scalar::EPSILON;
        DistanceInterval::new(
            self.lower - rounding * self.lower.abs(),
            self.upper + rounding * self.upper.abs(),
        )
    }
}

// --- piece lengths ------------------------------------------------------------

/// The certified length of `piece` in its station measure: a point where
/// it is read from stored data, else an interval no wider than
/// `tolerance`; an offset read through `table` (built when `None`).
///
/// # Errors
///
/// The refusals of the [module documentation](self#refusals), and the
/// evaluators' own.
pub(crate) fn piece_length(
    piece: &StationPiece<'_>,
    table: Option<&OffsetTable>,
    tolerance: Scalar,
) -> GeomResult<DistanceInterval> {
    match piece.curve {
        StationCurve::Offset(offset) => {
            if offset.measure_is_exact() {
                return Ok(DistanceInterval::point(piece.length()));
            }
            let built;
            let table = match table {
                Some(table) => table,
                None => {
                    built = offset.table()?;
                    &built
                }
            };
            let a = offset.base_measure(piece.start, table)?;
            let b = offset.base_measure(piece.end, table)?;
            offset_length(&offset, a.min(b), a.max(b), tolerance)
        }
        curve if curve.measure_is_exact() => Ok(DistanceInterval::point(piece.length())),
        curve => curve_span(curve, piece.start, piece.end, tolerance),
    }
}

/// The certified arc length of `curve` (not an offset) between its station
/// measures `m0 <= m1`: exact where its measure is read from stored data,
/// else by chords between the native parameters the kernel reads them at.
///
/// # Errors
///
/// As [`piece_length`].
pub(crate) fn curve_span(
    curve: StationCurve<'_>,
    m0: Scalar,
    m1: Scalar,
    tolerance: Scalar,
) -> GeomResult<DistanceInterval> {
    if curve.measure_is_exact() {
        return Ok(DistanceInterval::point(m1 - m0));
    }
    let native = Native::of(curve, m0, m1)?;
    let (t0, t1) = (native.parameter(m0)?, native.parameter(m1)?);
    let breaks = native.breaks(1);
    chord_turning(
        t0.min(t1),
        t0.max(t1),
        &breaks,
        tolerance,
        &Measured {
            point: &|t| native.point(t),
            velocity: &|t| native.derivative(t),
            second: &|a, b| Ok(native.bounds(a, b)?.second),
            scale: native.extent(),
        },
    )
}

/// The certified length of `offset` between its base measures `va <= vb`
/// (into the base piece, in the direction of travel), in its own
/// convention.
fn offset_length(
    offset: &StationOffset<'_>,
    va: Scalar,
    vb: Scalar,
    tolerance: Scalar,
) -> GeomResult<DistanceInterval> {
    if vb <= va {
        return Ok(DistanceInterval::point(0.0));
    }
    if let Some(model) = PlanModel::of(offset)? {
        return model.length(va, vb, tolerance);
    }
    chord_offset_length(offset, va, vb, tolerance)
}

// --- closed-form points -------------------------------------------------------

/// A base curve evaluated in closed form at its native parameter, 2D lifted
/// to `z = 0`; a polyline's segment as a line in its measure.
enum Native<'c> {
    Two(Cow<'c, Curve2>),
    Three(Cow<'c, Curve3>),
}

impl<'c> Native<'c> {
    /// The closed-form curve under `curve` for the measures `[m0, m1]`.
    fn of(curve: StationCurve<'c>, m0: Scalar, m1: Scalar) -> GeomResult<Self> {
        let mid = 0.5 * (m0 + m1);
        match curve {
            StationCurve::Two(Curve2::Polyline(polyline)) => {
                let (origin, unit) = segment(&polyline.points, polyline.closed, mid, |p| {
                    Vec3::new(p.x, p.y, 0.0)
                })?;
                Ok(Self::Two(Cow::Owned(Curve2::Line(Line2 {
                    origin: Point2::new(origin.x, origin.y),
                    direction: Vec2::new(unit.x, unit.y),
                }))))
            }
            StationCurve::Three(Curve3::Polyline(polyline)) => {
                let (origin, unit) = segment(&polyline.points, polyline.closed, mid, |p| p)?;
                Ok(Self::Three(Cow::Owned(Curve3::Line(Line3 {
                    origin,
                    direction: unit,
                }))))
            }
            StationCurve::Two(
                curve @ (Curve2::Line(_)
                | Curve2::Circle(_)
                | Curve2::Ellipse(_)
                | Curve2::BSpline(_)),
            ) => Ok(Self::Two(Cow::Borrowed(curve))),
            StationCurve::Three(
                curve @ (Curve3::Line(_)
                | Curve3::Circle(_)
                | Curve3::Ellipse(_)
                | Curve3::BSpline(_)),
            ) => Ok(Self::Three(Cow::Borrowed(curve))),
            _ => Err(uncertified()),
        }
    }

    fn point(&self, t: Scalar) -> GeomResult<Point3> {
        match self {
            Self::Two(curve) => {
                let p = evaluate2(curve, t)?;
                Ok(Point3::new(p.x, p.y, 0.0))
            }
            Self::Three(curve) => evaluate3(curve, t),
        }
    }

    fn derivative(&self, t: Scalar) -> GeomResult<Vec3> {
        match self {
            Self::Two(curve) => {
                let d = derivative2(curve, t)?;
                Ok(Vec3::new(d.x, d.y, 0.0))
            }
            Self::Three(curve) => derivative3(curve, t),
        }
    }

    fn second_derivative(&self, t: Scalar) -> GeomResult<Vec3> {
        match self {
            Self::Two(curve) => {
                let d = second_derivative2(curve, t)?;
                Ok(Vec3::new(d.x, d.y, 0.0))
            }
            Self::Three(curve) => second_derivative3(curve, t),
        }
    }

    /// The largest coordinate magnitude of the curve's data: a B-spline's
    /// control points, a conic's centre plus its axes, a line's origin.
    /// What a point's rounding is relative to.
    fn extent(&self) -> Scalar {
        let flat = |v: Vec2| v.abs().max_element();
        let out = match self {
            Self::Two(curve) => match curve.as_ref() {
                Curve2::Line(line) => flat(line.origin),
                Curve2::Circle(c) => {
                    flat(c.frame.origin) + c.radius * (flat(c.frame.x) + flat(c.frame.y))
                }
                Curve2::Ellipse(e) => {
                    flat(e.frame.origin)
                        + e.semi_axis_x * flat(e.frame.x)
                        + e.semi_axis_y * flat(e.frame.y)
                }
                Curve2::BSpline(spline) => spline
                    .control_points
                    .iter()
                    .map(|p| flat(*p))
                    .fold(0.0, Scalar::max),
                _ => Scalar::INFINITY,
            },
            Self::Three(curve) => match curve.as_ref() {
                Curve3::Line(line) => line.origin.abs().max_element(),
                Curve3::Circle(c) => {
                    c.frame.origin.abs().max_element()
                        + c.radius * (c.frame.x.abs().max_element() + c.frame.y.abs().max_element())
                }
                Curve3::Ellipse(e) => {
                    e.frame.origin.abs().max_element()
                        + e.semi_axis_x * e.frame.x.abs().max_element()
                        + e.semi_axis_y * e.frame.y.abs().max_element()
                }
                Curve3::BSpline(spline) => spline
                    .control_points
                    .iter()
                    .map(|p| p.abs().max_element())
                    .fold(0.0, Scalar::max),
                _ => Scalar::INFINITY,
            },
        };
        out.max(1.0)
    }

    /// Certified suprema of `|B'|`, `|B''|`, `|B'''|` over `[a, b]`.
    fn bounds(&self, a: Scalar, b: Scalar) -> GeomResult<CurveDerivativeBounds> {
        match self {
            Self::Two(curve) => curve_derivative_bounds2(curve, a, b),
            Self::Three(curve) => curve_derivative_bounds3(curve, a, b),
        }
        .filter(|d| d.first.is_finite() && d.second.is_finite() && d.third.is_finite())
        .ok_or_else(uncertified)
    }

    /// Where the curve may fail to be `C^k`.
    fn breaks(&self, k: usize) -> Vec<Scalar> {
        match self {
            Self::Two(curve) => continuity_breaks2(curve, k),
            Self::Three(curve) => continuity_breaks3(curve, k),
        }
    }

    /// `dm/dt`, the station measure's rate in the native parameter, where
    /// it is a constant: a line's `|d|`, a circle's radius.
    fn rate(&self) -> Option<Scalar> {
        match self {
            Self::Two(curve) => match curve.as_ref() {
                Curve2::Line(line) => Some(line.direction.length()),
                Curve2::Circle(circle) => Some(circle.radius),
                _ => None,
            },
            Self::Three(curve) => match curve.as_ref() {
                Curve3::Line(line) => Some(line.direction.length()),
                Curve3::Circle(circle) => Some(circle.radius),
                _ => None,
            },
        }
    }

    /// The native parameter the kernel reads the station measure `m` at:
    /// `m / rate` on a line or a circle; on an ellipse or a B-spline, the
    /// curve's start, end and seams (corner knots) exactly where `m` is
    /// within the arc-length tolerance of them, the arc-length inverse
    /// elsewhere (an ellipse's measure may run past one turn, into the
    /// next).
    fn parameter(&self, m: Scalar) -> GeomResult<Scalar> {
        if let Some(rate) = self.rate() {
            return Ok(m / rate);
        }
        let (length, lo, hi, seams) = match self {
            Self::Two(curve) => {
                let domain = domain2(curve);
                (
                    station_length2(curve)?,
                    domain.start.min(domain.end),
                    domain.start.max(domain.end),
                    super::station_seams2(curve)?,
                )
            }
            Self::Three(curve) => {
                let domain = domain3(curve);
                (
                    station_length3(curve)?,
                    domain.start.min(domain.end),
                    domain.start.max(domain.end),
                    super::station_seams3(curve)?,
                )
            }
        };
        let length = length.ok_or_else(uncertified)?;
        let closed = match self {
            Self::Two(curve) => matches!(curve.as_ref(), Curve2::Ellipse(_)),
            Self::Three(curve) => matches!(curve.as_ref(), Curve3::Ellipse(_)),
        };
        let (turn, m) = if closed && m > length + slack(length) {
            (hi - lo, m - length)
        } else {
            (0.0, m)
        };
        let t = if m.abs() <= slack(0.0) {
            lo
        } else if (m - length).abs() <= slack(length) {
            hi
        } else if let Some(seam) = seams
            .iter()
            .find(|seam| (seam.distance - m).abs() <= slack(m))
        {
            seam.parameter
        } else {
            match self {
                Self::Two(curve) => parameter_at_arc_length2(curve, lo, m)?,
                Self::Three(curve) => parameter_at_arc_length3(curve, lo, m)?,
            }
        };
        Ok(turn + t)
    }
}

/// A polyline's segment holding the measure `mid`, as the origin and unit
/// direction of the line whose own measure is the polyline's there.
fn segment<P: Copy>(
    points: &[P],
    closed: bool,
    mid: Scalar,
    lift: impl Fn(P) -> Vec3,
) -> GeomResult<(Point3, Vec3)> {
    let count = points.len();
    let spans = if closed {
        count
    } else {
        count.saturating_sub(1)
    };
    let mut run = 0.0;
    let mut found = None;
    for i in 0..spans {
        let (a, b) = (lift(points[i]), lift(points[(i + 1) % count]));
        let length = (b - a).length();
        if length > 0.0 {
            found = Some((a, (b - a) / length, run));
            if mid <= run + length {
                break;
            }
        }
        run += length;
    }
    let (start, unit, at) = found.ok_or_else(uncertified)?;
    Ok((start - unit * at, unit))
}

/// `|u'|` and `|u''|` of a unit vector `u = w / |w|`, from a lower bound
/// `r` of `|w|` and suprema of `|w'|`, `|w''|`: `|w'| / r` and `|w''| / r +
/// 3 |w'|^2 / r^2`; infinite where `r` is not positive.
fn unit_bounds(r: Scalar, w1: Scalar, w2: Scalar) -> (Scalar, Scalar) {
    if r.is_nan() || r <= 0.0 {
        return (Scalar::INFINITY, Scalar::INFINITY);
    }
    (w1 / r, w2 / r + 3.0 * w1 * w1 / (r * r))
}

/// A unit vector `w / |w|` and its derivative, from `w` and `w'`.
fn unit_rate(w: Vec3, dw: Vec3) -> GeomResult<(Vec3, Vec3)> {
    let r = w.length();
    if r.is_nan() || r <= 0.0 {
        return Err(uncertified());
    }
    let u = w / r;
    Ok((u, (dw - u * u.dot(dw)) / r))
}

/// The certified length of an offset of a closed-form base between its
/// base measures `va < vb`, by chords and turning.
fn chord_offset_length(
    offset: &StationOffset<'_>,
    va: Scalar,
    vb: Scalar,
    tolerance: Scalar,
) -> GeomResult<DistanceInterval> {
    let base = offset.base();
    let measure = |v: Scalar| {
        if base.reversed {
            base.end - v
        } else {
            base.start + v
        }
    };
    let (ma, mb) = (measure(va), measure(vb));
    let native = Native::of(base.curve, ma.min(mb), ma.max(mb))?;
    let law = offset.law();
    let rate = native.rate();
    if !law.is_constant() && rate.is_none() {
        return Err(uncertified());
    }
    let travel = if base.reversed { -1.0 } else { 1.0 };
    let placement = base.placement;
    let length = offset.base_length();
    // The base distance into the piece at a native parameter, and its rate;
    // read only by a law that varies, which needs a line or a circle.
    let into = |t: Scalar| match rate {
        Some(rate) => {
            let m = t * rate;
            if base.reversed {
                base.end - m
            } else {
                m - base.start
            }
        }
        None => 0.0,
    };
    let into_rate = travel * rate.unwrap_or(0.0);
    let move_point = |p: Point3| placement.map_or(p, |rigid| rigid.transform_point3(p));
    let move_vector = |v: Vec3| placement.map_or(v, |rigid| rigid.transform_vector3(v));
    let point = |t: Scalar| -> GeomResult<Point3> {
        let p = move_point(native.point(t)?);
        let d = move_vector(native.derivative(t)?);
        let section = reference_up_frame(p, travel * d)?;
        Ok(section.point + offset.displacement(&section, into(t))?)
    };
    // `P' = B' + D'`, the displacement's frame differentiated through the
    // unit vectors it is built of.
    let velocity = |t: Scalar| -> GeomResult<Vec3> {
        let d = move_vector(native.derivative(t)?);
        let (w, dw) = (
            travel * d,
            travel * move_vector(native.second_derivative(t)?),
        );
        let (tangent, dtangent) = unit_rate(w, dw)?;
        let (lateral, dlateral) = unit_rate(Vec3::Z.cross(w), Vec3::Z.cross(dw))?;
        let up = (
            tangent.cross(lateral),
            dtangent.cross(lateral) + tangent.cross(dlateral),
        );
        let displacement = match law {
            OffsetLaw::Planar { distance } => distance * dlateral,
            OffsetLaw::Directed {
                distance,
                reference_direction,
            } => {
                let v = reference_direction.normalize();
                distance * unit_rate(v.cross(w), v.cross(dw))?.1
            }
            OffsetLaw::Linear { start, end, frame } => {
                let axes = match frame {
                    OffsetFrame::Section => [(lateral, dlateral), up, (tangent, dtangent)],
                    OffsetFrame::Plan => {
                        let flat = |v: Vec3| Vec3::new(v.x, v.y, 0.0);
                        let (along, dalong) = unit_rate(flat(w), flat(dw))?;
                        [
                            (Vec3::Z.cross(along), Vec3::Z.cross(dalong)),
                            (Vec3::Z, Vec3::ZERO),
                            (along, dalong),
                        ]
                    }
                    _ => return Err(uncertified()),
                };
                let u = (into(t) / length).clamp(0.0, 1.0);
                let pairs = [
                    (start.lateral, end.lateral),
                    (start.vertical, end.vertical),
                    (start.longitudinal, end.longitudinal),
                ];
                pairs
                    .iter()
                    .zip(axes)
                    .map(|(&(a, b), (axis, daxis))| {
                        (b - a) / length * into_rate * axis + (a + (b - a) * u) * daxis
                    })
                    .fold(Vec3::ZERO, |acc, v| acc + v)
            }
            _ => return Err(uncertified()),
        };
        Ok(d + displacement)
    };
    let (dv, sizes) = law_rates(law, rate.unwrap_or(0.0) / length)?;
    let second = |a: Scalar, b: Scalar| -> GeomResult<Scalar> {
        let bounds = native.bounds(a, b)?;
        let h = b - a;
        let d = move_vector(native.derivative(0.5 * (a + b))?);
        let reach = 0.5 * h * bounds.second;
        let flat = Vec3::new(d.x, d.y, 0.0).length();
        let tangent = unit_bounds(d.length() - reach, bounds.second, bounds.third);
        let lateral = unit_bounds(flat - reach, bounds.second, bounds.third);
        let up = (
            tangent.0 + lateral.0,
            tangent.1 + 2.0 * tangent.0 * lateral.0 + lateral.1,
        );
        let displacement = match law {
            OffsetLaw::Planar { distance } => distance.abs() * lateral.1,
            OffsetLaw::Directed {
                distance,
                reference_direction,
            } => {
                let v = reference_direction.normalize();
                let r = v.cross(d).length() - reach;
                distance.abs() * unit_bounds(r, bounds.second, bounds.third).1
            }
            OffsetLaw::Linear { frame, .. } => {
                let axes = match frame {
                    OffsetFrame::Section => [lateral, up, tangent],
                    OffsetFrame::Plan => [lateral, (0.0, 0.0), lateral],
                    _ => return Err(uncertified()),
                };
                (0..3)
                    .map(|k| 2.0 * dv[k] * axes[k].0 + sizes[k] * axes[k].1)
                    .sum()
            }
            _ => return Err(uncertified()),
        };
        Ok(bounds.second + displacement)
    };
    let mut breaks = native.breaks(2);
    breaks.extend(native.breaks(1));
    let reach = sizes.iter().sum::<Scalar>()
        + match law {
            OffsetLaw::Planar { distance } | OffsetLaw::Directed { distance, .. } => distance.abs(),
            _ => 0.0,
        };
    let shift = placement.map_or(0.0, |rigid| rigid.translation.abs().max_element());
    let (ta, tb) = (native.parameter(ma)?, native.parameter(mb)?);
    chord_turning(
        ta.min(tb),
        ta.max(tb),
        &breaks,
        tolerance,
        &Measured {
            point: &point,
            velocity: &velocity,
            second: &second,
            scale: native.extent() + reach + shift,
        },
    )
}

/// A law's components `(lateral, vertical, longitudinal)`: their rates in
/// the native parameter (`rate` the base distance per parameter, divided
/// by the base length) and their largest sizes.
fn law_rates(law: OffsetLaw, rate: Scalar) -> GeomResult<([Scalar; 3], [Scalar; 3])> {
    Ok(match law {
        OffsetLaw::Linear { start, end, .. } => {
            let pairs = [
                (start.lateral, end.lateral),
                (start.vertical, end.vertical),
                (start.longitudinal, end.longitudinal),
            ];
            (
                pairs.map(|(a, b)| (b - a).abs() * rate.abs()),
                pairs.map(|(a, b)| a.abs().max(b.abs())),
            )
        }
        OffsetLaw::Planar { .. } | OffsetLaw::Directed { .. } => ([0.0; 3], [0.0; 3]),
        _ => return Err(uncertified()),
    })
}

/// A curve measured by chords: its point, its velocity and a bound on its
/// second derivative over a parameter interval, all in the native
/// parameter, and the coordinate magnitude a point's rounding is relative
/// to.
struct Measured<'a> {
    point: &'a dyn Fn(Scalar) -> GeomResult<Point3>,
    velocity: &'a dyn Fn(Scalar) -> GeomResult<Vec3>,
    second: &'a dyn Fn(Scalar, Scalar) -> GeomResult<Scalar>,
    scale: Scalar,
}

/// A sum of many positive terms, compensated (Neumaier): its error is a
/// few units in the last place of the sum, whatever the number of terms.
#[derive(Debug, Clone, Copy, Default)]
struct Compensated {
    sum: Scalar,
    carry: Scalar,
}

impl Compensated {
    fn add(&mut self, value: Scalar) {
        let total = self.sum + value;
        self.carry += if self.sum.abs() >= value.abs() {
            (self.sum - total) + value
        } else {
            (value - total) + self.sum
        };
        self.sum = total;
    }

    fn value(self) -> Scalar {
        self.sum + self.carry
    }
}

/// The certified length of `P` over the native parameters `[t0, t1]`, by
/// chords and turning (see the [module documentation](self)), split at
/// `breaks` (where `P'` may not be continuous) and halved until each panel
/// is within its share of `tolerance`.
fn chord_turning(
    t0: Scalar,
    t1: Scalar,
    breaks: &[Scalar],
    tolerance: Scalar,
    curve: &Measured<'_>,
) -> GeomResult<DistanceInterval> {
    let span = t1 - t0;
    if span.is_nan() || span <= 0.0 {
        return Ok(DistanceInterval::point(0.0));
    }
    let mut cuts = vec![t0];
    let mut inside: Vec<Scalar> = breaks
        .iter()
        .copied()
        .filter(|&t| t > t0 && t < t1)
        .collect();
    inside.sort_by(Scalar::total_cmp);
    inside.dedup();
    cuts.extend(inside);
    cuts.push(t1);
    let (mut chords, mut excesses) = (Compensated::default(), Compensated::default());
    let mut rounding = Rounding::new(curve.scale);
    let mut panels = 0_usize;
    for pair in cuts.windows(2) {
        let (a0, b0) = (pair[0], pair[1]);
        // Depth first, left half first: panels are accepted in order along
        // the curve, each starting where the last ended.
        let mut stack = vec![(a0, b0, (curve.point)(a0)?, (curve.point)(b0)?)];
        while let Some((a, b, pa, pb)) = stack.pop() {
            panels += 1;
            if panels > MAX_CERTIFIED_PANELS {
                return Err(budget());
            }
            let h = b - a;
            let mid = 0.5 * (a + b);
            let chord = (pb - pa).length();
            // Every tangent on the panel lies within `turn` of the
            // midpoint's, so the chord is at least `cos(turn)` of the arc.
            let bound = (curve.second)(a, b)?;
            let speed = (curve.velocity)(mid)?.length() - 0.5 * h * bound;
            let turn = if speed > 0.0 {
                0.5 * h * bound / speed + 1e-12
            } else {
                Scalar::INFINITY
            };
            // `c / cos(turn) - c`, without cancellation.
            let excess = if turn < 1.0 {
                let half = (0.5 * turn).sin();
                chord * (1.0 + 4.0 * Scalar::EPSILON) * 2.0 * half * half / turn.cos()
            } else {
                Scalar::INFINITY
            };
            if excess <= 0.5 * tolerance * h / span {
                chords.add(chord);
                excesses.add(excess);
                rounding.panel(pa, pb, chord);
                continue;
            }
            if !(mid > a && mid < b) {
                return Err(uncertified());
            }
            let pm = (curve.point)(mid)?;
            stack.push((mid, b, pm, pb));
            stack.push((a, mid, pa, pm));
        }
    }
    let error = rounding.bound();
    let sum = chords.value();
    finish(sum - error, sum + excesses.value() + error, tolerance)
}

/// The rounding of a chain of chords between computed points (#294).
///
/// Each point is off its curve by at most `e = POINT_ROUNDING * scale`.
/// Writing a chord `|D + d|`, `D` the true one and `d` the difference of
/// its ends' errors (`|d| <= 2e`): `|D + d| - |D| = d.u + q`, `u = D /
/// |D|`, with `0 <= q <= |d|^2 / (2 |D|)`. The parts `d.u` telescope along
/// the chain: their sum is at most `e (2 + sum |u_i - u_(i-1)|)`, the
/// directions' total variation, which the chords themselves bound. So the
/// rounding of the whole chain stays a few `e`, however many chords it has;
/// each chord's own norm adds `4 eps |c|`.
#[derive(Debug, Clone, Copy)]
struct Rounding {
    scale: Scalar,
    previous: Option<Vec3>,
    variation: Scalar,
    quadratic: Scalar,
    norms: Scalar,
}

impl Rounding {
    fn new(scale: Scalar) -> Self {
        Self {
            scale: scale.max(1.0),
            previous: None,
            variation: 0.0,
            quadratic: 0.0,
            norms: 0.0,
        }
    }

    /// An accepted panel from `pa` to `pb`, its computed chord `chord`.
    fn panel(&mut self, pa: Point3, pb: Point3, chord: Scalar) {
        self.scale = self
            .scale
            .max(pa.abs().max_element())
            .max(pb.abs().max_element());
        let e = POINT_ROUNDING * self.scale;
        if chord > 4.0 * e {
            let direction = (pb - pa) / chord;
            if let Some(previous) = self.previous {
                self.variation += (direction - previous).length();
            }
            self.previous = Some(direction);
            // The chord's second-order part, and the error of the computed
            // direction against the true one (at most `4 e / |D|`, twice).
            self.quadratic += 10.0 * e * e / (chord - 2.0 * e);
        } else {
            // Too short to have a direction: count its whole error.
            self.variation += 4.0;
            self.previous = None;
        }
        self.norms += 4.0 * Scalar::EPSILON * chord;
    }

    /// A bound on the difference of the computed and the true chord sums.
    fn bound(&self) -> Scalar {
        let e = POINT_ROUNDING * self.scale;
        e * (2.0 + self.variation) + self.quadratic + self.norms
    }
}

/// A compensated sum of panel bounds as an interval, widened for its own
/// rounding; refused when wider than `tolerance`.
fn finish(lower: Scalar, upper: Scalar, tolerance: Scalar) -> GeomResult<DistanceInterval> {
    let rounding = 4.0 * Scalar::EPSILON;
    let out = DistanceInterval::new(
        (lower * (1.0 - rounding)).max(0.0),
        upper * (1.0 + rounding),
    );
    if !(out.lower.is_finite() && out.upper.is_finite()) || out.width() > tolerance {
        return Err(uncertified());
    }
    Ok(out)
}

// --- speeds in closed form ------------------------------------------------------

/// Suprema of a function and its derivatives over a panel, `[|f|, |f'|,
/// ...]`.
type Jet<const N: usize> = [Scalar; N];

/// Leibniz: suprema of `(f g)^(k)`.
fn product<const N: usize>(f: Jet<N>, g: Jet<N>) -> Jet<N> {
    let binomial = |n: usize, k: usize| -> Scalar {
        (0..k).fold(1.0, |acc, i| acc * (n - i) as Scalar / (i + 1) as Scalar)
    };
    core::array::from_fn(|k| (0..=k).map(|j| binomial(k, j) * f[j] * g[k - j]).sum())
}

fn sum<const N: usize>(f: Jet<N>, g: Jet<N>) -> Jet<N> {
    core::array::from_fn(|k| f[k] + g[k])
}

/// Faa di Bruno to the third order: suprema of `(outer o p)^(k)` from the
/// suprema `outer` of the outer function's derivatives and `p` of the
/// inner's.
fn composed(outer: Jet<4>, p: Jet<4>) -> Jet<4> {
    [
        outer[0],
        outer[1] * p[1],
        outer[2] * p[1] * p[1] + outer[1] * p[2],
        outer[3] * p[1].powi(3) + 3.0 * outer[2] * p[1] * p[2] + outer[1] * p[3],
    ]
}

/// A displacement component linear in the base distance over `[0,
/// length]`: its suprema `[max(|a|, |b|), |b - a| / length, 0, 0]`.
fn linear(pair: (Scalar, Scalar), length: Scalar) -> Jet<4> {
    [
        pair.0.abs().max(pair.1.abs()),
        (pair.1 - pair.0).abs() / length,
        0.0,
        0.0,
    ]
}

/// Where a plan-parameterised base's curvature comes from.
#[derive(Debug, Clone, Copy)]
enum Curvature<'a> {
    /// The same everywhere: a line (zero) or a circle (`+-1 / r`).
    Constant(Scalar),
    /// `sign * law(s - shift)`: an intrinsic curve or a chain's intrinsic
    /// piece, `sign` the start frame's handedness.
    Law {
        law: &'a CurvatureLaw,
        shift: Scalar,
        sign: Scalar,
    },
}

impl Curvature<'_> {
    fn value(&self, s: Scalar) -> GeomResult<Scalar> {
        match self {
            Self::Constant(k) => Ok(*k),
            Self::Law { law, shift, sign } => crate::frenet::law_value(law, s - shift)
                .map(|k| sign * k)
                .ok_or_else(uncertified),
        }
    }

    /// Suprema of `|k|`, `|k'|`, `|k''|` over the base measures `[lo, hi]`.
    fn jet(&self, lo: Scalar, hi: Scalar) -> GeomResult<Jet<3>> {
        match self {
            Self::Constant(k) => Ok([k.abs(), 0.0, 0.0]),
            Self::Law { law, shift, .. } => {
                let (lo, hi) = (lo - shift, hi - shift);
                let first = law.derivative();
                let second = first.derivative();
                let bound = |law: &CurvatureLaw| {
                    crate::chain::curvature_bound(law, lo, hi)
                        .filter(|b| b.is_finite())
                        .ok_or_else(uncertified)
                };
                Ok([bound(law)?, bound(&first)?, bound(&second)?])
            }
        }
    }
}

/// The handedness of an orthonormal 2D frame, `+1` or `-1`.
///
/// # Errors
///
/// A frame that is not orthonormal (its measure would not be its curve's
/// arc length), as uncertified.
fn handedness(frame: &Frame2) -> GeomResult<Scalar> {
    let skew = (frame.x.length() - 1.0)
        .abs()
        .max((frame.y.length() - 1.0).abs())
        .max(frame.x.dot(frame.y).abs());
    if skew.is_nan() || skew > 1e-12 {
        return Err(uncertified());
    }
    Ok(frame.x.perp_dot(frame.y).signum())
}

/// The curvature of a plan, parameterised by its arc length, over the
/// measures `[lo, hi]` (inside one piece of it).
fn plan_curvature(plan: &Curve2, lo: Scalar, hi: Scalar) -> GeomResult<Curvature<'_>> {
    match plan {
        Curve2::Line(_) | Curve2::Polyline(_) => Ok(Curvature::Constant(0.0)),
        Curve2::Circle(circle) => {
            let sign = handedness(&circle.frame)?;
            if !(circle.radius.is_finite() && circle.radius > 0.0) {
                return Err(uncertified());
            }
            Ok(Curvature::Constant(sign / circle.radius))
        }
        Curve2::Intrinsic(intrinsic) => Ok(Curvature::Law {
            law: &intrinsic.curvature,
            shift: 0.0,
            sign: handedness(&intrinsic.start)?,
        }),
        Curve2::Chain(chain) => {
            let sign = handedness(&chain.start)?;
            let mid = 0.5 * (lo + hi);
            match chain.piece_at(mid) {
                Some((_, ChainPiece2::Intrinsic { curvature, .. }, local)) => Ok(Curvature::Law {
                    law: curvature,
                    shift: mid - local,
                    sign,
                }),
                _ => Err(uncertified()),
            }
        }
        _ => Err(uncertified()),
    }
}

/// An offset beside a base parameterised by its plan arc length, its speed
/// in closed form (see the [module documentation](self)).
struct PlanModel<'a> {
    curvature: Curvature<'a>,
    /// The profile, where the frame leans with its grade and the law has a
    /// vertical or longitudinal part.
    grade: Option<&'a ElevationLaw>,
    lateral: (Scalar, Scalar),
    vertical: (Scalar, Scalar),
    longitudinal: (Scalar, Scalar),
    /// Whether a vertical rate lengthens the offset: beside a 2D base,
    /// measured in 3D.
    vertical_rate: bool,
    start: Scalar,
    end: Scalar,
    reversed: bool,
}

impl<'a> PlanModel<'a> {
    /// The model of `offset`, `None` where its base is not parameterised by
    /// its plan arc length with a curvature law (an intrinsic curve or a
    /// chain, or an elevated curve).
    fn of(offset: &StationOffset<'a>) -> GeomResult<Option<Self>> {
        let base = offset.base();
        let (lo, hi) = (base.start, base.end);
        let (curvature, elevation, vertical_rate) = match base.curve {
            StationCurve::Two(plan @ (Curve2::Intrinsic(_) | Curve2::Chain(_))) => {
                (plan_curvature(plan, lo, hi)?, None, true)
            }
            StationCurve::Three(Curve3::Elevated(elevated)) => (
                plan_curvature(&elevated.plan, lo, hi)?,
                Some(&elevated.elevation),
                false,
            ),
            StationCurve::Three(Curve3::Banked(_) | Curve3::Intrinsic(_)) => {
                return Err(uncertified())
            }
            _ => return Ok(None),
        };
        if let Some(rigid) = base.placement {
            if (rigid.transform_vector3(Vec3::Z) - Vec3::Z).length() > KEEPS_UP_TOLERANCE {
                return Err(uncertified());
            }
        }
        let zero = (0.0, 0.0);
        let (lateral, vertical, longitudinal, leans) = match offset.law() {
            OffsetLaw::Planar { distance } => ((distance, distance), zero, zero, false),
            OffsetLaw::Directed {
                distance,
                reference_direction,
            } => {
                let v = reference_direction.normalize();
                if v.x.hypot(v.y) > 1e-12 {
                    return Err(uncertified());
                }
                let a = distance * v.z.signum();
                ((a, a), zero, zero, false)
            }
            OffsetLaw::Linear { start, end, frame } => {
                let leans = match frame {
                    OffsetFrame::Section => true,
                    OffsetFrame::Plan => false,
                    _ => return Err(uncertified()),
                };
                (
                    (start.lateral, end.lateral),
                    (start.vertical, end.vertical),
                    (start.longitudinal, end.longitudinal),
                    leans,
                )
            }
            _ => return Err(uncertified()),
        };
        let grade = elevation.filter(|_| leans && (vertical != zero || longitudinal != zero));
        Ok(Some(Self {
            curvature,
            grade,
            lateral,
            vertical,
            longitudinal,
            vertical_rate,
            start: base.start,
            end: base.end,
            reversed: base.reversed,
        }))
    }

    fn length_of_base(&self) -> Scalar {
        self.end - self.start
    }

    /// The base's own measure at `v` into the piece.
    fn measure(&self, v: Scalar) -> Scalar {
        if self.reversed {
            self.end - v
        } else {
            self.start + v
        }
    }

    fn travel(&self) -> Scalar {
        if self.reversed {
            -1.0
        } else {
            1.0
        }
    }

    /// The forward speed `1 + b' - a k` and the speed `f` at `v`.
    fn speed(&self, v: Scalar) -> GeomResult<(Scalar, Scalar)> {
        let length = self.length_of_base();
        let s = self.measure(v);
        let u = (v / length).clamp(0.0, 1.0);
        let at = |pair: (Scalar, Scalar)| pair.0 + (pair.1 - pair.0) * u;
        let rate = |pair: (Scalar, Scalar)| (pair.1 - pair.0) / length;
        let k = self.travel() * self.curvature.value(s)?;
        let (a, da) = (at(self.lateral), rate(self.lateral));
        let (long, dlong) = (at(self.longitudinal), rate(self.longitudinal));
        let (vert, dvert) = (at(self.vertical), rate(self.vertical));
        let (b, db) = match self.grade {
            Some(law) => {
                let p = self.travel() * crate::elevation::elevation_grade(law, s)?;
                let q = crate::elevation::elevation_second(law, s)?;
                let cos = 1.0 / p.hypot(1.0);
                let sin = p * cos;
                let turn = q / (1.0 + p * p);
                (
                    long * cos - vert * sin,
                    dlong * cos - long * sin * turn - dvert * sin - vert * cos * turn,
                )
            }
            None => (long, dlong),
        };
        let forward = 1.0 + db - a * k;
        let across = da + b * k;
        let rise = if self.vertical_rate { dvert } else { 0.0 };
        let speed = (forward * forward + across * across + rise * rise).sqrt();
        if !speed.is_finite() {
            return Err(uncertified());
        }
        Ok((forward, speed))
    }

    /// Suprema over the panel `[va, vb]` of the forward speed's first
    /// derivative and of `|f''|` given a lower bound `floor` of the forward
    /// speed; `None` for `|f''|` while `floor` is not positive.
    fn bounds(&self, va: Scalar, vb: Scalar) -> GeomResult<(Scalar, Jet<3>, Jet<3>)> {
        let length = self.length_of_base();
        let (sa, sb) = (self.measure(va), self.measure(vb));
        let (lo, hi) = (sa.min(sb), sa.max(sb));
        let k = self.curvature.jet(lo, hi)?;
        let k: Jet<4> = [k[0], k[1], k[2], 0.0];
        let a = linear(self.lateral, length);
        let long = linear(self.longitudinal, length);
        let vert = linear(self.vertical, length);
        let b = match self.grade {
            Some(law) => {
                let z = crate::elevation::profile_bounds4(law, lo, hi).ok_or_else(uncertified)?;
                let p: Jet<4> = z;
                sum(
                    product(long, composed(COS_OF_ATAN, p)),
                    product(vert, composed(SIN_OF_ATAN, p)),
                )
            }
            None => long,
        };
        let ak = product(a, k);
        let bk = product(b, k);
        let forward = [1.0 + b[1] + ak[0], b[2] + ak[1], b[3] + ak[2]];
        let across = [a[1] + bk[0], a[2] + bk[1], a[3] + bk[2]];
        Ok((forward[1], forward, across))
    }

    /// The certified length between the base distances `va < vb`.
    fn length(&self, va: Scalar, vb: Scalar, tolerance: Scalar) -> GeomResult<DistanceInterval> {
        let span = vb - va;
        let (mut lower, mut upper) = (Compensated::default(), Compensated::default());
        let mut panels = 0_usize;
        let mut stack = vec![(va, vb)];
        while let Some((a, b)) = stack.pop() {
            panels += 1;
            if panels > MAX_CERTIFIED_PANELS {
                return Err(budget());
            }
            let h = b - a;
            let mid = 0.5 * (a + b);
            let (forward_mid, speed_mid) = self.speed(mid)?;
            let (rate, forward, across) = self.bounds(a, b)?;
            let floor = forward_mid - 0.5 * h * rate;
            let error = if floor > 0.0 {
                let g1 = 2.0 * (forward[0] * forward[1] + across[0] * across[1]);
                let g2 = 2.0
                    * (forward[1] * forward[1]
                        + forward[0] * forward[2]
                        + across[1] * across[1]
                        + across[0] * across[2]);
                let gmin = floor * floor;
                let bend = g2 / (2.0 * floor) + g1 * g1 / (4.0 * gmin * floor);
                h * h * h / 24.0 * bend + VALUE_ROUNDING * h * speed_mid
            } else {
                Scalar::INFINITY
            };
            if error.is_finite() && 2.0 * error <= 0.5 * tolerance * h / span {
                lower.add(h * speed_mid - error);
                upper.add(h * speed_mid + error);
                continue;
            }
            if !(mid > a && mid < b) {
                return Err(uncertified());
            }
            stack.push((mid, b));
            stack.push((a, mid));
        }
        finish(lower.value(), upper.value(), tolerance)
    }
}

// --- seams ------------------------------------------------------------------------

/// A seam of a station basis whose frame may jump, with its certified
/// distance (#294).
#[derive(Debug, Clone, Copy, PartialEq)]
pub(crate) struct CertifiedSeam {
    /// Where the kernel reads it: the distance a station on it is read at.
    pub(crate) distance: Scalar,
    /// Where it truly lies.
    pub(crate) interval: DistanceInterval,
}

/// The certified seams of an atomic curve whose frame may jump there,
/// `seams` its seams (`station_seams2`/`station_seams3`): exact ones at
/// their distance, a B-spline's corner knot by chords from its start.
pub(crate) fn atomic_seams(
    curve: StationCurve<'_>,
    seams: Vec<StationSeam>,
    tolerance: Scalar,
) -> GeomResult<Vec<CertifiedSeam>> {
    seams
        .into_iter()
        .filter(|seam| !seam.smooth)
        .map(|seam| {
            let interval = if seam.exact {
                DistanceInterval::point(seam.distance)
            } else {
                curve_span(curve, 0.0, seam.distance, tolerance)?
            };
            Ok(CertifiedSeam {
                distance: seam.distance,
                interval,
            })
        })
        .collect()
}

#[cfg(test)]
mod tests {
    use super::*;

    /// The suprema the leaning section's jets use are upper bounds of
    /// `cos(atan p)`, `sin(atan p)` and their first three derivatives:
    /// sampled densely, and the derivatives checked against differences.
    #[test]
    fn the_grade_suprema_bound_their_functions() {
        let derivatives = |p: Scalar| {
            let u = 1.0 + p * p;
            (
                [
                    u.powf(-0.5),
                    -p * u.powf(-1.5),
                    (2.0 * p * p - 1.0) * u.powf(-2.5),
                    p * (9.0 - 6.0 * p * p) * u.powf(-3.5),
                ],
                [
                    p * u.powf(-0.5),
                    u.powf(-1.5),
                    -3.0 * p * u.powf(-2.5),
                    (12.0 * p * p - 3.0) * u.powf(-3.5),
                ],
            )
        };
        for i in -200_000..=200_000 {
            let p = Scalar::from(i) * 1e-4;
            let (cos, sin) = derivatives(p);
            for k in 0..4 {
                assert!(cos[k].abs() <= COS_OF_ATAN[k], "cos^({k}) at {p}");
                assert!(sin[k].abs() <= SIN_OF_ATAN[k], "sin^({k}) at {p}");
            }
        }
        let h = 1e-5;
        for p in [-3.0, -0.4, 0.0, 0.3, 2.5] {
            let (at, _) = derivatives(p);
            let (up, _) = derivatives(p + h);
            let (down, _) = derivatives(p - h);
            for k in 0..3 {
                let difference = (up[k] - down[k]) / (2.0 * h);
                assert!(
                    (difference - at[k + 1]).abs() < 1e-6,
                    "cos^({}) at {p}",
                    k + 1
                );
            }
        }
    }

    /// A unit vector's derivative bounds hold for a curve whose `|w|` dips.
    #[test]
    fn the_unit_vector_bounds_hold() {
        let w = |t: Scalar| Vec3::new(1.0 + 0.5 * t.cos(), 0.8 * t.sin(), 0.3 * t);
        let dw = |t: Scalar| Vec3::new(-0.5 * t.sin(), 0.8 * t.cos(), 0.3);
        let (w1, w2) = (1.0, 1.0);
        let r = 0.45;
        let (u1, u2) = unit_bounds(r, w1, w2);
        let h = 1e-4;
        for i in 0..2000 {
            let t = Scalar::from(i) * 3e-3;
            let (_, d) = unit_rate(w(t), dw(t)).unwrap();
            assert!(d.length() <= u1);
            let (_, up) = unit_rate(w(t + h), dw(t + h)).unwrap();
            let (_, down) = unit_rate(w(t - h), dw(t - h)).unwrap();
            assert!(((up - down) / (2.0 * h)).length() <= u2);
        }
    }

    #[test]
    fn a_compensated_sum_of_many_terms_stays_exact() {
        let mut sum = Compensated::default();
        for _ in 0..1_000_000 {
            sum.add(0.1);
        }
        assert!((sum.value() - 100_000.0).abs() <= 1e-9);
    }
}
