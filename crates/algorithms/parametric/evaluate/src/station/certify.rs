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
//! # Chords and their deficit
//!
//! Where the measured curve `P(t)` is evaluated in closed form at its
//! native parameter -- a line, a circle, an ellipse or a B-spline, and any
//! offset of one, 2D or 3D, measured by arc length -- each panel `[a, b]`,
//! `h = b - a`, is bounded by its chord `c = |P(b) - P(a)|` and certified
//! bounds `M2`, `M3` on `|P''|`, `|P'''|` there:
//!
//! - `c` is a lower bound of the arc length (a chord is never longer);
//! - `|P'(v)| >= m = |P'(mid)| - M2 h / 2` everywhere on the panel, so the
//!   tangent turns by at most `beta = h M2 / (2 m)` from the midpoint's and
//!   the arc length is at most `c / cos(beta)`;
//! - the deficit `L - c` is bounded from both sides through the velocity
//!   across the chord, a zero-mean function whose square integral `M2`
//!   bounds as `h^3`, and `M3` with the panel's own velocities to within
//!   `h^4` (#307; see [`panels`]), so a panel's width falls as `h^4`.
//!
//! The panels are refined adaptively: the widest is halved first, until
//! their widths sum to an eighth of the tolerance ([`panels::adapt`]), so
//! an ordinary curve certifies to `1e-8` relative and below within the
//! budget. `M2` and `M3` come from the certified derivative suprema of
//! [`crate::bound`] (`|B''|` to `|B''''|` of the base) and, for an offset,
//! the frame's: a unit vector `w / |w|` has `|u'| <= |w'| / r`, `|u''| <=
//! |w''| / r + 3 |w'|^2 / r^2` and `|u'''| <= |w'''| / r + 9 |w'| |w''| /
//! r^2 + 15 |w'|^3 / r^3` where `r` bounds `|w|` from below (`|B'|` for the
//! tangent, `|B'_xy|` for the reference-up lateral, `|V x B'|` for a 3D
//! offset's direction), and Leibniz for each law component times its axis.
//! Each point carries a rounding error of [`POINT_ROUNDING`] relative to
//! its largest coordinate, which widens the chain of chords.
//!
//! An offset by distances along an ellipse or a B-spline (#307) reads its
//! law in the base's arc length, a quadrature: the base is measured by
//! chords first, and the law at a native parameter is read through that
//! certified length (see [`beside`]), the measured curve being the offset
//! whose law runs linearly in the base's true arc length between the
//! points the kernel reads as the piece's ends.
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
//! a k) > 0` (a cusp is refused). Its panels are refined adaptively too.
//!
//! # Refusals
//!
//! By a typed [`GeomError::UnsupportedInput`] naming
//! [`UNCERTIFIED_LENGTH`]: an offset of a banked curve, of a 3D intrinsic
//! curve, of a chain's parametric piece or of an elevated curve with an
//! intrinsic profile in a leaning frame; a 3D offset direction that is not
//! vertical beside a plan-parameterised base; a base placed in a frame
//! that tilts `+Z` there; a skewed frame; and a tolerance the rounding of
//! the points already exceeds. A tolerance the budget of
//! [`MAX_CERTIFIED_PANELS`] panels does not reach is
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
use super::{station_length2, station_length3, StationSeam, KEEPS_UP_TOLERANCE};
use crate::arc_parameter::{
    parameter_at_arc_length2, parameter_at_arc_length3, ARC_LENGTH_TOLERANCE,
};
use crate::bound::{
    continuity_breaks2, continuity_breaks3, curve_derivative_bounds2, curve_derivative_bounds3,
    curve_fourth_bound2, curve_fourth_bound3, CurveDerivativeBounds,
};
use crate::curve::{
    derivative2, derivative3, domain2, domain3, evaluate2, evaluate3, second_derivative2,
    second_derivative3,
};
use beside::chord_offset_length;
use panels::{adapt, finish, product, Compensated, Jet, Measured, PANEL_SHARE, VELOCITY_ROUNDING};

mod beside;
mod panels;

/// The input a length no certified bound here reaches is refused for.
pub const UNCERTIFIED_LENGTH: &str =
    "a certified length of a station piece no bound here reaches: an offset of a banked curve, a \
     3D intrinsic curve, a chain's parametric piece or an elevated curve with an intrinsic \
     profile in a leaning frame; a 3D offset direction that is not vertical beside a \
     plan-measured or intrinsic base; a base placed in a frame that tilts +Z; a skewed frame; or \
     a tolerance below the points' rounding";

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
    // `P''` continuous on every panel, for its third-derivative bound.
    let mut breaks = native.breaks(2);
    breaks.extend(native.breaks(1));
    native.measured(|measured| measured.chord_turning(t0.min(t1), t0.max(t1), &breaks, tolerance))
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

    /// A certified supremum of `|B''''|` over `[a, b]`, where one is known.
    fn fourth(&self, a: Scalar, b: Scalar) -> Option<Scalar> {
        match self {
            Self::Two(curve) => curve_fourth_bound2(curve, a, b),
            Self::Three(curve) => curve_fourth_bound3(curve, a, b),
        }
        .filter(|d| d.is_finite())
    }

    /// `run` given the curve measured by chords: its points exact but for
    /// rounding, `P''` and `P'''` bounded by its own suprema.
    fn measured<R>(&self, run: impl FnOnce(&Measured<'_>) -> GeomResult<R>) -> GeomResult<R> {
        let point = |t: Scalar| Ok((self.point(t)?, 0.0));
        let velocity = |t: Scalar| {
            let d = self.derivative(t)?;
            Ok((d, VELOCITY_ROUNDING * d.length()))
        };
        let bounds = |a: Scalar, b: Scalar| {
            let d = self.bounds(a, b)?;
            Ok((d.second, Some(d.third)))
        };
        run(&Measured {
            point: &point,
            velocity: &velocity,
            bounds: &bounds,
            scale: self.extent(),
        })
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

// --- speeds in closed form ------------------------------------------------------

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

    /// Suprema over the panel `[va, vb]` of the forward speed `A = 1 + b' -
    /// a k` and the cross speed `B = a' + b k` and their first two
    /// derivatives, and `sup |A'|` first, for the floor of `A`.
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

    /// The midpoint rule over the base distances `[a, b]`: the panel's
    /// bounds `(a, b, lower, upper)`, unbounded where its forward speed
    /// may vanish.
    fn panel(&self, a: Scalar, b: Scalar) -> GeomResult<(Scalar, Scalar, Scalar, Scalar)> {
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
        if !error.is_finite() {
            return Ok((a, b, Scalar::NEG_INFINITY, Scalar::INFINITY));
        }
        Ok((a, b, h * speed_mid - error, h * speed_mid + error))
    }

    /// The certified length between the base distances `va < vb`, its
    /// panels refined where they are widest ([`panels::adapt`]).
    fn length(&self, va: Scalar, vb: Scalar, tolerance: Scalar) -> GeomResult<DistanceInterval> {
        let refined = adapt(
            vec![self.panel(va, vb)?],
            PANEL_SHARE * tolerance,
            |p| p.3 - p.2,
            |&(a, b, ..)| {
                let mid = 0.5 * (a + b);
                if !(mid > a && mid < b) {
                    return Err(uncertified());
                }
                Ok((self.panel(a, mid)?, self.panel(mid, b)?))
            },
        )?;
        let (mut lower, mut upper) = (Compensated::default(), Compensated::default());
        for &(_, _, lo, hi) in &refined.panels {
            lower.add(lo);
            upper.add(hi);
        }
        finish(lower.value(), upper.value(), tolerance).map_err(|error| {
            if refined.exhausted {
                budget()
            } else {
                error
            }
        })
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
}
