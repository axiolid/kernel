//! Certified bounds on a section frame law's derivatives in station
//! distance (#247, ADR 0082 amendment).
//!
//! A run of sections placed at stations ([the parent module](super)) puts
//! a section point at `A(s) = O(s) + x(s) L(s) + y(s) U(s) + z(s) T(s)`:
//! the basis point `O` and the section frame (tangent `T`, lateral `L`,
//! up `U`) at the distance `s`, with local coordinates linear in `s`
//! between two stations. A chord of `A` is within `h^2 / 8 sup |A''|` of
//! it (`crate::bound`), and
//!
//! `|A''| <= |O''| + |x| |L''| + |y| |U''| + |z| |T''| + 2 (|x'| |L'| + |y'| |U'| + |z'| |T'|)`,
//!
//! so a run is bounded once `|O'|, |O''|` and the first two derivatives of
//! each axis are, over a span of distance. [`piece_frame_bounds`] returns
//! them for a span inside one [`StationPiece`] that crosses no seam of its
//! curve, or `None` where they are not certified.
//!
//! # Where the bounds come from
//!
//! Everything is derived in the piece's own station measure `s` from the
//! curve's tangent law: `k1 >= sup |T'|`, `k2 >= sup |T''|`.
//!
//! - A line and a polyline's segment: `T` constant, `k1 = k2 = 0`.
//! - A circle of radius `r` with an orthonormal frame, by arc length:
//!   `|O''| = k1 = 1 / r`, `k2 = 1 / r^2`.
//! - An elevated curve (plan a line or a circle, profile of closed-form
//!   pieces) and a banked curve, by plan distance, which is their
//!   parameter: with `D_j >= sup |c^(j)|` (`crate::elevated`) and `sigma`
//!   a lower bound of `|c'|` over the span (`|c'(mid)| - h/2 D_2`), the
//!   unit tangent has `k1 = D_2 / sigma` and `k2 = 2 D_3 / sigma + 6 D_2^2
//!   / sigma^2`, the bounds a swept disk's wall rests on (#252).
//!
//! The axes follow from `T` through one rule: a unit vector `w / |w|` with
//! `|w| >= rho` has `|u'| <= |w'| / rho` and `|u''| <= 2 |w''| / rho + 6
//! |w'|^2 / rho^2` (`super::certify`). `rho` is the horizontal length of
//! `T` at the middle of the span less `h/2 k1`, a lower bound over it.
//!
//! - A 2D curve's planar frame: `L` is `T` turned a quarter, `U = +Z`.
//! - The reference-up frame: `R = (T x Z) / |T x Z|` with `|T x Z| =
//!   |T_xy| >= rho` gives `m1 = k1 / rho`, `m2 = 2 k2 / rho + 6 k1^2 /
//!   rho^2`; `L = -R`, and `U = R x T` has `|U'| <= m1 + k1`, `|U''| <=
//!   m2 + 2 m1 k1 + k2`.
//! - A banked curve's rolled section: `n` the horizontal left normal
//!   (bounds `m1`, `m2` as `R`), `u = T x n` (bounds `k1 + m1`, `k2 + 2 k1
//!   m1 + m2`), `L = cos(phi) n + sin(phi) u` and `U = -sin(phi) n +
//!   cos(phi) u` for the roll `phi`. With `a_j = hypot(|n^(j)|, |u^(j)|)`
//!   (Cauchy-Schwarz over `cos`, `sin`) and `p_j >= sup |phi^(j)|`, both
//!   have `|X'| <= p1 + a1` and `|X''| <= p2 + p1^2 + 2 p1 a1 + a2`. The
//!   roll is `asin(g)`, `g = D / b` (or `D sqrt(1 + G^2) / b` under the
//!   vertical-rise convention, `G` the grade), so `|phi'| <= |g'| / c` and
//!   `|phi''| <= |g''| / c + |g| |g'|^2 / c^3` with `c = sqrt(1 - sup
//!   |g|^2)`, from the cant law's certified derivative suprema
//!   (`crate::banked`).
//! - The plan reading ([`FrameReading::Plan`]): `T_p = T_xy / |T_xy|` has
//!   the bounds `m1`, `m2`, `L_p` its quarter turn, `U = +Z`.
//!
//! A reversed piece negates `T` and `L`, which keeps every magnitude; a
//! placement that keeps `+Z` moves the frame rigidly (its stretch, a hair
//! above one, scales the bounds); one that tilts `+Z` reads the
//! reference-up frame of the moved tangent ([`SectionFrame::carried`]),
//! bounded as above with `rho` from the moved tangent, and a plan-measured
//! piece so placed is not bounded.
//!
//! # The measure
//!
//! [`FrameLawBounds::measure`] bounds how far the distance a frame is read
//! at may lie from the one asked for. A line, an elevated and a banked
//! curve read their parameter from the distance in closed form: `0`. A
//! circle with an orthonormal frame and a polyline read it by the
//! arc-length inverse ([`crate::arc_parameter`]); their speed is constant
//! on every panel (a polyline's panels split at its vertices), so the
//! Gauss-Legendre rule is exact up to rounding and the root find stops
//! within its tolerance: twice [`ARC_LENGTH_TOLERANCE`]` * max(1, s)`.
//!
//! # The signed bend (#361)
//!
//! The derivative suprema bound how far a section moves, not whether the
//! section planes keep advancing. [`FrameLawBounds::forward`] bounds the
//! forward rate `O'.N` of the basis point from below (one in arc length and
//! in plan distance, the horizontal part of the tangent under the plan
//! reading of an arc-length space curve, times a placement's smallest
//! column length), and [`FrameLaw::bend`] bounds it together with the
//! tangent's turning toward the lateral and the up axis, `L.N'` and `U.N'`:
//! a section point `(a, b)` advances at `O'.N - a L.N' - b U.N'` plus its
//! longitudinal rate, and the planes fold where that is not positive. The
//! `bend` module states the closed forms.
//!
//! # Read by quadrature (#316)
//!
//! An ellipse, a B-spline and an offset beside a line, a circle, an
//! ellipse or a B-spline read their frames through estimates (the
//! arc-length inverse, an offset's own length table, an offset's tangent by
//! a difference quotient), so [`piece_frame_bounds`] leaves them unbounded.
//! [`FrameLaw`] bounds them against their EXACT frame law instead: the
//! curve the length certificates of [`super::certify`] measure, framed by
//! its own tangent at its true arc length from the piece's start. Each read
//! frame stands at a native parameter where the exact curve's point and
//! velocity are known in closed form to a certified error; the read
//! point's and each read axis's distance from the exact frame there are
//! computed and carried as [`FrameLawBounds::point_error`] and
//! [`FrameLawBounds::axis_error`], and the true arc length to that
//! parameter, bounded by chord panels over the piece measured once, gives
//! the measure. Over a span between two reads the exact law is bounded in
//! its own arc length from the certified suprema of `|P''|` and `|P'''|`
//! (the base's `|B''|` to `|B''''|` beside an offset, and the frame's unit
//! vectors) and a floor of `|P'|`: `|O''| = |T'| <= M2 / sigma^2` and
//! `|T''| <= (M3 / sigma + 3 M2^2 / sigma^2) / sigma^2 + M2^2 / sigma^4`.
//! The axes follow from the tangent as above. On a composite, a piece after
//! one measured by quadrature starts where the kernel sums the lengths
//! before it; its distance from the certified joint
//! ([`super::CompositeBasis::joint_intervals`]) adds to the measure.
//!
//! # Intrinsic curves (#361)
//!
//! An intrinsic curve's heading is the closed-form integral of its
//! curvature law, so over a seam-free span `|T'| <= sup |k|` and `|T''| <=
//! sup |k'| + sup k^2` in arc length, from the law's own suprema
//! (`crate::bound::intrinsic_derivative_bounds2`), and the axes follow by
//! the rules above. An elevated curve over an intrinsic plan composes them
//! with the profile in plan distance, as over a line or a circle. Only the
//! POINT is a quadrature, with no error bound of its own, so
//! [`piece_frame_bounds`] leaves these curves unbounded and [`FrameLaw`]
//! compares every frame it reads with the exact one, as for the curves
//! read by quadrature: the exact point from a certified
//! [`IntrinsicPositionTable`](crate::bound::IntrinsicPositionTable), the
//! exact axes from the closed-form heading. The measure is exact. The
//! families are a 2D intrinsic curve, a level planar intrinsic space curve
//! and an elevated curve over an intrinsic plan with a closed-form
//! profile, each with an orthonormal start frame and any curvature law
//! with closed-form bounds of `k` to `k'''`.
//!
//! # One station (#405)
//!
//! A part placed at ONE station rests on one frame, not a run of them:
//! [`FrameLaw::error_at`] charges the frame read at a distance against the
//! exact frame law's there, its point and each axis ([`StationFrameError`]).
//!
//! The frame is read as [`StationPiece::section_on`] reads it: the piece a
//! composite names at `distance` from `side`, a seam of the curve inside it
//! read from the side named (#263). Two charges make the error:
//!
//! - **The read** (`e_p`, `e_a`). A slot read in closed form (a circle, a
//!   polyline, an elevated or banked curve over a line or a circle) has no
//!   error of its own but the rounding of the read: [`POINT_ROUNDING`](super::POINT_ROUNDING)
//!   of the point's largest coordinate (one at least), and `1e-14` per axis. An intrinsic curve's slot (a 2D curve, a level planar space
//!   curve, an elevated curve over an intrinsic plan) compares the read
//!   frame with the exact one, as [`FrameLaw::bounds`] does for a run: both
//!   sides of a seam of its law are compared, which covers the one named.
//! - **The distance's resolution** (`delta`, the bounds' `measure`). Where
//!   the distance is resolved by the arc-length inverse (a circle, a
//!   polyline) the read stands within `delta` of the distance asked for,
//!   and on a composite a piece after pieces measured by quadrature adds
//!   its start's certified distance from the true joint. The exact frame at
//!   the distance asked for is then within `delta |O'|` of the read point
//!   and `delta max(|T'|, |L'|, |U'|)` of each read axis, the suprema taken
//!   over `[u - 2 delta, u + 2 delta]` on the side read: where the station
//!   is ON a seam it is clipped to the side named. A curve whose distance is
//!   exact (a line, an elevated or banked curve, an intrinsic curve) has
//!   `delta = 0`.
//!
//! Not charged, `None`: a slot read by quadrature or unbounded, a composite
//! piece placed inexactly, a seam or a joint within `2 delta` of the
//! station that it is not on (the exact frame might stand on its other
//! side), and a result that is not finite.
//!
//! # Not bounded
//!
//! `None`, never a guess: a chain, and an offset beside one of them, beside
//! an intrinsic, an elevated or a banked curve; a banked curve over an
//! intrinsic or chain plan; an intrinsic space curve with torsion or out of
//! a level plane; an intrinsic curve with a stretched or skewed start
//! frame; a circle with a skewed frame; an elevated or banked curve with
//! another plan or an intrinsic profile; a span that crosses a seam, a vertex of a polyline
//! or a piece of a law; a tangent too steep for a horizontal lower bound;
//! a cant whose roll reaches a quarter turn; a curve read by quadrature
//! placed in a frame that tilts `+Z`, and a placed offset piece.
//! [`piece_frame_bounds`] alone also leaves an ellipse, a B-spline and
//! every offset unbounded: they are [`FrameLaw`]'s.
//!
//! [`SectionFrame::carried`]: super::SectionFrame::carried

use std::cell::OnceCell;

use axiolid_contracts::GeomResult;
use axiolid_core::{Scalar, Transform3, Vec2, Vec3};
use axiolid_curve::{Curve2, Curve3};

use super::certify::{checked_tolerance, reads_intrinsic, IntrinsicReader, PieceReader};
use super::offset::OffsetTable;
use super::{CompositeBasis, DistanceConvention, StationCurve, StationPiece, KEEPS_UP_TOLERANCE};
use crate::arc_parameter::ARC_LENGTH_TOLERANCE;
use crate::bound::ROUNDING;

mod bend;
mod error;

pub use bend::BendBounds;
pub use error::StationFrameError;

/// Which section frame a run places its sections in.
#[non_exhaustive]
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum FrameReading {
    /// The curve's section frame ([`super::SectionFrame`]).
    Section,
    /// The plan frame ([`super::SectionFrame::plan`]).
    Plan,
}

/// Certified suprema, over a span of station distance, of the derivatives
/// of a section frame law in that distance; see the
/// [module documentation](self).
#[non_exhaustive]
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct FrameLawBounds {
    /// `sup |O'|` and `sup |O''|` of the basis point.
    pub point: [Scalar; 2],
    /// `sup |T'|` and `sup |T''|` of the unit tangent.
    pub tangent: [Scalar; 2],
    /// `sup |L'|` and `sup |L''|` of the lateral axis.
    pub lateral: [Scalar; 2],
    /// `sup |U'|` and `sup |U''|` of the up axis.
    pub up: [Scalar; 2],
    /// How far the distance a frame is read at may lie from the distance
    /// asked for (see [the measure](self#the-measure)).
    pub measure: Scalar,
    /// How far a frame's point read at a distance may lie from the exact
    /// frame law's point at the distance it truly stands at: zero where
    /// the frame is read in closed form, a certified bound where it is
    /// read by quadrature ([`FrameLaw`], #316).
    pub point_error: Scalar,
    /// How far each of a read frame's axes may lie from the exact frame
    /// law's, likewise.
    pub axis_error: Scalar,
    /// A certified lower bound over the span of `O'.N`, `N` the frame's
    /// tangent: one in arc length and in plan distance, the horizontal
    /// floor `rho` under the plan reading of an arc-length space curve,
    /// times a placement's smallest column length, less rounding (#361).
    pub forward: Scalar,
}

impl FrameLawBounds {
    /// Whether the frame is constant and the point moves on a straight
    /// line: every second derivative, and every axis's first, is zero.
    #[must_use]
    pub fn is_straight(&self) -> bool {
        self.point[1] == 0.0
            && self.tangent == [0.0; 2]
            && self.lateral == [0.0; 2]
            && self.up == [0.0; 2]
    }
}

/// How a curve's frame is built from its tangent.
#[derive(Debug, Clone, Copy)]
enum Kind {
    /// A 2D curve's frame: the tangent, its quarter turn, `+Z`.
    Planar,
    /// The reference-up frame against `+Z`.
    ReferenceUp,
    /// A banked curve's rolled section, with the roll's derivative bounds.
    Banked([Scalar; 2]),
}

/// A curve's tangent law over a span.
#[derive(Debug, Clone, Copy)]
struct TangentLaw {
    point: [Scalar; 2],
    tangent: [Scalar; 2],
    /// The unit tangent at the middle of the span.
    middle: Vec3,
    kind: Kind,
    measure: Scalar,
}

/// The measure slack of a distance read by the arc-length inverse on a
/// speed constant per panel.
fn numerical_measure(hi: Scalar) -> Scalar {
    2.0 * ARC_LENGTH_TOLERANCE * hi.abs().max(1.0)
}

fn unit(v: Vec3) -> Option<Vec3> {
    let length = v.length();
    (length.is_finite() && length > 0.0).then(|| v / length)
}

fn lift(v: Vec2) -> Vec3 {
    Vec3::new(v.x, v.y, 0.0)
}

/// Whether a conic frame's axes are orthonormal to rounding, so a circle's
/// speed is its radius everywhere.
fn orthonormal(x: Vec3, y: Vec3) -> bool {
    let near = 1e-12;
    (x.length() - 1.0).abs() <= near && (y.length() - 1.0).abs() <= near && x.dot(y).abs() <= near
}

fn straight(direction: Vec3, kind: Kind, measure: Scalar) -> Option<TangentLaw> {
    Some(TangentLaw {
        point: [1.0, 0.0],
        tangent: [0.0, 0.0],
        middle: unit(direction)?,
        kind,
        measure,
    })
}

fn circle(radius: Scalar, middle: Vec3, kind: Kind, hi: Scalar) -> Option<TangentLaw> {
    if !(radius.is_finite() && radius > 0.0) {
        return None;
    }
    let k = 1.0 / radius;
    Some(TangentLaw {
        point: [1.0, k],
        tangent: [k, k * k],
        middle: unit(middle)?,
        kind,
        measure: numerical_measure(hi),
    })
}

/// The segment of a polyline holding the measure span `[lo, hi]`, as its
/// direction; `None` across a vertex.
fn polyline_segment<P: Copy>(
    points: &[P],
    closed: bool,
    lo: Scalar,
    hi: Scalar,
    difference: impl Fn(P, P) -> Vec3,
) -> Option<Vec3> {
    let count = points.len();
    let spans: Vec<(usize, usize)> = if closed {
        (0..count).map(|i| (i, (i + 1) % count)).collect()
    } else {
        (0..count.saturating_sub(1)).map(|i| (i, i + 1)).collect()
    };
    let total: Scalar = spans
        .iter()
        .map(|&(i, j)| difference(points[i], points[j]).length())
        .sum();
    let slack = numerical_measure(total);
    let mut start = 0.0;
    for (i, j) in spans {
        let step = difference(points[i], points[j]);
        let end = start + step.length();
        if step.length() > 0.0 && lo >= start - slack && hi <= end + slack {
            return Some(step);
        }
        start = end;
    }
    None
}

fn tangent_law2(curve: &Curve2, lo: Scalar, hi: Scalar) -> Option<TangentLaw> {
    match curve {
        Curve2::Line(line) => straight(lift(line.direction), Kind::Planar, 0.0),
        Curve2::Circle(c) => {
            if !orthonormal(lift(c.frame.x), lift(c.frame.y)) {
                return None;
            }
            let middle = crate::curve::derivative2(curve, 0.5 * (lo + hi) / c.radius).ok()?;
            circle(c.radius, lift(middle), Kind::Planar, hi)
        }
        Curve2::Polyline(polyline) => {
            let step = polyline_segment(&polyline.points, polyline.closed, lo, hi, |a, b| {
                lift(b - a)
            })?;
            straight(step, Kind::Planar, numerical_measure(hi))
        }
        _ => None,
    }
}

/// An elevated or banked curve's point law in plan distance, from its
/// certified derivative suprema and its derivative at the middle.
fn plan_measured(
    bounds: crate::bound::CurveDerivativeBounds,
    middle: Vec3,
    lo: Scalar,
    hi: Scalar,
    kind: Kind,
) -> Option<TangentLaw> {
    let (d2, d3) = (bounds.second, bounds.third);
    let sigma = middle.length() - 0.5 * (hi - lo) * d2;
    if !(sigma.is_finite() && sigma > 0.0) {
        return None;
    }
    Some(TangentLaw {
        point: [bounds.first, d2],
        tangent: [
            d2 / sigma,
            2.0 * d3 / sigma + 6.0 * d2 * d2 / (sigma * sigma),
        ],
        middle: unit(middle)?,
        kind,
        measure: 0.0,
    })
}

/// Whether an elevated curve's point is closed-form with certified
/// derivatives: a line or circle plan and a profile of closed-form pieces.
fn closed_form_elevated(curve: &axiolid_curve::Elevated3) -> bool {
    let plan = match curve.plan.as_ref() {
        Curve2::Line(_) => true,
        Curve2::Circle(c) => orthonormal(lift(c.frame.x), lift(c.frame.y)),
        _ => false,
    };
    plan && crate::elevated::certifies_elevated(curve)
}

fn tangent_law3(curve: &Curve3, lo: Scalar, hi: Scalar) -> Option<TangentLaw> {
    let mid = 0.5 * (lo + hi);
    match curve {
        Curve3::Line(line) => straight(line.direction, Kind::ReferenceUp, 0.0),
        Curve3::Circle(c) => {
            if !orthonormal(c.frame.x, c.frame.y) {
                return None;
            }
            let middle = crate::curve::derivative3(curve, mid / c.radius).ok()?;
            circle(c.radius, middle, Kind::ReferenceUp, hi)
        }
        Curve3::Polyline(polyline) => {
            let step = polyline_segment(&polyline.points, polyline.closed, lo, hi, |a, b| b - a)?;
            straight(step, Kind::ReferenceUp, numerical_measure(hi))
        }
        Curve3::Elevated(elevated) => {
            if !closed_form_elevated(elevated) {
                return None;
            }
            let bounds = crate::elevated::elevated_derivative_bounds(elevated, lo, hi)?;
            let middle = crate::elevated::elevated_derivative(elevated, mid).ok()?;
            plan_measured(bounds, middle, lo, hi, Kind::ReferenceUp)
        }
        Curve3::Banked(banked) => {
            if !closed_form_elevated(&banked.base) {
                return None;
            }
            let bounds = crate::elevated::banked_derivative_bounds(banked, lo, hi)?;
            let roll = crate::banked::roll_bounds(banked, lo, hi)?;
            let middle = crate::banked::banked_derivative(banked, mid).ok()?;
            plan_measured(bounds, middle, lo, hi, Kind::Banked(roll))
        }
        _ => None,
    }
}

/// `sup |u'|` and `sup |u''|` of a unit vector `w / |w|` with `|w'| <=
/// w1`, `|w''| <= w2` and `|w| >= rho`.
fn normalised(w1: Scalar, w2: Scalar, rho: Scalar) -> [Scalar; 2] {
    [w1 / rho, 2.0 * w2 / rho + 6.0 * w1 * w1 / (rho * rho)]
}

/// The largest factor a placement's linear part stretches a length by, or
/// `None` for one whose columns are not orthogonal to rounding.
fn placement_stretch(rigid: axiolid_core::Transform3) -> Option<Scalar> {
    let m = rigid.matrix3;
    let (sx, sy, sz) = (m.x_axis.length(), m.y_axis.length(), m.z_axis.length());
    let eps = 32.0 * f64::EPSILON;
    let orthogonal = m.x_axis.dot(m.y_axis).abs() <= eps * sx * sy
        && m.x_axis.dot(m.z_axis).abs() <= eps * sx * sz
        && m.y_axis.dot(m.z_axis).abs() <= eps * sy * sz;
    orthogonal.then(|| sx.max(sy).max(sz) * (1.0 + 3.0 * eps))
}

/// The smallest factor a placement's linear part shrinks a length by; call
/// it only where [`placement_stretch`] holds (orthogonal columns).
fn placement_shrink(rigid: axiolid_core::Transform3) -> Scalar {
    let m = rigid.matrix3;
    let (sx, sy, sz) = (m.x_axis.length(), m.y_axis.length(), m.z_axis.length());
    sx.min(sy).min(sz) * (1.0 - 3.0 * 32.0 * f64::EPSILON)
}

/// Certified bounds on the derivatives, in station distance, of the
/// section frame `piece` presents over `[from, to]` (distances into the
/// piece, in the composite's direction), read as `reading`; see the
/// [module documentation](self).
///
/// The span must lie inside the piece and cross no seam of its curve
/// (`StationPiece`'s seams, a polyline's vertices, a piece of an
/// elevation, cant or pivot law); a caller splits there.
///
/// `None` where the frame law is not certified: an ellipse, a B-spline,
/// an intrinsic curve (only [`FrameLaw`] certifies those, #361), a chain,
/// an offset, a skewed circle, an elevated or
/// banked curve whose plan is not a line or a circle or whose profile is
/// not of closed-form pieces, a span across a seam, a tangent too steep to
/// bound its horizontal part, a roll reaching a quarter turn, and a
/// plan-measured piece under a placement that tilts `+Z`.
#[must_use]
pub fn piece_frame_bounds(
    piece: &StationPiece<'_>,
    from: Scalar,
    to: Scalar,
    reading: FrameReading,
) -> Option<FrameLawBounds> {
    let (u0, u1) = (from.min(to), from.max(to));
    let length = piece.length();
    let slack = numerical_measure(length);
    if !(u0.is_finite() && u1.is_finite()) || u0 < -slack || u1 > length + slack {
        return None;
    }
    let (u0, u1) = (u0.max(0.0), u1.min(length));
    let (m0, m1) = if piece.reversed {
        (piece.end - u1, piece.end - u0)
    } else {
        (piece.start + u0, piece.start + u1)
    };
    let law = match piece.curve {
        StationCurve::Two(curve) => tangent_law2(curve, m0, m1)?,
        StationCurve::Three(curve) => tangent_law3(curve, m0, m1)?,
        StationCurve::Offset(_) => return None,
    };
    finish(piece, law, m1 - m0, reading, [0.0; 2])
}

/// The bounds of a piece read by quadrature (#316) over `[from, to]` into
/// it, from `reader`'s certified frame law; `None` where it is not
/// bounded. An offset piece reads through `table`.
pub(crate) fn read_frame_bounds(
    reader: &PieceReader<'_>,
    piece: &StationPiece<'_>,
    from: Scalar,
    to: Scalar,
    table: Option<&OffsetTable>,
    reading: FrameReading,
) -> Option<FrameLawBounds> {
    let (u0, u1) = (from.min(to), from.max(to));
    let length = piece.length();
    let slack = numerical_measure(length);
    if !(u0.is_finite() && u1.is_finite()) || u0 < -slack || u1 > length + slack {
        return None;
    }
    let read = reader.read(piece, u0, u1, table).ok()??;
    let law = TangentLaw {
        point: read.point,
        tangent: read.tangent,
        middle: unit(read.middle)?,
        kind: if read.planar {
            Kind::Planar
        } else {
            Kind::ReferenceUp
        },
        measure: read.measure,
    };
    finish(
        piece,
        law,
        u1 - u0 + 2.0 * read.measure,
        reading,
        [read.point_error, read.axis_error],
    )
}

/// The bounds of `law` over a span `span` long, as `piece` presents it
/// (reversed, placed) and read as `reading`, with a read point's and a
/// read axis's errors `errors` (#316).
fn finish(
    piece: &StationPiece<'_>,
    mut law: TangentLaw,
    span: Scalar,
    reading: FrameReading,
    errors: [Scalar; 2],
) -> Option<FrameLawBounds> {
    if piece.reversed {
        law.middle = -law.middle;
    }
    let mut scale = 1.0;
    let mut shrink = 1.0;
    if let Some(rigid) = piece.placement {
        scale = placement_stretch(rigid)?;
        shrink = placement_shrink(rigid);
        law.middle = unit(rigid.transform_vector3(law.middle))?;
        if (rigid.transform_vector3(Vec3::Z) - Vec3::Z).length() > KEEPS_UP_TOLERANCE {
            if piece.convention() == DistanceConvention::PlanDistance {
                return None;
            }
            law.kind = Kind::ReferenceUp;
        }
    }
    let [k1, k2] = law.tangent;
    let horizontal = || -> Option<Scalar> {
        let rho = law.middle.x.hypot(law.middle.y) - 0.5 * span * k1;
        (rho.is_finite() && rho > 0.0).then_some(rho)
    };
    let tangent = [k1, k2];
    let (lateral, up) = match law.kind {
        Kind::Planar => (tangent, [0.0, 0.0]),
        Kind::ReferenceUp => {
            let [r1, r2] = normalised(k1, k2, horizontal()?);
            ([r1, r2], [r1 + k1, r2 + 2.0 * r1 * k1 + k2])
        }
        Kind::Banked([p1, p2]) => {
            let [n1, n2] = normalised(k1, k2, horizontal()?);
            let (s1, s2) = (k1 + n1, k2 + 2.0 * k1 * n1 + n2);
            let (a1, a2) = (n1.hypot(s1), n2.hypot(s2));
            let axis = [p1 + a1, p2 + p1 * p1 + 2.0 * p1 * a1 + a2];
            (axis, axis)
        }
    };
    let [point_error, mut axis_error] = errors;
    let mut forward = 1.0;
    let (tangent, lateral, up) = match reading {
        FrameReading::Section => (tangent, lateral, up),
        // A 2D curve's frame is already its plan frame.
        FrameReading::Plan if matches!(law.kind, Kind::Planar) => (tangent, lateral, up),
        FrameReading::Plan => {
            let rho = horizontal()?;
            // A tangent `e` off has its plan direction `2 e / rho` off.
            if axis_error > 0.0 {
                axis_error = 2.0 * axis_error / rho;
            }
            let flat = normalised(k1, k2, rho);
            // `O'` is the horizontal unit vector in plan distance, and the
            // tangent's horizontal part in arc length.
            if piece.convention() != DistanceConvention::PlanDistance {
                forward = rho;
            }
            (flat, flat, [0.0, 0.0])
        }
    };
    let grow = |v: [Scalar; 2]| v.map(|x| x * scale * (1.0 + ROUNDING));
    let out = FrameLawBounds {
        point: grow(law.point),
        tangent: grow(tangent),
        lateral: grow(lateral),
        up: grow(up),
        measure: law.measure,
        point_error: point_error * scale * (1.0 + ROUNDING),
        axis_error: axis_error * scale * (1.0 + ROUNDING),
        forward: forward * shrink * (1.0 - ROUNDING),
    };
    let finite = [out.point, out.tangent, out.lateral, out.up]
        .iter()
        .flatten()
        .chain([
            &out.measure,
            &out.point_error,
            &out.axis_error,
            &out.forward,
        ])
        .all(|v| v.is_finite());
    finite.then_some(out)
}

/// The tangent law of a piece of an intrinsic curve over `[lo, hi]` of its
/// measure (#361): the curvature law's closed-form suprema, `|T'| <= sup
/// |k|` and `|T''| <= sup |k'| + sup k^2` in arc length, composed with the
/// profile in plan distance on an elevated curve (the line and circle
/// plans' path). `None` across a seam of the law, where `k'` may jump: the
/// caller splits.
fn intrinsic_tangent_law(
    reader: &IntrinsicReader<'_>,
    lo: Scalar,
    hi: Scalar,
) -> Option<TangentLaw> {
    let (bounds, middle) = reader.derivatives(lo, hi)?;
    if reader.is_elevated() {
        return plan_measured(bounds, middle, lo, hi, Kind::ReferenceUp);
    }
    if !(bounds.third.is_finite() && bounds.second.is_finite()) {
        return None;
    }
    Some(TangentLaw {
        point: [bounds.first, bounds.second],
        tangent: [bounds.second, bounds.third],
        middle: unit(middle)?,
        kind: Kind::Planar,
        measure: 0.0,
    })
}

/// The bounds of a piece of an intrinsic curve (#361) over `[from, to]`
/// into it: the closed-form tangent law, with the read frames compared
/// against the exact ones. `None` where they are not bounded.
fn intrinsic_frame_bounds(
    reader: &IntrinsicReader<'_>,
    piece: &StationPiece<'_>,
    from: Scalar,
    to: Scalar,
    reading: FrameReading,
) -> Option<FrameLawBounds> {
    let (u0, u1) = (from.min(to), from.max(to));
    let length = piece.length();
    let slack = numerical_measure(length);
    if !(u0.is_finite() && u1.is_finite()) || u0 < -slack || u1 > length + slack {
        return None;
    }
    let (u0, u1) = (u0.max(0.0), u1.min(length));
    let (m0, m1) = if piece.reversed {
        (piece.end - u1, piece.end - u0)
    } else {
        (piece.start + u0, piece.start + u1)
    };
    let law = intrinsic_tangent_law(reader, m0, m1)?;
    let errors = reader.errors(piece, u0, u1).ok()??;
    finish(piece, law, m1 - m0, reading, errors)
}

/// A piece's reader, built when first asked for.
enum Slot<'c> {
    /// A piece whose frame is read in closed form.
    Closed,
    /// A piece read by quadrature, with its reader.
    Read(Box<PieceReader<'c>>),
    /// A piece of an intrinsic curve (#361), with its reader.
    Intrinsic(Box<IntrinsicReader<'c>>),
    /// A piece read by quadrature that is not bounded.
    Unbounded,
}

/// What a [`FrameLaw`] is the law along.
// Large only by `StationCurve`'s offset variant; one value lives per law.
#[allow(clippy::large_enum_variant)]
enum LawBasis<'b, 'c> {
    /// One atomic curve, possibly placed; `whole` the curve as one piece
    /// where it is read by quadrature.
    Atomic {
        curve: StationCurve<'c>,
        placement: Option<Transform3>,
        whole: Option<StationPiece<'c>>,
    },
    Composite(&'b CompositeBasis<'c>),
}

/// Certified bounds on the frame law along a station basis (#316): an
/// atomic curve or a [`CompositeBasis`], over spans that cross no seam.
///
/// A piece read in closed form is bounded by [`piece_frame_bounds`]. A
/// piece read by quadrature -- an ellipse, a B-spline, an offset beside a
/// line, a circle, an ellipse or a B-spline -- is bounded against its
/// exact frame law, its own length measured once by chord panels to the
/// tolerance (see [the module documentation](self#read-by-quadrature-316)):
/// the measure is how far a read's true distance lies from the one asked
/// for, and the read point's and axes' distances from the exact frame at
/// that true distance are [`FrameLawBounds::point_error`] and
/// [`FrameLawBounds::axis_error`]. On a composite, a piece after one
/// measured by quadrature starts where the kernel sums the lengths before
/// it, which [`CompositeBasis::joint_intervals`] certifies: its distance
/// from the true joint is added to the measure.
///
/// A 2D intrinsic curve, a level planar intrinsic space curve and an
/// elevated curve over an intrinsic plan are bounded by their curvature
/// law's closed-form suprema, each frame read compared with the exact one
/// (see [the module documentation](self#intrinsic-curves-361)).
///
/// Not bounded (`None`): every family [`piece_frame_bounds`] leaves
/// unbounded other than those above (a chain, an intrinsic space curve
/// with torsion or tilted, a banked curve over an intrinsic or chain plan,
/// an offset beside a curve parameterised by its plan length or beside a
/// banked curve, a skewed circle, a steep tangent), a curve read by
/// quadrature placed in a frame that tilts `+Z`, a placed offset piece,
/// and a piece after joints whose distances are not certified at the
/// tolerance.
pub struct FrameLaw<'b, 'c> {
    basis: LawBasis<'b, 'c>,
    tolerance: Scalar,
    /// Each joint's distance from where it truly lies, built when first
    /// needed; `None` where the joints are not certified.
    joints: OnceCell<Option<Vec<Scalar>>>,
    slots: Vec<OnceCell<Slot<'c>>>,
}

impl core::fmt::Debug for FrameLaw<'_, '_> {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        f.debug_struct("FrameLaw")
            .field("tolerance", &self.tolerance)
            .finish_non_exhaustive()
    }
}

/// A span of the basis, as the piece that holds it sees it.
struct Located<'b, 'c> {
    /// The piece's slot.
    slot: usize,
    piece: StationPiece<'c>,
    /// The span as distances into the piece.
    from: Scalar,
    to: Scalar,
    table: Option<&'b OffsetTable>,
    /// The composite the piece is in, whose joints add to the measure.
    composite: Option<&'b CompositeBasis<'c>>,
}

/// Whether a piece's frame is read by quadrature.
fn by_quadrature(curve: StationCurve<'_>) -> bool {
    matches!(curve, StationCurve::Offset(_)) || !curve.measure_is_exact()
}

/// Whether a piece's frame needs a reader: read by quadrature, or an
/// intrinsic curve whose point is compared with certified positions.
fn needs_reader(curve: StationCurve<'_>) -> bool {
    by_quadrature(curve) || reads_intrinsic(curve)
}

impl<'b, 'c> FrameLaw<'b, 'c> {
    /// The frame law along the atomic `curve`, placed by `placement` if
    /// any, its quadratures certified to `tolerance`.
    ///
    /// # Errors
    ///
    /// A tolerance that is not finite and positive, and a curve read by
    /// quadrature whose length a station cannot be measured against, by
    /// name.
    pub fn atomic(
        curve: StationCurve<'c>,
        placement: Option<Transform3>,
        tolerance: Scalar,
    ) -> GeomResult<Self> {
        let tolerance = checked_tolerance(tolerance)?;
        let whole = if needs_reader(curve) {
            let piece = StationPiece::whole(curve)?;
            Some(match placement {
                Some(rigid) => piece.placed(rigid, true),
                None => piece,
            })
        } else {
            None
        };
        Ok(Self {
            basis: LawBasis::Atomic {
                curve,
                placement,
                whole,
            },
            tolerance,
            joints: OnceCell::new(),
            slots: vec![OnceCell::new()],
        })
    }

    /// The frame law along `basis`, its quadratures certified to
    /// `tolerance`.
    ///
    /// # Errors
    ///
    /// A tolerance that is not finite and positive, by name.
    pub fn composite(basis: &'b CompositeBasis<'c>, tolerance: Scalar) -> GeomResult<Self> {
        let tolerance = checked_tolerance(tolerance)?;
        Ok(Self {
            basis: LawBasis::Composite(basis),
            tolerance,
            joints: OnceCell::new(),
            slots: basis.pieces().iter().map(|_| OnceCell::new()).collect(),
        })
    }

    /// Certified bounds of the frame law over the distances `[a, b]` along
    /// the basis, a span that crosses no seam, read as `reading`; `None`
    /// where they are not certified (see the [type documentation](Self)).
    #[must_use]
    pub fn bounds(&self, a: Scalar, b: Scalar, reading: FrameReading) -> Option<FrameLawBounds> {
        let at = self.locate(a, b)?;
        let mut bounds = self.read(at.slot, &at.piece, at.from, at.to, at.table, reading)?;
        if let Some(composite) = at.composite {
            bounds.measure += self.joint_error(composite, at.slot)?;
        }
        Some(bounds)
    }

    /// The piece holding the span `[a, b]` of the basis and the span as
    /// distances into it; `None` for a span that lies in no piece.
    fn locate(&self, a: Scalar, b: Scalar) -> Option<Located<'b, 'c>> {
        match &self.basis {
            LawBasis::Atomic {
                curve,
                placement,
                whole,
            } => match whole {
                Some(whole) => Some(Located {
                    slot: 0,
                    piece: *whole,
                    from: a,
                    to: b,
                    table: None,
                    composite: None,
                }),
                None => {
                    let mut piece = StationPiece::between(*curve, a, b).ok()?;
                    if let Some(rigid) = placement {
                        piece = piece.placed(*rigid, true);
                    }
                    Some(Located {
                        slot: 0,
                        piece,
                        from: 0.0,
                        to: b - a,
                        table: None,
                        composite: None,
                    })
                }
            },
            LawBasis::Composite(composite) => {
                let composite: &'b CompositeBasis<'c> = composite;
                let slack = numerical_measure(b);
                let starts = composite.starts();
                // The piece holding the span's middle: a span narrower than
                // the slack that starts at a joint is the later piece's.
                let mid = 0.5 * (a + b);
                composite
                    .pieces()
                    .iter()
                    .enumerate()
                    .find_map(|(k, piece)| {
                        let (start, end) = (starts[k], starts[k] + piece.length());
                        let holds = a >= start - slack && b <= end + slack;
                        (holds && mid >= start && mid <= end).then(|| Located {
                            slot: k,
                            piece: *piece,
                            from: a - start,
                            to: b - start,
                            table: composite.tables()[k].as_ref(),
                            composite: Some(composite),
                        })
                    })
            }
        }
    }

    /// Piece `k`'s reader, built when first asked for.
    fn slot(&self, k: usize, piece: &StationPiece<'c>, table: Option<&OffsetTable>) -> &Slot<'c> {
        self.slots[k].get_or_init(|| {
            if reads_intrinsic(piece.curve) {
                return match IntrinsicReader::new(piece, self.tolerance) {
                    Ok(Some(reader)) => Slot::Intrinsic(Box::new(reader)),
                    _ => Slot::Unbounded,
                };
            }
            if !by_quadrature(piece.curve) {
                return Slot::Closed;
            }
            match PieceReader::new(piece, table, self.tolerance) {
                Ok(Some(reader)) => Slot::Read(Box::new(reader)),
                _ => Slot::Unbounded,
            }
        })
    }

    /// The bounds of piece `k` over `[from, to]` into it.
    fn read(
        &self,
        k: usize,
        piece: &StationPiece<'c>,
        from: Scalar,
        to: Scalar,
        table: Option<&OffsetTable>,
        reading: FrameReading,
    ) -> Option<FrameLawBounds> {
        match self.slot(k, piece, table) {
            Slot::Closed => piece_frame_bounds(piece, from, to, reading),
            Slot::Read(reader) => read_frame_bounds(reader, piece, from, to, table, reading),
            Slot::Intrinsic(reader) => intrinsic_frame_bounds(reader, piece, from, to, reading),
            Slot::Unbounded => None,
        }
    }

    /// How far piece `k`'s start, where the kernel sums the lengths before
    /// it, may lie from where it truly lies: zero after pieces measured
    /// from stored data alone.
    fn joint_error(&self, composite: &CompositeBasis<'c>, k: usize) -> Option<Scalar> {
        if composite.pieces()[..k]
            .iter()
            .all(|piece| piece.curve.measure_is_exact())
        {
            return Some(0.0);
        }
        let errors = self.joints.get_or_init(|| {
            let joints = composite.joint_intervals(self.tolerance).ok()?;
            let starts = composite.starts();
            Some(
                joints
                    .iter()
                    .enumerate()
                    .map(|(j, joint)| {
                        let at = starts[j + 1];
                        (at - joint.lower).max(joint.upper - at).max(0.0) * (1.0 + ROUNDING)
                    })
                    .collect(),
            )
        });
        errors.as_ref()?.get(k.checked_sub(1)?).copied()
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::station::{station_section2, station_section3, SectionFrame};
    use axiolid_core::{Frame2, Point2, Point3};
    use axiolid_curve::{
        BankConvention, Banked3, CantLaw, CantPiece, Circle2, Circle3, Elevated3, ElevationLaw,
        Line2,
    };

    /// A frame's axis (or point) read as a vector.
    type Axis = fn(&SectionFrame) -> Vec3;

    /// Second differences of each axis and the point, over a lattice of
    /// the span, never above the bounds; first differences likewise.
    fn assert_bounded(
        section: impl Fn(Scalar) -> SectionFrame,
        bounds: &FrameLawBounds,
        lo: Scalar,
        hi: Scalar,
    ) {
        let h = 1e-4 * (hi - lo);
        for i in 1..20 {
            let s = lo + (hi - lo) * i as Scalar / 20.0;
            let (a, b, c) = (section(s - h), section(s), section(s + h));
            let first = |f: fn(&SectionFrame) -> Vec3| ((f(&c) - f(&a)) / (2.0 * h)).length();
            let second =
                |f: fn(&SectionFrame) -> Vec3| ((f(&c) - 2.0 * f(&b) + f(&a)) / (h * h)).length();
            let slack = |bound: Scalar| bound * (1.0 + 1e-3) + 1e-6;
            let axes: [(Axis, [Scalar; 2], &str); 4] = [
                (|f| f.point - Point3::ZERO, bounds.point, "point"),
                (|f| f.tangent, bounds.tangent, "tangent"),
                (|f| f.lateral, bounds.lateral, "lateral"),
                (|f| f.up, bounds.up, "up"),
            ];
            for (f, [b1, b2], name) in axes {
                assert!(first(f) <= slack(b1), "{name}' {} > {b1} at {s}", first(f));
                assert!(
                    second(f) <= slack(b2),
                    "{name}'' {} > {b2} at {s}",
                    second(f)
                );
            }
        }
    }

    fn whole3(curve: &Curve3, lo: Scalar, hi: Scalar) -> StationPiece<'_> {
        StationPiece::between(StationCurve::Three(curve), lo, hi).unwrap()
    }

    #[test]
    fn a_plan_circle_bounds_its_frame() {
        let curve = Curve2::Circle(Circle2 {
            frame: Frame2 {
                origin: Point2::new(1.0, 2.0),
                x: Vec2::X,
                y: Vec2::Y,
            },
            radius: 3.0,
        });
        let piece = StationPiece::between(StationCurve::Two(&curve), 1.0, 5.0).unwrap();
        let bounds = piece_frame_bounds(&piece, 0.0, 4.0, FrameReading::Section).unwrap();
        assert_bounded(
            |s| station_section2(&curve, 1.0 + s).unwrap(),
            &bounds,
            0.0,
            4.0,
        );
        assert!((bounds.point[1] - 1.0 / 3.0).abs() <= 1e-12);
    }

    /// A certified lower bound of `O'.N` in station distance, at sampled
    /// stations of a piece: `forward` never exceeds the forward rate the
    /// frames show.
    fn assert_forward(
        section: impl Fn(Scalar) -> SectionFrame,
        bounds: &FrameLawBounds,
        lo: Scalar,
        hi: Scalar,
    ) {
        let h = 1e-5 * (hi - lo);
        for i in 1..20 {
            let s = lo + (hi - lo) * i as Scalar / 20.0;
            let rate =
                (section(s + h).point - section(s - h).point).dot(section(s).tangent) / (2.0 * h);
            assert!(
                bounds.forward <= rate + 1e-6,
                "{} > {rate} at {s}",
                bounds.forward
            );
        }
    }

    #[test]
    fn a_plan_circle_advances_at_unit_rate() {
        let curve = Curve2::Circle(Circle2 {
            frame: Frame2 {
                origin: Point2::ZERO,
                x: Vec2::X,
                y: Vec2::Y,
            },
            radius: 10.0,
        });
        let piece = StationPiece::between(StationCurve::Two(&curve), 0.0, 20.0).unwrap();
        for reading in [FrameReading::Section, FrameReading::Plan] {
            let bounds = piece_frame_bounds(&piece, 0.0, 20.0, reading).unwrap();
            assert!(
                (1.0 - bounds.forward).abs() <= 2.0 * ROUNDING,
                "{}",
                bounds.forward
            );
            assert_forward(|s| station_section2(&curve, s).unwrap(), &bounds, 0.0, 20.0);
        }
    }

    #[test]
    fn a_graded_curve_advances_at_least_at_unit_rate_in_plan_distance() {
        let curve = Curve3::Elevated(Elevated3 {
            plan: Box::new(Curve2::Line(Line2 {
                origin: Point2::ZERO,
                direction: Vec2::X,
            })),
            elevation: ElevationLaw::Polynomial {
                coefficients: vec![0.0, 0.1],
            },
        });
        let piece = whole3(&curve, 0.0, 20.0);
        for reading in [FrameReading::Section, FrameReading::Plan] {
            let bounds = piece_frame_bounds(&piece, 0.0, 20.0, reading).unwrap();
            assert!(
                (1.0 - bounds.forward).abs() <= 2.0 * ROUNDING,
                "{}",
                bounds.forward
            );
        }
    }

    #[test]
    fn a_tilted_line_advances_at_its_horizontal_part_in_the_plan_reading() {
        let tilt: Scalar = 0.6;
        let curve = Curve3::Line(axiolid_curve::Line3 {
            origin: Point3::ZERO,
            direction: Vec3::new(tilt.cos(), 0.0, tilt.sin()),
        });
        let piece = whole3(&curve, 0.0, 10.0);
        let section = piece_frame_bounds(&piece, 0.0, 10.0, FrameReading::Section).unwrap();
        assert!((1.0 - section.forward).abs() <= 2.0 * ROUNDING);
        let plan = piece_frame_bounds(&piece, 0.0, 10.0, FrameReading::Plan).unwrap();
        assert!(plan.forward <= tilt.cos());
        assert!(
            tilt.cos() - plan.forward <= 4.0 * ROUNDING,
            "{}",
            plan.forward
        );
    }

    #[test]
    fn a_placement_scales_the_forward_rate_by_its_smallest_column() {
        let curve = Curve2::Line(Line2 {
            origin: Point2::ZERO,
            direction: Vec2::X,
        });
        let squeeze = Transform3::from_cols(
            Vec3::new(0.5, 0.0, 0.0),
            Vec3::new(0.0, 0.5, 0.0),
            Vec3::new(0.0, 0.0, 0.5),
            Point3::ZERO,
        );
        let piece = StationPiece::between(StationCurve::Two(&curve), 0.0, 10.0)
            .unwrap()
            .placed(squeeze, true);
        let bounds = piece_frame_bounds(&piece, 0.0, 10.0, FrameReading::Section).unwrap();
        assert!(bounds.forward <= 0.5);
        assert!(0.5 - bounds.forward <= 1e-12, "{}", bounds.forward);
    }

    #[test]
    fn a_tilted_circle_bounds_its_reference_up_frame() {
        let x = Vec3::new(1.0, 0.0, 0.0);
        let y = Vec3::new(0.0, 0.8, 0.6);
        let curve = Curve3::Circle(Circle3 {
            frame: axiolid_core::Frame3 {
                origin: Point3::ZERO,
                x,
                y,
                z: x.cross(y),
            },
            radius: 5.0,
        });
        // Away from where the tangent is steepest.
        let piece = whole3(&curve, 0.5, 2.5);
        for reading in [FrameReading::Section, FrameReading::Plan] {
            let bounds = piece_frame_bounds(&piece, 0.0, 2.0, reading).unwrap();
            assert_bounded(
                |s| {
                    let f = station_section3(&curve, 0.5 + s).unwrap();
                    match reading {
                        FrameReading::Section => f,
                        FrameReading::Plan => f.plan().unwrap(),
                    }
                },
                &bounds,
                0.0,
                2.0,
            );
        }
    }

    #[test]
    fn a_vertical_curve_bounds_its_leaning_frame() {
        let curve = Curve3::Elevated(Elevated3::new(
            Curve2::Circle(Circle2 {
                frame: Frame2 {
                    origin: Point2::ZERO,
                    x: Vec2::X,
                    y: Vec2::Y,
                },
                radius: 40.0,
            }),
            ElevationLaw::Polynomial {
                coefficients: vec![1.0, 0.03, -0.0004],
            },
        ));
        let piece = whole3(&curve, 0.0, 30.0);
        let bounds = piece_frame_bounds(&piece, 2.0, 28.0, FrameReading::Section).unwrap();
        assert_bounded(|s| station_section3(&curve, s).unwrap(), &bounds, 2.0, 28.0);
        assert_eq!(bounds.measure, 0.0);
    }

    #[test]
    fn a_banked_transition_bounds_its_rolled_frame() {
        // A long, gentle transition and a short, steep one, where the
        // roll's rate dominates the frame's turning.
        let cases = [
            (BankConvention::TangentRotation, 60.0, 10.0, 50.0),
            (BankConvention::VerticalRise, 60.0, 10.0, 50.0),
            (BankConvention::TangentRotation, 8.0, 1.0, 7.0),
            (BankConvention::VerticalRise, 8.0, 1.0, 7.0),
        ];
        for (convention, length, lo, hi) in cases {
            let curve = Curve3::Banked(Banked3::new(
                Elevated3::new(
                    Curve2::Circle(Circle2 {
                        frame: Frame2 {
                            origin: Point2::ZERO,
                            x: Vec2::X,
                            y: Vec2::Y,
                        },
                        radius: 300.0,
                    }),
                    ElevationLaw::constant_grade(60.0, 0.02),
                ),
                CantLaw::new(vec![CantPiece::cosine(length, 0.0, 0.15)]),
                CantLaw::zero(length),
                1.5,
                convention,
            ));
            let piece = whole3(&curve, 0.0, length);
            let bounds = piece_frame_bounds(&piece, lo, hi, FrameReading::Section).unwrap();
            assert_bounded(|s| station_section3(&curve, s).unwrap(), &bounds, lo, hi);
        }
    }

    #[test]
    fn lines_are_straight_and_numerical_families_unbounded() {
        let line = Curve2::Line(Line2 {
            origin: Point2::ZERO,
            direction: Vec2::new(3.0, 4.0),
        });
        let piece = StationPiece::between(StationCurve::Two(&line), 0.0, 10.0).unwrap();
        let bounds = piece_frame_bounds(&piece, 0.0, 10.0, FrameReading::Section).unwrap();
        assert!(bounds.is_straight());
        let ellipse = Curve3::Ellipse(axiolid_curve::Ellipse3 {
            frame: axiolid_core::Frame3 {
                origin: Point3::ZERO,
                x: Vec3::X,
                y: Vec3::Y,
                z: Vec3::Z,
            },
            semi_axis_x: 2.0,
            semi_axis_y: 1.0,
        });
        let piece = whole3(&ellipse, 0.0, 1.0);
        assert!(piece_frame_bounds(&piece, 0.0, 1.0, FrameReading::Section).is_none());
    }
}
