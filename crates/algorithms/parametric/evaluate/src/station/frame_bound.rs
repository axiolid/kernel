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
//! # Not bounded
//!
//! `None`, never a guess: an ellipse, a B-spline, an intrinsic curve, a
//! chain and an offset (their frames rest on a quadrature, an estimate); a
//! circle with a skewed frame; an elevated or banked curve with another
//! plan or an intrinsic profile; a span that crosses a seam, a vertex of a
//! polyline or a piece of a law; a tangent too steep for a horizontal
//! lower bound; a cant whose roll reaches a quarter turn.
//!
//! [`SectionFrame::carried`]: super::SectionFrame::carried

use axiolid_core::{Scalar, Vec2, Vec3};
use axiolid_curve::{Curve2, Curve3};

use super::{DistanceConvention, StationCurve, StationPiece, KEEPS_UP_TOLERANCE};
use crate::arc_parameter::ARC_LENGTH_TOLERANCE;
use crate::bound::ROUNDING;

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
/// an intrinsic curve, a chain, an offset, a skewed circle, an elevated or
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
    let mut law = match piece.curve {
        StationCurve::Two(curve) => tangent_law2(curve, m0, m1)?,
        StationCurve::Three(curve) => tangent_law3(curve, m0, m1)?,
        StationCurve::Offset(_) => return None,
    };
    if piece.reversed {
        law.middle = -law.middle;
    }
    let mut scale = 1.0;
    if let Some(rigid) = piece.placement {
        scale = placement_stretch(rigid)?;
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
        let rho = law.middle.x.hypot(law.middle.y) - 0.5 * (m1 - m0) * k1;
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
    let (tangent, lateral, up) = match reading {
        FrameReading::Section => (tangent, lateral, up),
        // A 2D curve's frame is already its plan frame.
        FrameReading::Plan if matches!(law.kind, Kind::Planar) => (tangent, lateral, up),
        FrameReading::Plan => {
            let flat = normalised(k1, k2, horizontal()?);
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
    };
    let finite = [out.point, out.tangent, out.lateral, out.up]
        .iter()
        .flatten()
        .chain([&out.measure])
        .all(|v| v.is_finite());
    finite.then_some(out)
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
