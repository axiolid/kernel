//! Certified derivative and chord bounds (#232).
//!
//! A chord of a curve, or a triangle of a surface, is only as close to the
//! exact geometry as the geometry's second derivative allows, and a sample
//! at the midpoint does not show that: a midpoint can sit on the chord while
//! the curve bulges either side of it. This module bounds the derivatives
//! from the data itself, so the numbers it returns are upper bounds, never
//! estimates.
//!
//! # The chord bound
//!
//! For a curve `c` that is `C^1` with a bounded second derivative on
//! `[a, b]`, the linear interpolant `l` of its ends satisfies
//! `|c(t) - l(t)| <= (b - a)^2 / 8 * sup |c''|`: componentwise,
//! `c - l = -int G(t, s) c''(s) ds` with the non-negative Green's function
//! of `-d^2/dt^2` on `[a, b]`, whose integral is `(t - a)(b - t) / 2 <=
//! (b - a)^2 / 8`. `l(t)` lies on the chord, so every point of the arc is
//! within that of the chord. A circular arc of at most half a turn uses its
//! exact sagitta instead.
//!
//! # Where the derivative bounds come from
//!
//! - Conics: `c'' = -(r_x cos t X + r_y sin t Y)`, bounded in closed form
//!   over the interval (the frame is used as given, so a skewed frame is
//!   accounted for, not assumed away).
//! - B-splines: the derivative of a B-spline is a B-spline on the derivative
//!   control polygon, and a B-spline lies in the convex hull of the control
//!   points that act on its span, so each derivative is bounded by the
//!   largest of those control points. A rational spline `c = A / w` is
//!   bounded through the quotient rule, with `w` bounded below by its
//!   smallest active weight (weights are positive) and `c` bounded by the
//!   control points' spread about their centre; the centre is subtracted
//!   first, which leaves every derivative unchanged.
//! - Elementary surfaces: closed forms over the parameter box. B-spline
//!   surfaces: the derivative control nets of the tensor product.
//! - Implicit curves (ADR 0077): cell by cell, the solved parameter's
//!   slope and bend through the implicit function theorem, with interval
//!   bounds of the field's partials over the stretch's box; a bridge cell
//!   by its cubic's control points (`implicit`, #249).
//!
//! A B-spline is only as smooth as its knots allow. A knot of multiplicity
//! `m` leaves a degree `d` spline `C^(d - m)` there, so a bound that needs a
//! continuous `k`-th derivative refuses an interval with such a knot
//! strictly inside: the derivative bound would hold on each side but the
//! Taylor argument does not cross the jump. [`continuity_breaks2`] and
//! [`continuity_breaks3`] name those knots so a caller can split there.
//!
//! Everything is computed in floating point from the data; the bounds are
//! inflated by [`ROUNDING`] relative so rounding in their own arithmetic
//! cannot make them undershoot.

use axiolid_core::{Frame2, Frame3, Scalar, Vec2, Vec3};
use axiolid_curve::{BSplineCurve, Curve2, Curve3};
use axiolid_surface::{BSplineSurface, Surface};

use crate::curve::span_in;
use crate::nurbs::SplineAxis;

mod implicit;

/// Relative inflation applied to every returned bound, so the rounding of
/// the bound's own arithmetic cannot make it undershoot.
pub const ROUNDING: Scalar = 1e-12;

fn inflate(value: Scalar) -> Scalar {
    value * (1.0 + ROUNDING)
}

/// Largest operator norm of a frame's linear part, bounded through its
/// Gram matrix: `sigma^2 = lambda_max(F^T F) <= 1 + |F^T F - I|_F`. Exact
/// for an orthonormal frame, and never below the true stretch for a skewed
/// or scaled one.
fn gram_stretch(columns: &[Vec3]) -> Scalar {
    let mut off: Scalar = 0.0;
    for (i, a) in columns.iter().enumerate() {
        for (j, b) in columns.iter().enumerate() {
            let target = if i == j { 1.0 } else { 0.0 };
            let entry = a.dot(*b) - target;
            off += entry * entry;
        }
    }
    (1.0 + off.sqrt()).sqrt()
}

/// The largest factor by which a 3D frame stretches a local vector.
#[must_use]
pub fn frame_stretch3(frame: &Frame3) -> Scalar {
    inflate(gram_stretch(&[frame.x, frame.y, frame.z]))
}

/// The largest factor by which a 2D frame stretches a local vector.
#[must_use]
pub fn frame_stretch2(frame: &Frame2) -> Scalar {
    let lift = |v: Vec2| Vec3::new(v.x, v.y, 0.0);
    inflate(gram_stretch(&[lift(frame.x), lift(frame.y)]))
}

// --- interval helpers -------------------------------------------------------

/// Whether `[lo, hi]` contains `offset + k * period` for some integer `k`.
fn contains_lattice(lo: Scalar, hi: Scalar, offset: Scalar, period: Scalar) -> bool {
    let k = ((lo - offset) / period).ceil();
    offset + k * period <= hi
}

/// `max cos^2` over `[lo, hi]`.
fn max_cos2(lo: Scalar, hi: Scalar) -> Scalar {
    if contains_lattice(lo, hi, 0.0, core::f64::consts::PI) {
        1.0
    } else {
        lo.cos().powi(2).max(hi.cos().powi(2))
    }
}

/// `min cos^2` over `[lo, hi]`.
fn min_cos2(lo: Scalar, hi: Scalar) -> Scalar {
    if contains_lattice(lo, hi, core::f64::consts::FRAC_PI_2, core::f64::consts::PI) {
        0.0
    } else {
        lo.cos().powi(2).min(hi.cos().powi(2))
    }
}

/// `max |cos|` over `[lo, hi]`.
fn max_abs_cos(lo: Scalar, hi: Scalar) -> Scalar {
    max_cos2(lo, hi).sqrt()
}

/// `max |sin|` over `[lo, hi]`.
fn max_abs_sin(lo: Scalar, hi: Scalar) -> Scalar {
    (1.0 - min_cos2(lo, hi)).max(0.0).sqrt()
}

/// `(min cos, max cos)` over `[lo, hi]`.
fn cos_range(lo: Scalar, hi: Scalar) -> (Scalar, Scalar) {
    let tau = core::f64::consts::TAU;
    let max = if contains_lattice(lo, hi, 0.0, tau) {
        1.0
    } else {
        lo.cos().max(hi.cos())
    };
    let min = if contains_lattice(lo, hi, core::f64::consts::PI, tau) {
        -1.0
    } else {
        lo.cos().min(hi.cos())
    };
    (min, max)
}

fn ordered(a: Scalar, b: Scalar) -> (Scalar, Scalar) {
    if a <= b {
        (a, b)
    } else {
        (b, a)
    }
}

/// `sup |c''|` of the conic `o + r_x cos t X + r_y sin t Y` over `[a, b]`.
///
/// `|c''|^2 = r_x^2 |X|^2 cos^2 + r_y^2 |Y|^2 sin^2 + 2 r_x r_y X.Y cos sin`;
/// the first two terms are `beta^2 + (alpha^2 - beta^2) cos^2`, extreme
/// where `cos^2` is, and the cross term is at most `r_x r_y |X.Y|`.
fn conic_second(rx: Scalar, ry: Scalar, x: Vec3, y: Vec3, a: Scalar, b: Scalar) -> Scalar {
    let (lo, hi) = ordered(a, b);
    let alpha2 = (rx * rx) * x.length_squared();
    let beta2 = (ry * ry) * y.length_squared();
    let cos2 = if alpha2 >= beta2 {
        max_cos2(lo, hi)
    } else {
        min_cos2(lo, hi)
    };
    let main = beta2 + (alpha2 - beta2) * cos2;
    let cross = (rx * ry * x.dot(y)).abs();
    (main.max(0.0) + cross).sqrt()
}

/// The exact sagitta bound of a circular arc of at most half a turn, when
/// the frame is orthonormal; `None` otherwise.
fn circle_sagitta(radius: Scalar, x: Vec3, y: Vec3, a: Scalar, b: Scalar) -> Option<Scalar> {
    let h = (b - a).abs();
    let orthonormal = (x.length_squared() - 1.0).abs() <= 1e-12
        && (y.length_squared() - 1.0).abs() <= 1e-12
        && x.dot(y).abs() <= 1e-12;
    (orthonormal && h <= core::f64::consts::PI).then(|| radius.abs() * (1.0 - (0.5 * h).cos()))
}

fn lift2(v: Vec2) -> Vec3 {
    Vec3::new(v.x, v.y, 0.0)
}

// --- B-spline derivative bounds ---------------------------------------------

/// Suprema of the first three derivatives of a curve over an interval.
///
/// Each is the supremum over the open pieces between the curve's own
/// breaks: at a knot where a derivative jumps, both one-sided values are
/// covered. Whether a Taylor argument may cross such a knot is the
/// caller's question; see [`continuity_breaks3`].
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct CurveDerivativeBounds {
    /// `sup |c'|`.
    pub first: Scalar,
    /// `sup |c''|`.
    pub second: Scalar,
    /// `sup |c'''|`.
    pub third: Scalar,
}

/// Derivative control polygon of one homogeneous axis: `d (P_{i+1} - P_i) /
/// (u_{i+d+1} - u_{i+1})`, the zero-width spans contributing nothing.
fn hodograph<const N: usize>(
    points: &[[Scalar; N]],
    weights: &[Scalar],
    knots: &[Scalar],
    degree: usize,
) -> (Vec<[Scalar; N]>, Vec<Scalar>) {
    let mut out_points = Vec::with_capacity(points.len().saturating_sub(1));
    let mut out_weights = Vec::with_capacity(weights.len().saturating_sub(1));
    for i in 0..points.len().saturating_sub(1) {
        let width = knots[i + degree + 1] - knots[i + 1];
        let factor = if width > 0.0 {
            degree as Scalar / width
        } else {
            0.0
        };
        out_points.push(core::array::from_fn(|k| {
            (points[i + 1][k] - points[i][k]) * factor
        }));
        out_weights.push((weights[i + 1] - weights[i]) * factor);
    }
    (out_points, out_weights)
}

fn norm<const N: usize>(p: &[Scalar; N]) -> Scalar {
    p.iter().map(|c| c * c).sum::<Scalar>().sqrt()
}

/// Rational derivative bounds from homogeneous ones: `A = w c`, so
/// `c' = (A' - w' c) / w`, `c'' = (A'' - 2 w' c' - w'' c) / w` and
/// `c''' = (A''' - 3 w'' c' - 3 w' c'' - w''' c) / w`.
fn quotient_bounds(a: [Scalar; 4], w: [Scalar; 4], radius: Scalar, w_min: Scalar) -> [Scalar; 3] {
    let first = (a[1] + w[1] * radius) / w_min;
    let second = (a[2] + 2.0 * w[1] * first + w[2] * radius) / w_min;
    let third = (a[3] + 3.0 * w[2] * first + 3.0 * w[1] * second + w[3] * radius) / w_min;
    [first, second, third]
}

/// `sup |c^(k)|`, `k = 1, 2, 3`, of a B-spline over `[a, b]`, and its
/// chord coefficient `K`: every point of the arc is within `(b - a)^2 / 8
/// K` of the chord. `None` for invalid data or a non-positive weight.
///
/// The chord coefficient is not `sup |c''|` for a rational spline. With
/// the centred homogeneous numerator `A = w c` and the linear
/// interpolants `l_A`, `l_w` of `A` and `w` over `[a, b]`, the point
/// `q = l_A / l_w` lies on the chord (a projective map keeps segments), and
/// `c - q = ((A - l_A) - q (w - l_w)) / w`. Both differences are
/// polynomial chord errors, so `|c - q| <= (b - a)^2 / 8 (|A''| + R
/// |w''|) / w_min` with `R` bounding `|q|` (the active points' spread):
/// no first-derivative terms, which the quotient rule for `c''` carries.
fn bspline_bounds<P, const N: usize>(
    curve: &BSplineCurve<P>,
    to: impl Fn(&P) -> [Scalar; N],
    a: Scalar,
    b: Scalar,
) -> Option<(CurveDerivativeBounds, Scalar)> {
    let count = curve.control_points.len();
    let axis = SplineAxis::new(
        &curve.knots,
        &curve.multiplicities,
        curve.degree,
        count,
        "B-spline curve",
    )
    .ok()?;
    let weights: Vec<Scalar> = match &curve.weights {
        Some(w) if w.len() == count => w.clone(),
        Some(_) => return None,
        None => vec![1.0; count],
    };
    if weights.iter().any(|w| !(w.is_finite() && *w > 0.0)) {
        return None;
    }
    let degree = axis.degree;
    let (lo, hi) = ordered(axis.clamp(a), axis.clamp(b));
    let first_span = span_in(&axis.knots, count, degree, lo);
    let last_span = span_in(&axis.knots, count, degree, hi);
    // Control points acting on the spans [first_span, last_span]; every
    // hodograph's active points fall in the same index window.
    let window = first_span - degree..=last_span;
    let points: Vec<[Scalar; N]> = curve.control_points.iter().map(&to).collect();
    if points.iter().flatten().any(|c| !c.is_finite()) {
        return None;
    }
    // Centre the active points so `|c|` (which the rational quotient rule
    // multiplies by weight derivatives) is their spread, not their offset.
    let mut low = [Scalar::INFINITY; N];
    let mut high = [Scalar::NEG_INFINITY; N];
    for p in &points[window.clone()] {
        for k in 0..N {
            low[k] = low[k].min(p[k]);
            high[k] = high[k].max(p[k]);
        }
    }
    let centre: [Scalar; N] = core::array::from_fn(|k| 0.5 * (low[k] + high[k]));
    let radius = points[window.clone()]
        .iter()
        .map(|p| norm(&core::array::from_fn::<Scalar, N, _>(|k| p[k] - centre[k])))
        .fold(0.0, Scalar::max);
    let w_min = weights[window.clone()]
        .iter()
        .copied()
        .fold(Scalar::INFINITY, Scalar::min);
    let mut hom: Vec<[Scalar; N]> = points
        .iter()
        .zip(&weights)
        .map(|(p, w)| core::array::from_fn(|k| (p[k] - centre[k]) * w))
        .collect();
    let mut hw = weights.clone();
    let mut knots: &[Scalar] = &axis.knots;
    let mut a_sup = [0.0; 4];
    let mut w_sup = [0.0; 4];
    for order in 1..=3 {
        let current_degree = degree + 1 - order;
        if current_degree == 0 || hom.len() < 2 {
            break;
        }
        let (next, next_w) = hodograph(&hom, &hw, knots, current_degree);
        knots = &knots[1..knots.len() - 1];
        let end = (*window.end()).min(next.len().saturating_sub(1));
        let start = (*window.start()).min(end);
        a_sup[order] = next[start..=end].iter().map(norm).fold(0.0, Scalar::max);
        w_sup[order] = next_w[start..=end]
            .iter()
            .map(|w| w.abs())
            .fold(0.0, Scalar::max);
        hom = next;
        hw = next_w;
    }
    let [first, second, third] = quotient_bounds(a_sup, w_sup, radius, w_min);
    let bounds = CurveDerivativeBounds {
        first: inflate(first),
        second: inflate(second),
        third: inflate(third),
    };
    let chord = inflate((a_sup[2] + radius * w_sup[2]) / w_min).min(bounds.second);
    (bounds.first.is_finite()
        && bounds.second.is_finite()
        && bounds.third.is_finite()
        && chord.is_finite())
    .then_some((bounds, chord))
}

/// Interior knots of a B-spline where it may fail to be `C^k`: those of
/// multiplicity above `degree - k`.
fn bspline_breaks<P>(curve: &BSplineCurve<P>, k: usize) -> Vec<Scalar> {
    let degree = usize::from(curve.degree);
    let interior = curve.knots.len().saturating_sub(2);
    curve
        .knots
        .iter()
        .zip(&curve.multiplicities)
        .skip(1)
        .take(interior)
        .filter(|(_, &m)| (m as usize) + k > degree)
        .map(|(&knot, _)| knot)
        .collect()
}

/// Parameters strictly inside the curve's domain where a 2D curve may fail
/// to be `C^k` (`k >= 1`): a polyline's vertices, a B-spline's knots of
/// multiplicity above `degree - k`, an arc-length chain's joins (where its
/// curvature may jump; the chain is `C^1` there, but each piece is bounded
/// on its own). Empty for the smooth families.
#[must_use]
pub fn continuity_breaks2(curve: &Curve2, k: usize) -> Vec<Scalar> {
    match curve {
        Curve2::BSpline(b) => bspline_breaks(b, k),
        Curve2::Polyline(p) if k >= 1 => polyline_breaks(p.points.len(), p.closed),
        // A chain is tangent-continuous at its joins by construction, but
        // each piece is bounded on its own, so the joins are named for
        // every `k`.
        Curve2::Chain(c) if k >= 1 => c.joins().unwrap_or_default(),
        // An implicit curve's cells are bounded one at a time (#249).
        Curve2::Implicit(c) if k >= 1 => implicit::joins(c),
        _ => Vec::new(),
    }
}

/// [`continuity_breaks2`] for a 3D curve.
#[must_use]
pub fn continuity_breaks3(curve: &Curve3, k: usize) -> Vec<Scalar> {
    match curve {
        Curve3::BSpline(b) => bspline_breaks(b, k),
        Curve3::Polyline(p) if k >= 1 => polyline_breaks(p.points.len(), p.closed),
        // A chain plan's joins and the profile's seams; curvature seams of
        // the plan for `k >= 2` (#252).
        Curve3::Elevated(e) => crate::elevated::elevated_breaks(e, k),
        Curve3::Banked(b) => crate::elevated::banked_breaks(b, k),
        _ => Vec::new(),
    }
}

fn polyline_breaks(count: usize, closed: bool) -> Vec<Scalar> {
    let segments = if closed {
        count
    } else {
        count.saturating_sub(1)
    };
    (1..segments).map(|k| k as Scalar).collect()
}

fn breaks_inside(breaks: &[Scalar], a: Scalar, b: Scalar) -> bool {
    let (lo, hi) = ordered(a, b);
    breaks.iter().any(|&t| t > lo && t < hi)
}

/// Certified derivative suprema of a 3D curve over `[a, b]`, or `None` for
/// a family this module cannot bound.
///
/// Covered: lines, circles, ellipses and B-splines (rational ones too).
#[must_use]
pub fn curve_derivative_bounds3(
    curve: &Curve3,
    a: Scalar,
    b: Scalar,
) -> Option<CurveDerivativeBounds> {
    if !(a.is_finite() && b.is_finite()) {
        return None;
    }
    match curve {
        Curve3::Line(line) => Some(CurveDerivativeBounds {
            first: inflate(line.direction.length()),
            second: 0.0,
            third: 0.0,
        }),
        Curve3::Circle(c) => Some(conic_bounds(c.radius, c.radius, c.frame.x, c.frame.y, a, b)),
        Curve3::Ellipse(e) => Some(conic_bounds(
            e.semi_axis_x,
            e.semi_axis_y,
            e.frame.x,
            e.frame.y,
            a,
            b,
        )),
        Curve3::BSpline(spline) => bspline_bounds(spline, |p| [p.x, p.y, p.z], a, b).map(|b| b.0),
        // In plan distance, plan and profile combined (#252).
        Curve3::Elevated(e) => crate::elevated::elevated_derivative_bounds(e, a, b),
        Curve3::Banked(banked) => crate::elevated::banked_derivative_bounds(banked, a, b),
        _ => None,
    }
}

/// Certified derivative suprema of a 2D curve over `[a, b]`; see
/// [`curve_derivative_bounds3`]. Also covers sinusoids.
#[must_use]
pub fn curve_derivative_bounds2(
    curve: &Curve2,
    a: Scalar,
    b: Scalar,
) -> Option<CurveDerivativeBounds> {
    if !(a.is_finite() && b.is_finite()) {
        return None;
    }
    match curve {
        Curve2::Line(line) => Some(CurveDerivativeBounds {
            first: inflate(line.direction.length()),
            second: 0.0,
            third: 0.0,
        }),
        Curve2::Circle(c) => Some(conic_bounds(
            c.radius,
            c.radius,
            lift2(c.frame.x),
            lift2(c.frame.y),
            a,
            b,
        )),
        Curve2::Ellipse(e) => Some(conic_bounds(
            e.semi_axis_x,
            e.semi_axis_y,
            lift2(e.frame.x),
            lift2(e.frame.y),
            a,
            b,
        )),
        Curve2::BSpline(spline) => bspline_bounds(spline, |p| [p.x, p.y], a, b).map(|b| b.0),
        // `(t, mean + a cos t + b sin t)`: every derivative's wave has
        // amplitude `sqrt(a^2 + b^2)`.
        Curve2::Sinusoid(w) => {
            let amplitude = w.cosine.hypot(w.sine);
            Some(CurveDerivativeBounds {
                first: inflate((1.0 + amplitude * amplitude).sqrt()),
                second: inflate(amplitude),
                third: inflate(amplitude),
            })
        }
        _ => None,
    }
}

/// A conic's derivatives cycle: `c' = -c'''` is `c''` a quarter turn on,
/// so `|c'|^2 <= max(r_x^2 |X|^2, r_y^2 |Y|^2) + r_x r_y |X.Y|` bounds
/// both; `c''` gets the interval's own tighter bound.
fn conic_bounds(
    rx: Scalar,
    ry: Scalar,
    x: Vec3,
    y: Vec3,
    a: Scalar,
    b: Scalar,
) -> CurveDerivativeBounds {
    let speed = ((rx * rx * x.length_squared()).max(ry * ry * y.length_squared())
        + (rx * ry * x.dot(y)).abs())
    .sqrt();
    CurveDerivativeBounds {
        first: inflate(speed),
        second: inflate(conic_second(rx, ry, x, y, a, b)),
        third: inflate(speed),
    }
}

/// Certified chord bound of a 2D curve over `[a, b]`.
///
/// An upper bound on the distance from any point of the curve between `a`
/// and `b` to the chord joining `c(a)` and `c(b)`. `None` when the family
/// has no derivative bound here, or when the curve may not be `C^1`
/// strictly inside the interval (a polyline vertex, a B-spline knot of full
/// multiplicity): split there first ([`continuity_breaks2`]).
#[must_use]
pub fn chord_bound2(curve: &Curve2, a: Scalar, b: Scalar) -> Option<Scalar> {
    if !(a.is_finite() && b.is_finite()) {
        return None;
    }
    if a == b {
        return Some(0.0);
    }
    if breaks_inside(&continuity_breaks2(curve, 1), a, b) {
        return None;
    }
    match curve {
        Curve2::Line(_) | Curve2::Polyline(_) => Some(0.0),
        Curve2::Circle(c) => Some(inflate(
            circle_sagitta(c.radius, lift2(c.frame.x), lift2(c.frame.y), a, b).unwrap_or_else(
                || {
                    taylor(
                        b - a,
                        conic_second(c.radius, c.radius, lift2(c.frame.x), lift2(c.frame.y), a, b),
                    )
                },
            ),
        )),
        Curve2::BSpline(spline) => bspline_bounds(spline, |p| [p.x, p.y], a, b)
            .map(|(_, chord)| inflate(taylor(b - a, chord))),
        Curve2::Chain(chain) => crate::chain::chain_chord_bound(chain, a, b),
        Curve2::Implicit(c) => {
            let (lo, hi) = ordered(a, b);
            implicit::chord_bound(c, lo, hi).map(inflate)
        }
        _ => curve_derivative_bounds2(curve, a, b).map(|d| inflate(taylor(b - a, d.second))),
    }
}

/// Certified chord bound of a 3D curve over `[a, b]`; see [`chord_bound2`].
#[must_use]
pub fn chord_bound3(curve: &Curve3, a: Scalar, b: Scalar) -> Option<Scalar> {
    if !(a.is_finite() && b.is_finite()) {
        return None;
    }
    if a == b {
        return Some(0.0);
    }
    if breaks_inside(&continuity_breaks3(curve, 1), a, b) {
        return None;
    }
    match curve {
        Curve3::Line(_) | Curve3::Polyline(_) => Some(0.0),
        Curve3::Circle(c) => Some(inflate(
            circle_sagitta(c.radius, c.frame.x, c.frame.y, a, b).unwrap_or_else(|| {
                taylor(
                    b - a,
                    conic_second(c.radius, c.radius, c.frame.x, c.frame.y, a, b),
                )
            }),
        )),
        Curve3::BSpline(spline) => bspline_bounds(spline, |p| [p.x, p.y, p.z], a, b)
            .map(|(_, chord)| inflate(taylor(b - a, chord))),
        // Plan and profile chord bounds at equal plan distance, composed
        // in quadrature (#252; `crate::elevated`).
        Curve3::Elevated(e) => crate::elevated::elevated_chord_bound(e, a, b),
        Curve3::Banked(banked) => crate::elevated::banked_chord_bound(banked, a, b),
        _ => curve_derivative_bounds3(curve, a, b).map(|d| inflate(taylor(b - a, d.second))),
    }
}

/// `h^2 / 8 * sup |c''|`.
fn taylor(h: Scalar, second: Scalar) -> Scalar {
    h * h * 0.125 * second
}

/// Whether [`crate::curve::flatten2`] certifies its chords for this family:
/// every point of the curve then lies within the chord tolerance of the
/// returned polyline (#232). Lines, polylines, circles, ellipses,
/// sinusoids, B-splines with positive weights and arc-length chains whose
/// pieces are all bounded ([`crate::chain::chain_certifies`]); other families are
/// flattened on their midpoint sagitta alone, which is a measurement, not a
/// bound.
#[must_use]
pub fn certifies_flattening2(curve: &Curve2) -> bool {
    match curve {
        Curve2::Line(_)
        | Curve2::Polyline(_)
        | Curve2::Circle(_)
        | Curve2::Ellipse(_)
        | Curve2::Sinusoid(_) => true,
        Curve2::BSpline(b) => spline_weights_positive(b.weights.as_deref()),
        Curve2::Chain(c) => crate::chain::chain_certifies(c),
        _ => false,
    }
}

/// [`certifies_flattening2`] for [`crate::curve::flatten3`].
#[must_use]
pub fn certifies_flattening3(curve: &Curve3) -> bool {
    match curve {
        Curve3::Line(_) | Curve3::Polyline(_) | Curve3::Circle(_) | Curve3::Ellipse(_) => true,
        Curve3::BSpline(b) => spline_weights_positive(b.weights.as_deref()),
        // Certified where both halves always bound: a certified plan and
        // a profile of closed-form pieces (#252). An intrinsic profile's
        // bound rests on a below-vertical certificate that may run out of
        // budget, and a banked curve's on its pivot forms: neither is
        // claimed for every span.
        Curve3::Elevated(e) => crate::elevated::certifies_elevated(e),
        Curve3::Banked(_) => false,
        _ => false,
    }
}

fn spline_weights_positive(weights: Option<&[Scalar]>) -> bool {
    weights.is_none_or(|w| w.iter().all(|w| w.is_finite() && *w > 0.0))
}

// --- surfaces ---------------------------------------------------------------

/// Suprema of a surface's partial derivatives over a parameter box.
///
/// `du`, `dv` bound `|S_u|`, `|S_v|`; `duu`, `duv`, `dvv` bound the second
/// partials. `smooth` is false when a B-spline knot line of full
/// multiplicity crosses the box: the second-derivative bounds then hold on
/// each side but not across, and a caller must fall back to a first-order
/// argument there.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct SurfaceBounds {
    /// `sup |S_u|`.
    pub du: Scalar,
    /// `sup |S_v|`.
    pub dv: Scalar,
    /// `sup |S_uu|`.
    pub duu: Scalar,
    /// `sup |S_uv|`.
    pub duv: Scalar,
    /// `sup |S_vv|`.
    pub dvv: Scalar,
    /// Whether the second-derivative bounds may be used across the box.
    pub smooth: bool,
    /// Coefficients `[A, B, C]` for linear interpolation over a parameter
    /// triangle in the box: every surface point over the triangle is within
    /// `1/2 sum_i l_i q(y_i - x)`, `q(d) = A d_u^2 + 2 B |d_u d_v| + C
    /// d_v^2`, of the triangle spanned by the corners' surface points.
    ///
    /// For a polynomial surface they are `duu, duv, dvv`. For a rational
    /// B-spline `S = A / w` they are `(|A_uu| + R |w_uu|) / w_min` and so on
    /// (centred numerator, `R` the active points' spread): the point
    /// `l_A / l_w` of the homogeneous interpolants is a convex combination
    /// of the corners' surface points, so it lies in their triangle, and
    /// `S - l_A / l_w = ((A - l_A) - (l_A / l_w)(w - l_w)) / w` carries no
    /// first-derivative terms, unlike the quotient rule for `S_uu`.
    pub interpolation: [Scalar; 3],
}

impl SurfaceBounds {
    fn inflated(self) -> Self {
        Self {
            du: inflate(self.du),
            dv: inflate(self.dv),
            duu: inflate(self.duu),
            duv: inflate(self.duv),
            dvv: inflate(self.dvv),
            smooth: self.smooth,
            interpolation: self.interpolation.map(inflate),
        }
    }

    fn finite(self) -> Option<Self> {
        [self.du, self.dv, self.duu, self.duv, self.dvv]
            .iter()
            .chain(&self.interpolation)
            .all(|v| v.is_finite())
            .then_some(self)
    }
}

/// Derivative bounds of one surface, prepared once for many boxes.
///
/// A B-spline surface's derivative control nets are built on construction,
/// so each [`Self::bounds`] call only scans the window of the box.
#[derive(Debug, Clone)]
pub struct SurfaceBoundOracle<'a> {
    surface: &'a Surface,
    spline: Option<SplineNets>,
}

/// Homogeneous derivative nets of a B-spline surface, centred.
#[derive(Debug, Clone)]
struct SplineNets {
    u_knots: Vec<Scalar>,
    v_knots: Vec<Scalar>,
    p: usize,
    q: usize,
    rows: usize,
    cols: usize,
    points: Vec<Vec<[Scalar; 3]>>,
    weights: Vec<Vec<Scalar>>,
    u_breaks: Vec<Scalar>,
    v_breaks: Vec<Scalar>,
}

impl<'a> SurfaceBoundOracle<'a> {
    /// Prepare the bounds of `surface`, or `None` for a family (or data)
    /// this module cannot bound.
    #[must_use]
    pub fn new(surface: &'a Surface) -> Option<Self> {
        let spline = match surface {
            Surface::Plane(_)
            | Surface::Cylinder(_)
            | Surface::EllipticalCylinder(_)
            | Surface::Cone(_)
            | Surface::Sphere(_)
            | Surface::Torus(_) => None,
            Surface::BSpline(b) => Some(SplineNets::new(b)?),
            _ => return None,
        };
        Some(Self { surface, spline })
    }

    /// Bounds over the box `[u.0, u.1] x [v.0, v.1]`; `None` where the
    /// surface is not defined (a cone past its apex) or the data is not
    /// finite.
    #[must_use]
    pub fn bounds(&self, u: (Scalar, Scalar), v: (Scalar, Scalar)) -> Option<SurfaceBounds> {
        let (u0, u1) = ordered(u.0, u.1);
        let (v0, v1) = ordered(v.0, v.1);
        if ![u0, u1, v0, v1].iter().all(|x| x.is_finite()) {
            return None;
        }
        let raw = match self.surface {
            Surface::Plane(p) => SurfaceBounds {
                du: p.frame.x.length(),
                dv: p.frame.y.length(),
                duu: 0.0,
                duv: 0.0,
                dvv: 0.0,
                smooth: true,
                interpolation: [0.0; 3],
            },
            Surface::Cylinder(c) => {
                let s = frame_stretch3(&c.frame);
                SurfaceBounds {
                    du: s * c.radius.abs(),
                    dv: c.frame.z.length(),
                    duu: s * c.radius.abs(),
                    duv: 0.0,
                    dvv: 0.0,
                    smooth: true,
                    interpolation: [0.0; 3],
                }
            }
            Surface::EllipticalCylinder(c) => {
                let s = frame_stretch3(&c.frame);
                let r = c.semi_axis_x.abs().max(c.semi_axis_y.abs());
                SurfaceBounds {
                    du: s * r,
                    dv: c.frame.z.length(),
                    duu: s * r,
                    duv: 0.0,
                    dvv: 0.0,
                    smooth: true,
                    interpolation: [0.0; 3],
                }
            }
            Surface::Cone(c) => {
                let s = frame_stretch3(&c.frame);
                let slope = c.semi_angle.tan();
                let (r0, r1) = (c.radius + v0 * slope, c.radius + v1 * slope);
                if r0 < 0.0 || r1 < 0.0 {
                    return None;
                }
                let rho = r0.max(r1);
                SurfaceBounds {
                    du: s * rho,
                    dv: s * (1.0 + slope * slope).sqrt(),
                    duu: s * rho,
                    duv: s * slope.abs(),
                    dvv: 0.0,
                    smooth: true,
                    interpolation: [0.0; 3],
                }
            }
            Surface::Sphere(sphere) => {
                let s = frame_stretch3(&sphere.frame) * sphere.radius.abs();
                let cos = max_abs_cos(v0, v1);
                SurfaceBounds {
                    du: s * cos,
                    dv: s,
                    duu: s * cos,
                    duv: s * max_abs_sin(v0, v1),
                    dvv: s,
                    smooth: true,
                    interpolation: [0.0; 3],
                }
            }
            Surface::Torus(t) => {
                let s = frame_stretch3(&t.frame);
                let (cmin, cmax) = cos_range(v0, v1);
                let r = t.minor_radius.abs();
                let ring = (t.major_radius + r * cmin)
                    .abs()
                    .max((t.major_radius + r * cmax).abs());
                SurfaceBounds {
                    du: s * ring,
                    dv: s * r,
                    duu: s * ring,
                    duv: s * r * max_abs_sin(v0, v1),
                    dvv: s * r,
                    smooth: true,
                    interpolation: [0.0; 3],
                }
            }
            Surface::BSpline(_) => self.spline.as_ref()?.bounds(u0, u1, v0, v1)?,
            _ => return None,
        };
        let raw = if matches!(self.surface, Surface::BSpline(_)) {
            raw
        } else {
            SurfaceBounds {
                interpolation: [raw.duu, raw.duv, raw.dvv],
                ..raw
            }
        };
        raw.inflated().finite()
    }
}

impl SplineNets {
    fn new(b: &BSplineSurface) -> Option<Self> {
        let rows = b.control_points.len();
        let cols = b.control_points.first()?.len();
        if b.control_points.iter().any(|row| row.len() != cols) {
            return None;
        }
        let ua = SplineAxis::new(&b.u_knots, &b.u_multiplicities, b.u_degree, rows, "u").ok()?;
        let va = SplineAxis::new(&b.v_knots, &b.v_multiplicities, b.v_degree, cols, "v").ok()?;
        let weights: Vec<Vec<Scalar>> = match &b.weights {
            Some(net) => {
                if net.len() != rows || net.iter().any(|row| row.len() != cols) {
                    return None;
                }
                net.clone()
            }
            None => vec![vec![1.0; cols]; rows],
        };
        if weights
            .iter()
            .flatten()
            .any(|w| !(w.is_finite() && *w > 0.0))
        {
            return None;
        }
        let points: Vec<Vec<[Scalar; 3]>> = b
            .control_points
            .iter()
            .map(|row| row.iter().map(|p| [p.x, p.y, p.z]).collect())
            .collect();
        if points.iter().flatten().flatten().any(|c| !c.is_finite()) {
            return None;
        }
        let breaks = |knots: &[Scalar], mults: &[u32], degree: u16| {
            let interior = knots.len().saturating_sub(2);
            knots
                .iter()
                .zip(mults)
                .skip(1)
                .take(interior)
                .filter(|(_, &m)| m as usize >= usize::from(degree))
                .map(|(&k, _)| k)
                .collect::<Vec<_>>()
        };
        Some(Self {
            u_breaks: breaks(&b.u_knots, &b.u_multiplicities, b.u_degree),
            v_breaks: breaks(&b.v_knots, &b.v_multiplicities, b.v_degree),
            u_knots: ua.knots,
            v_knots: va.knots,
            p: ua.degree,
            q: va.degree,
            rows,
            cols,
            points,
            weights,
        })
    }

    /// Bounds over a box: centre the active net, take its derivative nets
    /// in the window and apply the quotient rule.
    #[allow(clippy::needless_range_loop)]
    fn bounds(&self, u0: Scalar, u1: Scalar, v0: Scalar, v1: Scalar) -> Option<SurfaceBounds> {
        let clamp = |x: Scalar, knots: &[Scalar], d: usize, n: usize| x.clamp(knots[d], knots[n]);
        let (cu0, cu1) = (
            clamp(u0, &self.u_knots, self.p, self.rows),
            clamp(u1, &self.u_knots, self.p, self.rows),
        );
        let (cv0, cv1) = (
            clamp(v0, &self.v_knots, self.q, self.cols),
            clamp(v1, &self.v_knots, self.q, self.cols),
        );
        let su0 = span_in(&self.u_knots, self.rows, self.p, cu0);
        let su1 = span_in(&self.u_knots, self.rows, self.p, cu1);
        let sv0 = span_in(&self.v_knots, self.cols, self.q, cv0);
        let sv1 = span_in(&self.v_knots, self.cols, self.q, cv1);
        let (r0, r1) = (su0 - self.p, su1);
        let (c0, c1) = (sv0 - self.q, sv1);

        // Centre the active points.
        let mut low = [Scalar::INFINITY; 3];
        let mut high = [Scalar::NEG_INFINITY; 3];
        for row in &self.points[r0..=r1] {
            for p in &row[c0..=c1] {
                for k in 0..3 {
                    low[k] = low[k].min(p[k]);
                    high[k] = high[k].max(p[k]);
                }
            }
        }
        let centre: [Scalar; 3] = core::array::from_fn(|k| 0.5 * (low[k] + high[k]));
        let mut radius: Scalar = 0.0;
        let mut w_min = Scalar::INFINITY;
        // The window plus the rows and columns a second hodograph reads.
        let (wr0, wr1) = (r0, r1.min(self.rows - 1));
        let (wc0, wc1) = (c0, c1.min(self.cols - 1));
        let mut hom: Vec<Vec<[Scalar; 4]>> = Vec::with_capacity(wr1 + 1 - wr0);
        for i in wr0..=wr1 {
            let mut row = Vec::with_capacity(wc1 + 1 - wc0);
            for j in wc0..=wc1 {
                let p = self.points[i][j];
                let w = self.weights[i][j];
                let d: [Scalar; 3] = core::array::from_fn(|k| p[k] - centre[k]);
                radius = radius.max(norm(&d));
                w_min = w_min.min(w);
                row.push([d[0] * w, d[1] * w, d[2] * w, w]);
            }
            hom.push(row);
        }
        // Derivatives along u act on rows, along v on columns; the knot
        // slices are offset by the window start.
        let du_net =
            |net: &Vec<Vec<[Scalar; 4]>>, knots: &[Scalar], degree: usize, offset: usize| {
                let mut out = Vec::new();
                for i in 0..net.len().saturating_sub(1) {
                    let width = knots[offset + i + degree + 1] - knots[offset + i + 1];
                    let f = if width > 0.0 {
                        degree as Scalar / width
                    } else {
                        0.0
                    };
                    out.push(
                        net[i]
                            .iter()
                            .zip(&net[i + 1])
                            .map(|(a, b)| core::array::from_fn(|k| (b[k] - a[k]) * f))
                            .collect::<Vec<[Scalar; 4]>>(),
                    );
                }
                out
            };
        let dv_net =
            |net: &Vec<Vec<[Scalar; 4]>>, knots: &[Scalar], degree: usize, offset: usize| {
                net.iter()
                    .map(|row| {
                        (0..row.len().saturating_sub(1))
                            .map(|j| {
                                let width = knots[offset + j + degree + 1] - knots[offset + j + 1];
                                let f = if width > 0.0 {
                                    degree as Scalar / width
                                } else {
                                    0.0
                                };
                                core::array::from_fn(|k| (row[j + 1][k] - row[j][k]) * f)
                            })
                            .collect::<Vec<[Scalar; 4]>>()
                    })
                    .collect::<Vec<_>>()
            };
        let sup = |net: &Vec<Vec<[Scalar; 4]>>| {
            let mut a: Scalar = 0.0;
            let mut w: Scalar = 0.0;
            for p in net.iter().flatten() {
                a = a.max((p[0] * p[0] + p[1] * p[1] + p[2] * p[2]).sqrt());
                w = w.max(p[3].abs());
            }
            (a, w)
        };
        // The k-th u-hodograph uses the knots of the (k-1)-th trimmed at
        // both ends; in full-vector indexing its knot offset grows by one.
        let (p, q) = (self.p, self.q);
        let nu = if p >= 1 {
            du_net(&hom, &self.u_knots, p, wr0)
        } else {
            Vec::new()
        };
        let nuu = if p >= 2 {
            du_net(&nu, &self.u_knots, p - 1, wr0 + 1)
        } else {
            Vec::new()
        };
        let nv = if q >= 1 {
            dv_net(&hom, &self.v_knots, q, wc0)
        } else {
            Vec::new()
        };
        let nvv = if q >= 2 {
            dv_net(&nv, &self.v_knots, q - 1, wc0 + 1)
        } else {
            Vec::new()
        };
        let nuv = if p >= 1 && q >= 1 {
            dv_net(&nu, &self.v_knots, q, wc0)
        } else {
            Vec::new()
        };
        let (au, wu) = sup(&nu);
        let (av, wv) = sup(&nv);
        let (auu, wuu) = sup(&nuu);
        let (avv, wvv) = sup(&nvv);
        let (auv, wuv) = sup(&nuv);
        if !(w_min.is_finite() && w_min > 0.0) {
            return None;
        }
        // S = A / w: S_u = (A_u - w_u S) / w, S_uu = (A_uu - 2 w_u S_u -
        // w_uu S) / w, S_uv = (A_uv - w_u S_v - w_v S_u - w_uv S) / w.
        let su = (au + wu * radius) / w_min;
        let sv = (av + wv * radius) / w_min;
        let suu = (auu + 2.0 * wu * su + wuu * radius) / w_min;
        let svv = (avv + 2.0 * wv * sv + wvv * radius) / w_min;
        let suv = (auv + wu * sv + wv * su + wuv * radius) / w_min;
        // Evaluation clamps to the domain, which kinks the surface at its
        // edge, so a box reaching past it is not smooth either.
        let inside = cu0 == u0 && cu1 == u1 && cv0 == v0 && cv1 == v1;
        let smooth = inside
            && !breaks_inside(&self.u_breaks, cu0, cu1)
            && !breaks_inside(&self.v_breaks, cv0, cv1);
        Some(SurfaceBounds {
            du: su,
            dv: sv,
            duu: suu,
            duv: suv,
            dvv: svv,
            smooth,
            interpolation: [
                (auu + radius * wuu) / w_min,
                (auv + radius * wuv) / w_min,
                (avv + radius * wvv) / w_min,
            ],
        })
    }
}

#[cfg(test)]
mod tests;
