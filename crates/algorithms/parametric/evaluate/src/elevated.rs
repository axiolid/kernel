//! Derivatives, certified bounds and seams of elevated and banked curves in
//! plan distance (#252).
//!
//! An [`Elevated3`] is `c(d) = (p(d), z(d))`: a plan `p` read by its own
//! arc length `d` and a height `z` written against that same `d`. Its
//! parameter is plan distance (ADR 0082), so its derivatives are taken in
//! `d`: `c' = (p', z')` with `|p'| = 1` (not a unit vector in 3D: its length
//! is `sqrt(1 + g^2)`, `g` the grade), `c'' = (k n, z'')` with `k` the plan
//! curvature and `n` the plan's left normal, `c''' = (k' n - k^2 p', z''')`.
//! A [`Banked3`]'s point is its centreline raised by the pivot `e(d)`, so it
//! adds `e^(j)` to every `z^(j)`.
//!
//! # The certified 3D chord bound
//!
//! Over a span `[a, b]` of width `h`, let `l` be the linear interpolant of
//! the ends, `l(d) = c(a) + (d - a) / h (c(b) - c(a))`. It lies on the 3D
//! chord. Componentwise, `c - l = -int G(d, s) c''(s) ds` with the
//! non-negative Green's function of `-d^2/dd^2` on `[a, b]`, whose integral
//! is at most `h^2 / 8` ([`crate::bound`]'s chord argument). Applied to the
//! horizontal and the vertical part separately, at the SAME `d`:
//!
//! - `|p(d) - l_p(d)| <= P = h^2 / 8 sup |p''| = h^2 / 8 sup |k|`, the plan's
//!   own chord bound read at equal plan distance;
//! - `|z(d) - l_z(d)| <= Z`, the profile's chord bound
//!   ([`elevation_chord_bound`](crate::elevation::elevation_chord_bound)),
//!   which is the same Taylor argument on `z`;
//!
//! and the two parts are orthogonal (`p` horizontal, `z` vertical), so
//!
//! ```text
//! |c(d) - l(d)| <= sqrt(P^2 + Z^2).
//! ```
//!
//! Every point of the curve is within that of the chord. The bound only
//! needs `c` to be `C^1` with a bounded second derivative on the span: a
//! curvature seam of the plan (a clothoid meeting an arc) does not stop it,
//! a join of a chain or a seam of the profile (where the grade may jump)
//! does, and those are cut first ([`elevated_breaks`]). A flattener that
//! accepts a span once `sqrt(P^2 + Z^2) <= tol` has split the budget
//! between plan and profile in quadrature: `P <= tol cos(a)`,
//! `Z <= tol sin(a)` for the angle `a = atan(Z / P)` the span itself sets.
//! A banked curve's vertical part is `z + e`, bounded by `Z + h^2 / 8 sup
//! |e''|`.
//!
//! The plan's chord bound is NOT the plan's own `chord_bound2`: that bounds
//! the distance from the plan arc to the plan chord, which for a circle is
//! the exact sagitta at the arc's own nearest chord point, not at equal
//! plan distance, and the vertical part could then be read at a different
//! `d`. The Taylor bound at equal `d` is what composes.
//!
//! # Derivative bounds
//!
//! For a sweep's certificate (the tube wall of `axiolid-mesh-compile`'s
//! deviation report) the first three derivatives are bounded the same way:
//! `|c^(j)| <= sqrt(P_j^2 + Z_j^2)`, with `P_1` the plan frame's stretch,
//! `P_2 = sup |k|`, `P_3 = sup |k'| + sup k^2` (in arc length
//! `p''' = k' n - k^2 t`; a chain's parametric piece is bounded through its
//! own derivative suprema), and `Z_j` from
//! [`crate::elevation`]'s per-law closed forms. They need `C^2` within the
//! span, so curvature seams of the plan are named as breaks for `k >= 2`.

use axiolid_contracts::{GeomError, GeomResult};
use axiolid_core::{Interval, Scalar, Vec3};
use axiolid_curve::{Banked3, Curve2, Curve3, Elevated3};

use crate::arc_length::{plan_derivative, plan_second_derivative};
use crate::bound::{frame_stretch2, CurveDerivativeBounds, ROUNDING};

fn invalid(detail: &str) -> GeomError {
    GeomError::InvalidInput(detail.to_owned())
}

/// Plan length an elevated curve's parameter runs over: unbounded for a
/// line, one full turn for a circle, the stored length of an intrinsic
/// curve or a chain; `None` for a plan that cannot carry an elevation law
/// or has no positive finite length.
#[must_use]
pub fn elevated_span(curve: &Elevated3) -> Option<Scalar> {
    let span = match curve.plan.as_ref() {
        Curve2::Line(_) => return Some(Scalar::INFINITY),
        Curve2::Circle(c) => core::f64::consts::TAU * c.radius,
        Curve2::Intrinsic(i) => i.length,
        Curve2::Chain(c) => c.length()?,
        _ => return None,
    };
    (span.is_finite() && span > 0.0).then_some(span)
}

/// Domain of an elevated curve in plan distance: `[0, L]`, `L` from
/// [`elevated_span`] (infinite for a line plan), or the empty interval.
pub(crate) fn elevated_domain(curve: &Elevated3) -> Interval {
    elevated_span(curve).map_or(
        Interval {
            start: 0.0,
            end: 0.0,
        },
        |end| Interval { start: 0.0, end },
    )
}

/// Derivative of an elevated curve in plan distance: `(p', z')`. Not unit:
/// its length is `sqrt(1 + g^2)` on an orthonormal plan frame.
///
/// # Errors
///
/// The plan's and the elevation law's refusals, by name.
pub fn elevated_derivative(curve: &Elevated3, d: Scalar) -> GeomResult<Vec3> {
    if !d.is_finite() {
        return Err(invalid("plan distance must be finite"));
    }
    let plan = plan_derivative(&curve.plan, d)?;
    let grade = crate::elevation::elevation_grade(&curve.elevation, d)?;
    let out = Vec3::new(plan.x, plan.y, grade);
    if out.is_finite() {
        Ok(out)
    } else {
        Err(invalid("elevated curve derivative is not finite"))
    }
}

/// Second derivative of an elevated curve in plan distance:
/// `(k n, z'')`, the plan's curvature vector and the profile's `z''`.
///
/// # Errors
///
/// As [`elevated_derivative`].
pub fn elevated_second_derivative(curve: &Elevated3, d: Scalar) -> GeomResult<Vec3> {
    if !d.is_finite() {
        return Err(invalid("plan distance must be finite"));
    }
    let plan = plan_second_derivative(&curve.plan, d)?;
    let bend = crate::elevation::elevation_second(&curve.elevation, d)?;
    let out = Vec3::new(plan.x, plan.y, bend);
    if out.is_finite() {
        Ok(out)
    } else {
        Err(invalid("elevated curve second derivative is not finite"))
    }
}

/// `(sup |p'|, sup |p''|, sup |p'''|)` of the plan over `[lo, hi]` in plan
/// distance, scaled by its frame's stretch; `p'''` is infinite across a
/// curvature seam. `None` for a plan family or span not bounded here.
fn plan_bounds(plan: &Curve2, lo: Scalar, hi: Scalar) -> Option<[Scalar; 3]> {
    let slack = |length: Scalar| 1e-9 * length.max(1.0);
    let out = match plan {
        Curve2::Line(line) => {
            let length = line.direction.length();
            if !(length.is_finite() && length > 0.0) {
                return None;
            }
            [1.0, 0.0, 0.0]
        }
        Curve2::Circle(c) => {
            if !(c.radius.is_finite() && c.radius > 0.0) {
                return None;
            }
            let stretch = frame_stretch2(&c.frame);
            let k = 1.0 / c.radius;
            [stretch, stretch * k, stretch * k * k]
        }
        Curve2::Intrinsic(i) => {
            if !(lo >= 0.0 && hi <= i.length + slack(i.length)) {
                return None;
            }
            let stretch = frame_stretch2(&i.start);
            let k = crate::chain::curvature_bound(&i.curvature, lo, hi)?;
            let seamed = i
                .curvature
                .seams_within(i.length)
                .iter()
                .any(|&seam| seam > lo && seam < hi);
            let rate = if seamed {
                Scalar::INFINITY
            } else {
                crate::chain::curvature_bound(&i.curvature.derivative(), lo, hi)?
            };
            [stretch, stretch * k, stretch * (rate + k * k)]
        }
        Curve2::Chain(chain) => {
            let (p2, p3) = crate::chain::chain_plan_bounds(chain, lo, hi)?;
            [frame_stretch2(&chain.start), p2, p3]
        }
        _ => return None,
    };
    (out[0].is_finite() && out[1].is_finite()).then_some(out)
}

/// Where an elevated curve may fail to be `C^k` (`k >= 1`), ascending: a
/// chain plan's joins and every seam of the elevation law (the grade may
/// jump there), and for `k >= 2` the curvature seams of an intrinsic plan
/// or of a chain's intrinsic pieces.
#[must_use]
pub fn elevated_breaks(curve: &Elevated3, k: usize) -> Vec<Scalar> {
    if k == 0 {
        return Vec::new();
    }
    let mut out = crate::bound::continuity_breaks2(&curve.plan, 1);
    out.extend(crate::elevation::elevation_seams(&curve.elevation));
    if k >= 2 {
        match curve.plan.as_ref() {
            Curve2::Intrinsic(i) => out.extend(i.curvature.seams_within(i.length)),
            Curve2::Chain(chain) => out.extend(crate::chain::chain_curvature_seams(chain)),
            _ => {}
        }
    }
    out.retain(|b| b.is_finite());
    out.sort_by(Scalar::total_cmp);
    out.dedup();
    out
}

/// [`elevated_breaks`] of a banked curve's point path, plus the seams of
/// its pivot law (its rate may jump there). The cant law moves only the
/// section frame, not the point, so its seams are not named -- except
/// where a held-rail pivot piece derives the point from the cant (#279):
/// the cant's seams there are the pivot's.
#[must_use]
pub fn banked_breaks(curve: &Banked3, k: usize) -> Vec<Scalar> {
    if k == 0 {
        return Vec::new();
    }
    let mut out = elevated_breaks(&curve.base, k);
    out.extend(curve.pivot.seams());
    out.extend(crate::banked::derived_cant_seams(curve));
    out.sort_by(Scalar::total_cmp);
    out.dedup();
    out
}

fn inside(breaks: &[Scalar], lo: Scalar, hi: Scalar) -> bool {
    breaks.iter().any(|&b| b > lo && b < hi)
}

/// Certified 3D chord bound of an elevated curve over plan distances
/// `[a, b]`: `sqrt(P^2 + Z^2)`, `P = h^2 / 8 sup |k|` the plan's and `Z`
/// the profile's chord bound at equal plan distance (see the
/// [module documentation](self) for the derivation).
///
/// `None` -- unbounded, never a guess -- across a chain join or a profile
/// seam, and where either half is not bounded.
#[must_use]
pub fn elevated_chord_bound(curve: &Elevated3, a: Scalar, b: Scalar) -> Option<Scalar> {
    let (lo, hi) = (a.min(b), a.max(b));
    if !(lo.is_finite() && hi.is_finite()) {
        return None;
    }
    if lo == hi {
        return Some(0.0);
    }
    if inside(&elevated_breaks(curve, 1), lo, hi) {
        return None;
    }
    let h = hi - lo;
    let plan = h * h * 0.125 * plan_bounds(&curve.plan, lo, hi)?[1];
    let profile = crate::elevation::elevation_chord_bound(&curve.elevation, lo, hi)?;
    let bound = plan.hypot(profile) * (1.0 + ROUNDING);
    bound.is_finite().then_some(bound)
}

/// [`elevated_chord_bound`] of a banked curve's point path: the vertical
/// part is `z + e`, bounded by `Z + h^2 / 8 sup |e''|`.
#[must_use]
pub fn banked_chord_bound(curve: &Banked3, a: Scalar, b: Scalar) -> Option<Scalar> {
    let (lo, hi) = (a.min(b), a.max(b));
    if !(lo.is_finite() && hi.is_finite()) {
        return None;
    }
    if lo == hi {
        return Some(0.0);
    }
    if lo < 0.0 || inside(&banked_breaks(curve, 1), lo, hi) {
        return None;
    }
    let h = hi - lo;
    let taylor = h * h * 0.125;
    let plan = taylor * plan_bounds(&curve.base.plan, lo, hi)?[1];
    let profile = crate::elevation::elevation_chord_bound(&curve.base.elevation, lo, hi)?;
    let pivot = taylor * crate::banked::pivot_bounds(curve, lo, hi)?[1];
    let bound = plan.hypot(profile + pivot) * (1.0 + ROUNDING);
    bound.is_finite().then_some(bound)
}

/// Certified suprema of `|c'|`, `|c''|`, `|c'''|` of an elevated curve in
/// plan distance over `[a, b]`: `sqrt(P_j^2 + Z_j^2)` (module
/// documentation). `None` where the curve may not be `C^2` strictly
/// inside the span, or where either half is not bounded.
#[must_use]
pub fn elevated_derivative_bounds(
    curve: &Elevated3,
    a: Scalar,
    b: Scalar,
) -> Option<CurveDerivativeBounds> {
    let (lo, hi) = (a.min(b), a.max(b));
    if !(lo.is_finite() && hi.is_finite()) || inside(&elevated_breaks(curve, 2), lo, hi) {
        return None;
    }
    let plan = plan_bounds(&curve.plan, lo, hi)?;
    let profile = crate::elevation::profile_bounds(&curve.elevation, lo, hi)?;
    combined(plan, profile)
}

/// [`elevated_derivative_bounds`] of a banked curve's point path, the
/// pivot's derivatives added to the profile's.
#[must_use]
pub fn banked_derivative_bounds(
    curve: &Banked3,
    a: Scalar,
    b: Scalar,
) -> Option<CurveDerivativeBounds> {
    let (lo, hi) = (a.min(b), a.max(b));
    if !(lo.is_finite() && hi.is_finite()) || lo < 0.0 || inside(&banked_breaks(curve, 2), lo, hi) {
        return None;
    }
    let plan = plan_bounds(&curve.base.plan, lo, hi)?;
    let profile = crate::elevation::profile_bounds(&curve.base.elevation, lo, hi)?;
    let pivot = crate::banked::pivot_bounds(curve, lo, hi)?;
    combined(
        plan,
        [
            profile[0] + pivot[0],
            profile[1] + pivot[1],
            profile[2] + pivot[2],
        ],
    )
}

fn combined(plan: [Scalar; 3], vertical: [Scalar; 3]) -> Option<CurveDerivativeBounds> {
    let bound = |j: usize| plan[j].hypot(vertical[j]) * (1.0 + ROUNDING);
    let out = CurveDerivativeBounds {
        first: bound(0),
        second: bound(1),
        third: bound(2),
    };
    (out.first.is_finite() && out.second.is_finite() && out.third.is_finite()).then_some(out)
}

/// Whether [`crate::curve::flatten3`] certifies an elevated curve's chords:
/// its plan is a line, a circle, an intrinsic curve whose curvature law is
/// bounded, or a chain whose pieces all are
/// ([`crate::chain::chain_certifies`]), and every piece of its profile is
/// a polynomial or a circular arc. A circular arc's bound exists wherever
/// the arc has heights, so no span the curve evaluates on is left to the
/// midpoint sagitta.
pub(crate) fn certifies_elevated(curve: &Elevated3) -> bool {
    fn closed_form(law: &axiolid_curve::ElevationLaw) -> bool {
        use axiolid_curve::ElevationLaw as L;
        match law {
            L::Polynomial { .. } | L::CircularArc { .. } => true,
            L::Piecewise { laws, .. } => laws.iter().all(closed_form),
            _ => false,
        }
    }
    let plan = match curve.plan.as_ref() {
        Curve2::Line(_) | Curve2::Circle(_) => true,
        Curve2::Intrinsic(i) => {
            crate::chain::curvature_bound(&i.curvature, 0.0, i.length.max(0.0)).is_some()
        }
        Curve2::Chain(chain) => crate::chain::chain_certifies(chain),
        _ => false,
    };
    plan && curve.elevation.is_well_formed() && closed_form(&curve.elevation)
}

/// Where an elevated or banked curve turns a corner: the seams (of its
/// profile, and of a banked curve's pivot) across which the grade of its
/// point path jumps by more than `angle` radians, ascending.
///
/// A chain plan is tangent-continuous at its joins by construction, so
/// only the vertical part can kink. Read from the laws on either side of
/// each seam, not from nearby samples. Empty for every other family.
///
/// # Errors
///
/// The laws' refusals at a seam, by name.
pub fn grade_corners3(curve: &Curve3, angle: Scalar) -> GeomResult<Vec<Scalar>> {
    let (elevated, banked) = match curve {
        Curve3::Elevated(e) => (e, None),
        Curve3::Banked(b) => (&b.base, Some(b)),
        _ => return Ok(Vec::new()),
    };
    let mut seams = crate::elevation::elevation_seams(&elevated.elevation);
    if let Some(banked) = banked {
        seams.extend(banked.pivot.seams());
        seams.extend(crate::banked::derived_cant_seams(banked));
    }
    seams.sort_by(Scalar::total_cmp);
    seams.dedup();
    let mut out = Vec::new();
    for seam in seams {
        let mut before = crate::elevation::elevation_grade_before(&elevated.elevation, seam)?;
        let mut after = crate::elevation::elevation_grade(&elevated.elevation, seam)?;
        if let Some(banked) = banked {
            let rate = |before: bool| {
                crate::banked::pivot_jet(banked, seam, before)
                    .map(|jet| jet[1])
                    .ok_or_else(|| invalid("banked curve: pivot law has no rate at a seam"))
            };
            before += rate(true)?;
            after += rate(false)?;
        }
        if (before.atan() - after.atan()).abs() > angle {
            out.push(seam);
        }
    }
    Ok(out)
}
