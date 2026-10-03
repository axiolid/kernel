//! Heights, grades and chord bounds of elevation laws, including the profile
//! given by its own curvature, which has no closed form (#238).
//!
//! # Closed forms are delegated, the intrinsic profile is integrated
//!
//! A polynomial or a circular vertical arc answers
//! [`ElevationLaw::height_at`] and [`ElevationLaw::grade_at`] exactly, and
//! this module calls those unchanged. An [`ElevationLaw::Intrinsic`]
//! profile is a planar curve in `(d, z)` whose curvature runs along its own
//! arc length `s`, so a height at PLAN distance `d` needs two numerical
//! steps, both stated here:
//!
//! 1. the position `(d(s), z(s))` is [`intrinsic_point`] on the profile
//!    curve, the Gauss-Legendre quadrature pinned against the Fresnel
//!    closed form (ADR 0060);
//! 2. `d(s) = d` is inverted by a Newton iteration safeguarded by a
//!    bracket, since `d'(s) = cos t(s)` is exact from the closed-form
//!    heading.
//!
//! # Accuracy contract
//!
//! The inversion stops once `|d(s) - d| <= INVERSION_TOLERANCE * max(1, d)`
//! and corrects the height to first order by the residual times the grade.
//! A reported height therefore carries the quadrature error of
//! [`intrinsic_point`] plus at most `|grade| * INVERSION_TOLERANCE *
//! max(1, d)`; a grade carries the curvature times that same arc-length
//! residual, divided by `cos^2 t`. Measured against an independent 50-digit
//! integration and root, a 150 m sag-to-crest clothoid agrees to within
//! `4e-15` m in height and `2e-17` in grade; the tests pin `1e-13` and
//! `1e-15`.
//!
//! # Where the plan distance stops increasing
//!
//! `d(s)` is monotone only while the profile's direction stays strictly
//! inside `(-pi/2, pi/2)`. That is certified, never assumed: over the
//! arc-length span that brackets the inverse, the heading is bounded by its
//! closed-form value at a panel start plus the total variation of curvature
//! over the panel, and panels are bisected until the bound stays below
//! vertical. The root is then unique in its bracket. A profile that
//! reaches vertical first, or that cannot be certified within the bisection
//! budget, is refused by name rather than read on a branch of `d(s)` that
//! runs backwards.

use core::f64::consts::FRAC_PI_2;

use axiolid_contracts::{GeomError, GeomResult};
use axiolid_core::{Frame2, Point2, Scalar, Vec2};
use axiolid_curve::{CurvatureLaw, ElevationLaw, Intrinsic2};

use crate::arc_length::intrinsic_point;

/// Relative stopping tolerance of the plan-distance inversion: the solve
/// ends once `|d(s) - d| <= INVERSION_TOLERANCE * max(1, d)`.
pub const INVERSION_TOLERANCE: Scalar = 1e-12;

/// Newton steps before the inversion is refused. Quadratic convergence from
/// `s = d` needs fewer than ten on any profile that stays below vertical;
/// the rest is bisection headroom for a nearly vertical one.
const MAX_INVERSION_STEPS: usize = 100;

/// Bracket enlargements before an inverse is refused as unbracketable.
const MAX_BRACKET_ROUNDS: usize = 64;

/// Bisection depth of the below-vertical certificate.
const MAX_CERTIFY_DEPTH: u32 = 24;

/// Panels the below-vertical certificate may inspect in total.
const MAX_CERTIFY_PANELS: usize = 1 << 16;

/// Least fraction of `cos t` at a panel start that its certified floor may
/// keep before the panel is bisected.
const ACCEPT: Scalar = 15.0 / 16.0;

fn invalid(detail: &str) -> GeomError {
    GeomError::InvalidInput(detail.to_owned())
}

/// Height at plan distance `distance`, for every law in the family.
///
/// Closed-form laws answer through [`ElevationLaw::height_at`]; an
/// [`ElevationLaw::Intrinsic`] piece, alone or inside a
/// [`ElevationLaw::Piecewise`], is integrated and inverted as the module
/// documentation states.
///
/// # Errors
///
/// [`GeomError::InvalidInput`] for a non-finite distance, a malformed law,
/// a distance outside a law's domain (past the vertical of a circular arc,
/// before the start of an intrinsic profile, or past where it turns
/// vertical), and [`GeomError::BudgetExceeded`] when the below-vertical
/// certificate or the inversion does not settle within its budget.
pub fn elevation_height(law: &ElevationLaw, distance: Scalar) -> GeomResult<Scalar> {
    match leaf(law, distance)? {
        (
            ElevationLaw::Intrinsic {
                height,
                grade,
                curvature,
            },
            local,
        ) => Ok(intrinsic_reading(*height, *grade, curvature, local)?.height),
        // Every closed-form law, including any later one: the curve crate's
        // match over its own enum is exhaustive and refuses what it cannot
        // answer, so nothing is guessed here.
        (closed, local) => closed
            .height_at(local)
            .ok_or_else(|| invalid("elevation law has no height at that distance")),
    }
}

/// Grade `dz/dd` at plan distance `distance`, for every law in the family.
///
/// # Errors
///
/// As [`elevation_height`].
pub fn elevation_grade(law: &ElevationLaw, distance: Scalar) -> GeomResult<Scalar> {
    match leaf(law, distance)? {
        (
            ElevationLaw::Intrinsic {
                height,
                grade,
                curvature,
            },
            local,
        ) => Ok(intrinsic_reading(*height, *grade, curvature, local)?.grade),
        (closed, local) => closed
            .grade_at(local)
            .ok_or_else(|| invalid("elevation law has no grade at that distance")),
    }
}

/// The innermost piece at `distance`, or a named refusal.
fn leaf(law: &ElevationLaw, distance: Scalar) -> GeomResult<(&ElevationLaw, Scalar)> {
    if !distance.is_finite() {
        return Err(invalid("plan distance must be finite"));
    }
    law.piece_at(distance)
        .ok_or_else(|| invalid("elevation law has no piece at that distance"))
}

/// Height and grade read off an intrinsic profile.
struct Reading {
    height: Scalar,
    grade: Scalar,
}

/// The profile as a plane curve: origin `(0, height)`, `x` along the start
/// direction `atan(grade)`.
fn profile_curve(
    height: Scalar,
    grade: Scalar,
    curvature: &CurvatureLaw,
    span: Scalar,
) -> Intrinsic2 {
    let norm = grade.hypot(1.0);
    let (sin0, cos0) = (grade / norm, 1.0 / norm);
    Intrinsic2::new(
        Frame2 {
            origin: Point2::new(0.0, height),
            x: Vec2::new(cos0, sin0),
            y: Vec2::new(-sin0, cos0),
        },
        curvature.clone(),
        span,
    )
}

/// Check the stored numbers of an intrinsic profile.
fn checked(height: Scalar, grade: Scalar, curvature: &CurvatureLaw) -> GeomResult<()> {
    if !(height.is_finite() && grade.is_finite()) {
        return Err(invalid(
            "intrinsic profile needs a finite start height and grade",
        ));
    }
    if !curvature.is_well_formed() {
        return Err(invalid("intrinsic profile curvature law is malformed"));
    }
    Ok(())
}

/// Height and grade of an intrinsic profile at plan distance `d >= 0`.
fn intrinsic_reading(
    height: Scalar,
    grade: Scalar,
    curvature: &CurvatureLaw,
    d: Scalar,
) -> GeomResult<Reading> {
    checked(height, grade, curvature)?;
    if d < 0.0 {
        return Err(invalid(
            "an intrinsic profile is defined forward from its start only",
        ));
    }
    if d == 0.0 {
        return Ok(Reading { height, grade });
    }
    let theta0 = grade.atan();
    let curve = profile_curve(height, grade, curvature, d);
    let upper = bracket(&curve, theta0, d)?;

    // d(s) <= s, so the root is at or past `d`; `upper` has d(upper) >= d.
    let (mut lo, mut hi) = (d, upper);
    let mut s = d;
    let tolerance = INVERSION_TOLERANCE * d.max(1.0);
    for _ in 0..MAX_INVERSION_STEPS {
        let point = intrinsic_point(&curve, s)?;
        let theta = theta0 + heading(&curve, s)?;
        let residual = point.x - d;
        if residual.abs() <= tolerance {
            let slope = theta.tan();
            // z(d) = z(s) - tan t * (d(s) - d) to first order.
            return Ok(Reading {
                height: point.y - slope * residual,
                grade: slope,
            });
        }
        if residual < 0.0 {
            lo = s;
        } else {
            hi = s;
        }
        // d'(s) = cos t(s), exact from the closed-form heading.
        let mut next = s - residual / theta.cos();
        if !(next > lo && next < hi) {
            next = 0.5 * (lo + hi);
        }
        if next == s {
            break;
        }
        s = next;
    }
    Err(GeomError::BudgetExceeded {
        resource: "intrinsic profile plan-distance inversion",
    })
}

/// Closed-form heading turned through `[0, s]` of the profile curve.
fn heading(curve: &Intrinsic2, s: Scalar) -> GeomResult<Scalar> {
    curve
        .heading_at(s)
        .ok_or_else(|| invalid("curvature law does not integrate to that arc length"))
}

/// An arc length `upper` with `d(upper) >= d` (to the inversion tolerance),
/// over whose `[0, upper]` the profile is certified below vertical, so
/// `d(s)` increases there and the root of `d(s) = d` in `[d, upper]` is
/// unique.
///
/// Steps forward as Newton does from below, `upper += (d - d(upper)) /
/// cos t(upper)` with `d(upper)` from the quadrature, certifying only each
/// extension. An extension that runs into the vertical is halved, since the
/// target may still lie before it; one that cannot shrink further refuses
/// with the certificate's reason.
fn bracket(curve: &Intrinsic2, theta0: Scalar, d: Scalar) -> GeomResult<Scalar> {
    cos_floor(curve, theta0, 0.0, d)?;
    let tolerance = INVERSION_TOLERANCE * d.max(1.0);
    let mut upper = d;
    for _ in 0..MAX_BRACKET_ROUNDS {
        let reached = intrinsic_point(curve, upper)?.x;
        if reached >= d - tolerance {
            return Ok(upper);
        }
        let cos = (theta0 + heading(curve, upper)?).cos();
        // Slightly long, so a concave d(s) -- approached from below without
        // ever crossing -- is still overtaken.
        let mut step = (d - reached) / cos * (1.0 + 1e-6) + Scalar::EPSILON * upper.max(1.0);
        loop {
            match cos_floor(curve, theta0, upper, upper + step) {
                Ok(_) => break,
                Err(error) => {
                    step *= 0.5;
                    if step <= Scalar::EPSILON * upper.max(1.0) {
                        return Err(error);
                    }
                }
            }
        }
        upper += step;
    }
    Err(GeomError::BudgetExceeded {
        resource: "intrinsic profile plan-distance bracket",
    })
}

/// Certified lower bound of `cos t(s)` over `[lo, hi]`, where
/// `t = theta0 + heading`: the certificate that the profile stays below
/// vertical there.
///
/// On a panel `[a, b]` the heading lies within `t(a) +- V`, `V` the total
/// variation bound of curvature over the panel, so `cos t >= cos(|t(a)| +
/// V)` there. A panel is accepted once that reach stays below `pi/2` and
/// keeps at least [`ACCEPT`] of `cos t(a)`, which keeps the floor close to
/// the true least `cos`; otherwise it is bisected, and at the depth limit it
/// is accepted on the first condition alone. A panel start at or past
/// vertical refuses, and so does a panel at the depth limit whose reach
/// still touches vertical: the profile meets vertical within the
/// certificate's resolution there.
fn cos_floor(curve: &Intrinsic2, theta0: Scalar, lo: Scalar, hi: Scalar) -> GeomResult<Scalar> {
    let vertical = || invalid("the profile turns vertical before reaching that plan distance");
    let mut floor: Scalar = 1.0;
    let mut stack = vec![(lo, hi, 0_u32)];
    let mut inspected = 0_usize;
    while let Some((a, b, depth)) = stack.pop() {
        inspected += 1;
        if inspected > MAX_CERTIFY_PANELS {
            return Err(GeomError::BudgetExceeded {
                resource: "intrinsic profile below-vertical certificate",
            });
        }
        let theta_a = theta0 + heading(curve, a)?;
        if theta_a.is_nan() || theta_a.abs() >= FRAC_PI_2 {
            return Err(vertical());
        }
        let shifted = curve
            .curvature
            .shifted(a)
            .ok_or_else(|| invalid("curvature law cannot be read from that arc length"))?;
        let variation = Intrinsic2::new(curve.start, shifted, b - a)
            .turning_variation_bound(b - a)
            .ok_or_else(|| invalid("curvature law does not integrate over the span"))?;
        let reach = theta_a.abs() + variation;
        let panel_floor = reach.cos();
        let below_vertical = reach < FRAC_PI_2 && panel_floor > 0.0;
        if below_vertical && (panel_floor >= ACCEPT * theta_a.cos() || depth >= MAX_CERTIFY_DEPTH) {
            floor = floor.min(panel_floor);
            continue;
        }
        if depth >= MAX_CERTIFY_DEPTH {
            return Err(vertical());
        }
        let mid = 0.5 * (a + b);
        stack.push((mid, b, depth + 1));
        stack.push((a, mid, depth + 1));
    }
    Ok(floor)
}

/// Certified bound on how far an elevation law strays from the chord of
/// its heights over the plan-distance span `[a, b]`:
/// `sup |z(d) - L(d)|`, `L` the straight line through `(a, z(a))` and
/// `(b, z(b))`.
///
/// This is the vertical half of flattening an elevated curve:
/// `(b - a)^2 / 8 * sup |z''|` with `z'' = k / cos^3 t` for a curve in
/// `(d, z)` of curvature `k`, bounded in closed form per law:
///
/// - polynomial: `sum i (i - 1) |c_i| m^(i - 2)`, `m = max(|a|, |b|)`;
/// - circular arc: `1 / (|R| cos^3 t)`, `cos t` least at the span end with
///   the larger `|sin t|`, since `sin t` is linear in `d`;
/// - intrinsic: the curvature supremum over a certified arc-length span
///   containing `[s(a), s(b)]`, over the cube of the certified `cos t`
///   floor there.
///
/// `None` -- unbounded, never a guess -- when the span is not finite, a
/// piecewise seam lies strictly inside it (the grade may jump there), the
/// span leaves a law's domain, or the law or its curvature law is a family
/// this function does not bound.
#[must_use]
pub fn elevation_chord_bound(law: &ElevationLaw, a: Scalar, b: Scalar) -> Option<Scalar> {
    if !(a.is_finite() && b.is_finite()) {
        return None;
    }
    let (lo, hi) = if a <= b { (a, b) } else { (b, a) };
    if lo == hi {
        return Some(0.0);
    }
    let second = second_derivative_bound(law, lo, hi)?;
    let h = hi - lo;
    let bound = h * h * 0.125 * second;
    // A few ulps of slack for the rounding in the bound's own arithmetic.
    Some(bound * (1.0 + 1e-12)).filter(|b| b.is_finite())
}

/// `sup |z''|` over `[lo, hi]`, or `None` when not bounded.
fn second_derivative_bound(law: &ElevationLaw, lo: Scalar, hi: Scalar) -> Option<Scalar> {
    match law {
        ElevationLaw::Polynomial { coefficients } => {
            let m = lo.abs().max(hi.abs());
            Some(
                coefficients
                    .iter()
                    .enumerate()
                    .skip(2)
                    .map(|(i, c)| {
                        let i_f = i as Scalar;
                        i_f * (i_f - 1.0) * c.abs() * m.powi(i as i32 - 2)
                    })
                    .sum(),
            )
            .filter(|s: &Scalar| s.is_finite())
        }
        ElevationLaw::Piecewise { breaks, laws } => {
            if laws.len() != breaks.len() + 1 {
                return None;
            }
            // The piece owning the span; it must contain both ends.
            let mid = 0.5 * (lo + hi);
            let index = breaks.partition_point(|b| *b <= mid);
            let start = if index == 0 { 0.0 } else { breaks[index - 1] };
            let lower = if index == 0 {
                Scalar::NEG_INFINITY
            } else {
                start
            };
            let upper = breaks.get(index).copied().unwrap_or(Scalar::INFINITY);
            if lo < lower || hi > upper {
                return None;
            }
            second_derivative_bound(laws.get(index)?, lo - start, hi - start)
        }
        ElevationLaw::CircularArc { grade, radius, .. } => {
            if !(grade.is_finite() && radius.is_finite() && *radius != 0.0) {
                return None;
            }
            let sin0 = grade / grade.hypot(1.0);
            let reach = (sin0 + lo / radius).abs().max((sin0 + hi / radius).abs());
            if reach.is_nan() || reach >= 1.0 {
                return None;
            }
            let cos = ((1.0 - reach) * (1.0 + reach)).sqrt();
            Some(1.0 / (radius.abs() * cos * cos * cos)).filter(|s| s.is_finite())
        }
        ElevationLaw::Intrinsic {
            height,
            grade,
            curvature,
        } => {
            checked(*height, *grade, curvature).ok()?;
            if lo < 0.0 {
                return None;
            }
            // s(lo) >= lo because d(s) <= s, and `bracket` returns an arc
            // length past s(hi); the span read is [s(lo), s(hi)] inside
            // [lo, upper], and the floor is certified over that.
            let theta0 = grade.atan();
            let curve = profile_curve(*height, *grade, curvature, hi);
            // Bracket a little past `hi`, so the bracket's own tolerance
            // cannot leave s(hi) outside it.
            let past = hi + 2.0 * INVERSION_TOLERANCE * hi.max(1.0);
            let upper = bracket(&curve, theta0, past).ok()?;
            let floor = cos_floor(&curve, theta0, lo, upper).ok()?;
            let kappa = curvature_sup(&curvature.shifted(lo)?, upper - lo)?;
            Some(kappa / (floor * floor * floor)).filter(|s| s.is_finite())
        }
        // A later elevation family this function does not know: unbounded.
        _ => None,
    }
}

/// Upper bound on `|k(s)|` over `[0, span]`, `span >= 0`.
fn curvature_sup(law: &CurvatureLaw, span: Scalar) -> Option<Scalar> {
    let polynomial = |coefficients: &[Scalar]| -> Scalar {
        coefficients
            .iter()
            .enumerate()
            .map(|(i, c)| c.abs() * span.powi(i as i32))
            .sum()
    };
    let sup = match law {
        CurvatureLaw::Constant { curvature } => curvature.abs(),
        // A linear law (the clothoid) peaks at an end: exact.
        CurvatureLaw::Polynomial { coefficients } if coefficients.len() <= 2 => {
            let start = coefficients.first().copied().unwrap_or(0.0);
            let rate = coefficients.get(1).copied().unwrap_or(0.0);
            start.abs().max((start + rate * span).abs())
        }
        CurvatureLaw::Polynomial { coefficients } => polynomial(coefficients),
        CurvatureLaw::Sinusoid {
            mean, amplitude, ..
        } => mean.abs() + amplitude.abs(),
        CurvatureLaw::Composite {
            polynomial: terms,
            harmonics,
        } => polynomial(terms) + harmonics.iter().map(|h| h.amplitude.abs()).sum::<Scalar>(),
        CurvatureLaw::Piecewise { breaks, laws } => {
            if !law.is_well_formed() {
                return None;
            }
            let mut sup: Scalar = 0.0;
            let mut start = 0.0;
            for (index, piece) in laws.iter().enumerate() {
                let end = breaks.get(index).copied().unwrap_or(span).min(span);
                sup = sup.max(curvature_sup(piece, (end - start).max(0.0))?);
                if end >= span {
                    break;
                }
                start = end;
            }
            sup
        }
        // A later curvature family this function does not know: unbounded.
        _ => return None,
    };
    Some(sup).filter(|s| s.is_finite())
}
