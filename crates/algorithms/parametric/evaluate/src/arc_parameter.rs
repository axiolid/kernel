//! Arc length along any evaluable curve, and its inverse.
//!
//! Most curves have no closed-form arc length: a cubic parabola
//! `y = A x^3` needs `int sqrt(1 + 9 A^2 x^4) dx`, an elliptic integral, and
//! a B-spline or an ellipse is no better. Yet some data names a position by
//! distance along such a curve: an arc-length trim selector, or a chain
//! piece read by arc length. This module discharges both with quadrature of
//! the speed `|c'(t)|` and a safeguarded Newton root find, the same bargain
//! [`crate::arc_length`] makes for intrinsic curves: the stored value stays
//! exact, only reading it needs numerical work.
//!
//! # Accuracy contract
//!
//! - [`arc_length2`] and [`arc_length3`] integrate the speed by adaptive
//!   8-point Gauss-Legendre, panels split at every knot and polyline vertex
//!   (the speed is only piecewise smooth there) and bisected until a panel
//!   and its two halves agree to the panel's share of
//!   `ARC_LENGTH_TOLERANCE * max(1, length)`. That agreement is an error
//!   ESTIMATE, as QUADPACK's is, not a proof; on a smooth speed the halves
//!   are far more accurate than the estimate says.
//! - [`parameter_at_arc_length2`] and [`parameter_at_arc_length3`] return a
//!   parameter whose arc length from `from` matches the request to the same
//!   tolerance, or a parameter bracket collapsed to rounding.
//! - Lines, intrinsic curves and chains are exact: their parameter is a
//!   multiple of arc length, so no quadrature is run.
//!
//! A request past the end of a bounded curve, a start parameter outside its
//! domain, a curve whose speed does not evaluate, and a quadrature that
//! needs more than [`MAX_PANELS`] panels are refused by name.

use axiolid_contracts::{GeomError, GeomResult};
use axiolid_core::{Interval, Scalar};
use axiolid_curve::{Curve2, Curve3};

use crate::curve::{derivative2, derivative3, domain2, domain3};

/// Relative tolerance of every arc length and arc-length inverse this module
/// returns: the absolute tolerance for a length `l` is
/// `ARC_LENGTH_TOLERANCE * max(1, |l|)`, in the curve's length unit.
pub const ARC_LENGTH_TOLERANCE: Scalar = 1e-12;

/// Most quadrature panels one arc length may use before it refuses.
pub const MAX_PANELS: usize = 1 << 14;

/// Root-find iterations before the inverse refuses.
const MAX_ITERATIONS: usize = 200;

/// Nodes of the 8-point Gauss-Legendre rule on `[-1, 1]`.
const NODES: [Scalar; 8] = [
    -0.960_289_856_497_536_2,
    -0.796_666_477_413_626_7,
    -0.525_532_409_916_328_9,
    -0.183_434_642_495_649_8,
    0.183_434_642_495_649_8,
    0.525_532_409_916_328_9,
    0.796_666_477_413_626_7,
    0.960_289_856_497_536_2,
];

/// Weights matching [`NODES`].
const WEIGHTS: [Scalar; 8] = [
    0.101_228_536_290_376_3,
    0.222_381_034_453_374_5,
    0.313_706_645_877_887_3,
    0.362_683_783_378_361_9,
    0.362_683_783_378_361_9,
    0.313_706_645_877_887_3,
    0.222_381_034_453_374_5,
    0.101_228_536_290_376_3,
];

fn invalid(detail: String) -> GeomError {
    GeomError::InvalidInput(detail)
}

/// How a curve family's parameter relates to its arc length.
enum Family {
    /// Arc length is `speed * (to - from)` exactly.
    Proportional(Scalar),
    /// Arc length needs quadrature; the parameter may run over `domain`
    /// (`None`: anywhere, as for a periodic conic) and the speed is only
    /// piecewise smooth across `breaks`.
    Numeric {
        domain: Option<(Scalar, Scalar)>,
        breaks: Vec<Scalar>,
    },
}

/// Every knot of a B-spline, every polyline vertex, every chain join.
const ALL_BREAKS: usize = u16::MAX as usize + 1;

fn bounded(domain: Interval) -> Option<(Scalar, Scalar)> {
    Some((domain.start.min(domain.end), domain.start.max(domain.end)))
}

fn family2(curve: &Curve2) -> GeomResult<Family> {
    Ok(match curve {
        Curve2::Line(line) => {
            let speed = line.direction.length();
            if !(speed.is_finite() && speed > 0.0) {
                return Err(invalid("line has no direction, so no arc length".into()));
            }
            Family::Proportional(speed)
        }
        Curve2::Intrinsic(_) | Curve2::Chain(_) => Family::Proportional(1.0),
        Curve2::Circle(_) | Curve2::Ellipse(_) | Curve2::Sinusoid(_) => Family::Numeric {
            domain: None,
            breaks: Vec::new(),
        },
        _ => Family::Numeric {
            domain: bounded(domain2(curve)),
            breaks: crate::bound::continuity_breaks2(curve, ALL_BREAKS),
        },
    })
}

fn family3(curve: &Curve3) -> GeomResult<Family> {
    Ok(match curve {
        Curve3::Line(line) => {
            let speed = line.direction.length();
            if !(speed.is_finite() && speed > 0.0) {
                return Err(invalid("line has no direction, so no arc length".into()));
            }
            Family::Proportional(speed)
        }
        Curve3::Intrinsic(_) => Family::Proportional(1.0),
        Curve3::Circle(_) | Curve3::Ellipse(_) => Family::Numeric {
            domain: None,
            breaks: Vec::new(),
        },
        _ => Family::Numeric {
            domain: bounded(domain3(curve)),
            breaks: crate::bound::continuity_breaks3(curve, ALL_BREAKS),
        },
    })
}

/// Signed arc length of a 2D curve from parameter `from` to `to`: negative
/// when `to < from`. See the [module documentation](self) for the
/// tolerance.
pub fn arc_length2(curve: &Curve2, from: Scalar, to: Scalar) -> GeomResult<Scalar> {
    let speed = |t| derivative2(curve, t).map(|d| d.length());
    arc_length(&family2(curve)?, &speed, from, to)
}

/// Signed arc length of a 3D curve from parameter `from` to `to`; see
/// [`arc_length2`].
pub fn arc_length3(curve: &Curve3, from: Scalar, to: Scalar) -> GeomResult<Scalar> {
    let speed = |t| derivative3(curve, t).map(|d| d.length());
    arc_length(&family3(curve)?, &speed, from, to)
}

/// The parameter at signed arc length `length` from parameter `from` along
/// a 2D curve: in increasing parameter for a positive length, decreasing
/// for a negative one.
///
/// Refuses by name a length past the end of a bounded curve, a `from`
/// outside its domain, and non-finite input. See the
/// [module documentation](self) for the tolerance.
pub fn parameter_at_arc_length2(
    curve: &Curve2,
    from: Scalar,
    length: Scalar,
) -> GeomResult<Scalar> {
    let speed = |t| derivative2(curve, t).map(|d| d.length());
    parameter_at(&family2(curve)?, &speed, from, length)
}

/// [`parameter_at_arc_length2`] for a 3D curve.
pub fn parameter_at_arc_length3(
    curve: &Curve3,
    from: Scalar,
    length: Scalar,
) -> GeomResult<Scalar> {
    let speed = |t| derivative3(curve, t).map(|d| d.length());
    parameter_at(&family3(curve)?, &speed, from, length)
}

fn tolerance_for(length: Scalar) -> Scalar {
    ARC_LENGTH_TOLERANCE * length.abs().max(1.0)
}

fn check_finite(from: Scalar, other: Scalar, what: &str) -> GeomResult<()> {
    if from.is_finite() && other.is_finite() {
        Ok(())
    } else {
        Err(invalid(format!(
            "arc length needs a finite start parameter and {what}"
        )))
    }
}

fn inside(domain: Option<(Scalar, Scalar)>, t: Scalar) -> bool {
    domain.is_none_or(|(lo, hi)| {
        let slack = 1e-12 * (hi - lo).abs().max(1.0);
        t >= lo - slack && t <= hi + slack
    })
}

fn arc_length<F>(family: &Family, speed: &F, from: Scalar, to: Scalar) -> GeomResult<Scalar>
where
    F: Fn(Scalar) -> GeomResult<Scalar>,
{
    check_finite(from, to, "a finite end parameter")?;
    match family {
        Family::Proportional(rate) => Ok(rate * (to - from)),
        Family::Numeric { domain, breaks } => {
            if !inside(*domain, from) || !inside(*domain, to) {
                return Err(invalid(format!(
                    "arc length from parameter {from} to {to} leaves the curve's domain"
                )));
            }
            let (lo, hi) = (from.min(to), from.max(to));
            let magnitude = integrate(speed, breaks, lo, hi)?;
            Ok(if to < from { -magnitude } else { magnitude })
        }
    }
}

/// Integral of the speed over `[lo, hi]`, split at the breaks inside.
fn integrate<F>(speed: &F, breaks: &[Scalar], lo: Scalar, hi: Scalar) -> GeomResult<Scalar>
where
    F: Fn(Scalar) -> GeomResult<Scalar>,
{
    if hi <= lo {
        return Ok(0.0);
    }
    let mut cuts = vec![lo];
    cuts.extend(breaks.iter().copied().filter(|&t| t > lo && t < hi));
    cuts.push(hi);
    // A first, coarse pass sizes the absolute tolerance to the length.
    let mut coarse = 0.0;
    for w in cuts.windows(2) {
        coarse += rule(speed, w[0], w[1])?;
    }
    let tolerance = tolerance_for(coarse);
    let mut total = 0.0;
    let mut panels = 0usize;
    for w in cuts.windows(2) {
        let share = tolerance * (w[1] - w[0]) / (hi - lo);
        total += adaptive(speed, w[0], w[1], share, &mut panels)?;
    }
    if !total.is_finite() {
        return Err(invalid(
            "curve speed does not integrate to a finite arc length".into(),
        ));
    }
    Ok(total)
}

/// One 8-point Gauss-Legendre panel.
fn rule<F>(speed: &F, a: Scalar, b: Scalar) -> GeomResult<Scalar>
where
    F: Fn(Scalar) -> GeomResult<Scalar>,
{
    let half = 0.5 * (b - a);
    let mid = 0.5 * (a + b);
    let mut sum = 0.0;
    for (node, weight) in NODES.iter().zip(WEIGHTS.iter()) {
        let value = speed(mid + half * node)?;
        if !value.is_finite() {
            return Err(invalid(format!(
                "curve speed is not finite at parameter {}",
                mid + half * node
            )));
        }
        sum += weight * value;
    }
    Ok(sum * half)
}

/// Bisect `[a, b]` until each panel agrees with its two halves to its share
/// of `tolerance`; panels are processed from a stack in increasing
/// parameter, so the sum is deterministic.
fn adaptive<F>(
    speed: &F,
    a: Scalar,
    b: Scalar,
    tolerance: Scalar,
    panels: &mut usize,
) -> GeomResult<Scalar>
where
    F: Fn(Scalar) -> GeomResult<Scalar>,
{
    let width = b - a;
    let mut stack = vec![(a, b, rule(speed, a, b)?)];
    let mut total = 0.0;
    while let Some((lo, hi, whole)) = stack.pop() {
        *panels += 1;
        if *panels > MAX_PANELS {
            return Err(GeomError::BudgetExceeded {
                resource: "arc-length quadrature panels",
            });
        }
        let mid = 0.5 * (lo + hi);
        let left = rule(speed, lo, mid)?;
        let right = rule(speed, mid, hi)?;
        let share = tolerance * (hi - lo) / width;
        // A panel too narrow to bisect has nothing left to refine.
        if (left + right - whole).abs() <= share || !(mid > lo && mid < hi) {
            total += left + right;
        } else {
            stack.push((mid, hi, right));
            stack.push((lo, mid, left));
        }
    }
    Ok(total)
}

fn parameter_at<F>(family: &Family, speed: &F, from: Scalar, length: Scalar) -> GeomResult<Scalar>
where
    F: Fn(Scalar) -> GeomResult<Scalar>,
{
    check_finite(from, length, "a finite arc length")?;
    let (domain, breaks) = match family {
        Family::Proportional(rate) => return Ok(from + length / rate),
        Family::Numeric { domain, breaks } => (*domain, breaks),
    };
    if !inside(domain, from) {
        return Err(invalid(format!(
            "arc length is measured from parameter {from}, outside the curve's domain"
        )));
    }
    if length == 0.0 {
        return Ok(from);
    }
    let direction = length.signum();
    let target = length.abs();
    let tolerance = tolerance_for(target);
    // `g(x)`: arc length from `from` to `from + direction * x`, x >= 0.
    let travelled = |x: Scalar| -> GeomResult<Scalar> {
        let to = from + direction * x;
        integrate(speed, breaks, from.min(to), from.max(to))
    };
    // Bracket the root: up to the domain end, or by doubling when the
    // curve runs on.
    let (mut hi, mut at_hi) = match domain {
        Some((lo, hi)) => {
            let reach = if direction > 0.0 {
                hi - from
            } else {
                from - lo
            }
            .max(0.0);
            let available = travelled(reach)?;
            if available < target - tolerance {
                return Err(invalid(format!(
                    "arc length {target} exceeds the {available} available along the curve \
                     from parameter {from}"
                )));
            }
            if available <= target {
                return Ok(from + direction * reach);
            }
            (reach, available)
        }
        None => {
            let mut reach: Scalar = 1.0;
            let mut grown = 0;
            loop {
                let available = travelled(reach)?;
                if available >= target {
                    break (reach, available);
                }
                grown += 1;
                if grown > 64 {
                    return Err(invalid(format!(
                        "arc length {target} is not reached along the curve from parameter {from}"
                    )));
                }
                reach *= 2.0;
            }
        }
    };
    let mut lo = 0.0;
    let mut at_lo = 0.0;
    // Start from the linear interpolant of the bracket.
    let mut x = hi * target / at_hi;
    for _ in 0..MAX_ITERATIONS {
        if !(x > lo && x < hi) {
            x = 0.5 * (lo + hi);
        }
        let g = travelled(x)?;
        let residual = g - target;
        if residual.abs() <= tolerance {
            // One last Newton step costs no quadrature and squares the
            // residual; it is kept only if it stays inside the bracket.
            let rate = speed(from + direction * x)?;
            let polished = x - residual / rate;
            let x = if rate.is_finite() && rate > 0.0 && polished >= lo && polished <= hi {
                polished
            } else {
                x
            };
            return Ok(from + direction * x);
        }
        if residual < 0.0 {
            lo = x;
            at_lo = g;
        } else {
            hi = x;
            at_hi = g;
        }
        let midpoint = 0.5 * (lo + hi);
        if !(midpoint > lo && midpoint < hi) {
            // The bracket has collapsed to rounding: the parameter is known
            // as well as a float can hold it.
            return Ok(from + direction * x);
        }
        let rate = speed(from + direction * x)?;
        let newton = if rate.is_finite() && rate > 0.0 {
            x - residual / rate
        } else {
            Scalar::NAN
        };
        x = if newton > lo && newton < hi {
            newton
        } else {
            // Fall back to the secant of the bracket, then bisection.
            let secant = lo + (target - at_lo) * (hi - lo) / (at_hi - at_lo);
            if secant > lo && secant < hi {
                secant
            } else {
                midpoint
            }
        };
    }
    Err(GeomError::BudgetExceeded {
        resource: "arc-length inverse iterations",
    })
}
