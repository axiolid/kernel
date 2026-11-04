//! Certified orientation predicates.
//!
//! `orient2d` answers "does C lie left of, right of, or exactly on the directed
//! line AB". The answer is a *sign*, and a sign drives topology, so a plausible
//! answer is not good enough: this module returns [`Certified`] and escalates
//! to exact arithmetic rather than guessing when the filter is inconclusive.
//!
//! The strategy is the standard filtered cascade:
//!
//! 1. Evaluate the determinant in plain f64 and compute a forward error bound.
//! 2. If the magnitude exceeds the bound, the sign is proven; return it.
//! 3. Otherwise evaluate the determinant exactly (see below).
//!
//! Step 1 succeeds for almost all real inputs, so the exact path is rare.
//!
//! The exact path first tries two-term expansions, which are exact when the
//! coordinate differences carry no rounding error (grid-aligned collinear
//! points, the common degenerate case), and otherwise runs a fixed-size dyadic
//! evaluation over the original
//! coordinates, as `orient3d`'s fallback is (#404): it forms no difference and
//! no rounded product, so it is exact for every finite input. An earlier
//! expansion-based evaluation overflowed beyond about 1e154, lost error terms
//! to underflow, and could misjudge coordinates of very different magnitudes,
//! so it was removed rather than guarded. The filter itself declines an input
//! whose differences are outside [1e-90, 1e90], where its products may
//! overflow or underflow, and so hands it to the exact path.

use axiolid_core::Point2;
use axiolid_guarantees::{Certified, Precision, Sign};

use crate::expansion::{two_diff, two_product, two_sum};
use crate::orient2_dyadic::orient2d_exact_dyadic;

/// Machine epsilon for binary64: the gap between 1.0 and the next value.
const EPSILON: f64 = f64::EPSILON / 2.0;

/// Relative error bound for the 2x2 determinant filter.
///
/// The determinant costs two products and one subtraction; propagating the
/// standard `(1 + eps)` model over that expression yields `3*eps + O(eps^2)`.
/// The `16.0 * EPSILON * EPSILON` term absorbs the higher-order remainder, so
/// the bound is a true upper bound rather than a first-order approximation.
const ORIENT2D_ERROR_FACTOR: f64 = (3.0 + 16.0 * EPSILON) * EPSILON;

/// Orientation of `c` relative to the directed line `a` -> `b`.
///
/// Returns `Sign::Positive` when `a`, `b`, `c` turn counter-clockwise,
/// `Sign::Negative` for clockwise, and `Sign::Zero` when the three points
/// are exactly collinear.
///
/// The result is [`Certified::Certain`] for every finite input: this function
/// escalates to exact arithmetic internally, so it never returns an unproven
/// sign. A caller may therefore use it to drive a topology decision directly.
/// A non-finite coordinate returns [`Certified::Uncertain`], not a guess.
#[must_use]
pub fn orient2d(a: Point2, b: Point2, c: Point2) -> Certified {
    match orient2d_filter(a, b, c) {
        Certified::Certain { sign, .. } => Certified::exact_sign(sign),
        // The filter could not prove a sign, so pay for exact arithmetic. This
        // is the whole point of the cascade: correctness does not depend on the
        // fast path being lucky. `Certified` is non-exhaustive: a future
        // variant we do not understand must escalate, never be assumed decisive.
        _ => orient2d_exact_expansion(a, b, c)
            .or_else(|| orient2d_exact_dyadic(a, b, c))
            .map_or(
                Certified::Uncertain {
                    attempted: Precision::Exact,
                },
                Certified::exact_sign,
            ),
    }
}

/// Smallest nonzero rounding tail the expansion stage accepts. With the
/// differences in [1e-90, 1e90], a tail at or above this keeps every product
/// of two parts (at least 1e-240) and its own error term, some 2^-106 smaller
/// still, far above the subnormal range.
const MIN_TAIL: f64 = 1.0e-120;

/// Exact sign through expansions, or `None` when that is not exact.
///
/// Each coordinate difference is `high + tail` exactly (`two_diff`), so the
/// determinant is the exact sum of the products of those parts, as in the end
/// stage of Shewchuk's `orient2dadapt`. When every difference is in the
/// filter's range and every nonzero tail is at least [`MIN_TAIL`], each of
/// the eight part products is exactly `hi + lo`, so the sign is that of an
/// exact sum of up to sixteen doubles. This settles the common degenerate
/// case, exactly collinear points on a coordinate grid, without the bit-level
/// dyadic accumulator.
fn orient2d_exact_expansion(a: Point2, b: Point2, c: Point2) -> Option<Sign> {
    let [acx, bcy, acy, bcx] = differences(a, b, c)?;
    let tails = [
        two_diff(a.x, c.x).1,
        two_diff(b.y, c.y).1,
        two_diff(a.y, c.y).1,
        two_diff(b.x, c.x).1,
    ];
    if tails
        .iter()
        .any(|t| *t != 0.0 && !(MIN_TAIL..=1.0e90).contains(&t.abs()))
    {
        return None;
    }
    let [tacx, tbcy, tacy, tbcx] = tails;
    if tails == [0.0; 4] {
        // The usual case: no difference rounded, so only the main products.
        let (left_high, left_low) = two_product(acx, bcy);
        let (right_high, right_low) = two_product(acy, bcx);
        return sign_of_exact_sum(&[left_low, -right_low, left_high, -right_high]);
    }
    // Terms of `left - right`; a pair with a zero tail is skipped outright.
    let mut terms = [0.0_f64; 16];
    let mut count = 0;
    // Interleaved left and right, so the near-equal products cancel early and
    // the running expansion stays short.
    let pairs = [
        (acx, bcy, false),
        (acy, bcx, true),
        (acx, tbcy, false),
        (acy, tbcx, true),
        (tacx, bcy, false),
        (tacy, bcx, true),
        (tacx, tbcy, false),
        (tacy, tbcx, true),
    ];
    for (x, y, negate) in pairs {
        if x == 0.0 || y == 0.0 {
            continue;
        }
        let (high, low) = two_product(x, y);
        let sign = if negate { -1.0 } else { 1.0 };
        terms[count] = sign * high;
        terms[count + 1] = sign * low;
        count += 2;
    }
    sign_of_exact_sum(&terms[..count])
}

/// Sign of the exact sum of finite doubles (Shewchuk's grow-expansion).
fn sign_of_exact_sum(terms: &[f64]) -> Option<Sign> {
    let mut expansion = [0.0_f64; 16];
    let mut len = 0;
    for &term in terms.iter().filter(|t| **t != 0.0) {
        let mut carry = term;
        let mut kept = 0;
        for index in 0..len {
            let (sum, error) = two_sum(carry, expansion[index]);
            carry = sum;
            // Zero errors are dropped, which keeps the expansion short when
            // the terms cancel, as they do for collinear points.
            if error != 0.0 {
                expansion[kept] = error;
                kept += 1;
            }
        }
        expansion[kept] = carry;
        len = kept + 1;
    }
    // Components are non-overlapping and increasing, so the last nonzero one
    // carries the sign of the whole sum.
    let top = expansion[..len].iter().rev().find(|v| **v != 0.0);
    Some(match top {
        None => Sign::Zero,
        Some(v) if v.is_nan() => return None,
        Some(v) if *v > 0.0 => Sign::Positive,
        Some(_) => Sign::Negative,
    })
}

/// The four coordinate differences of the determinant, or `None` when one
/// overflowed or is so small that a product of two may underflow: such a
/// determinant cannot be bounded in f64 and is left to the exact path.
fn differences(a: Point2, b: Point2, c: Point2) -> Option<[f64; 4]> {
    let d = [a.x - c.x, b.y - c.y, a.y - c.y, b.x - c.x];
    d.iter()
        .all(|v| *v == 0.0 || (1.0e-90..=1.0e90).contains(&v.abs()))
        .then_some(d)
}

/// The fast filter alone, exposed so the escalation can be observed and tested.
///
/// Returns [`Certified::Uncertain`] when the floating-point determinant is too
/// close to zero for its own error bound to exclude the opposite sign.
#[must_use]
pub fn orient2d_filter(a: Point2, b: Point2, c: Point2) -> Certified {
    let Some([acx, bcy, acy, bcx]) = differences(a, b, c) else {
        return Certified::Uncertain {
            attempted: Precision::F64,
        };
    };
    let left = acx * bcy;
    let right = acy * bcx;
    let determinant = left - right;

    // The bound scales with the operand magnitudes: a determinant of 1e-9 is
    // decisive for millimetre coordinates and noise for national-grid ones.
    let magnitude = left.abs() + right.abs();
    let error_bound = ORIENT2D_ERROR_FACTOR * magnitude;

    Certified::from_filter(determinant, error_bound, Precision::F64)
}
