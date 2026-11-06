//! Fixed-size exact dyadic fallback for `orient2d`.
//!
//! The expansion path of `orient2d` multiplies coordinate differences, so it
//! overflows near 1e154 and loses error terms to underflow near 1e-154. This
//! path never forms a difference: it expands the determinant over the original
//! coordinates and accumulates each product of two binary64 significands at its
//! exact bit position, so it is exact for every finite input (#404).

use axiolid_core::Point2;
use axiolid_guarantees::Sign;

use crate::orient3_dyadic::{binary64_parts, OrientMagnitude};

/// Smallest exponent of a product of two binary64 values: `2 * -1074`.
const MIN_PAIR_EXPONENT: i32 = -2148;

/// Exact sign of the orientation determinant, `None` for a non-finite input.
///
/// Expanding `(ax - cx)(by - cy) - (ay - cy)(bx - cx)` cancels the `cx * cy`
/// terms and leaves six signed products of one `x` and one `y` coordinate.
#[must_use]
pub(super) fn orient2d_exact_dyadic(a: Point2, b: Point2, c: Point2) -> Option<Sign> {
    // (negated, x coordinate, y coordinate)
    let terms = [
        (false, a.x, b.y),
        (true, a.x, c.y),
        (true, c.x, b.y),
        (true, a.y, b.x),
        (false, a.y, c.x),
        (false, c.y, b.x),
    ];
    let mut positive = OrientMagnitude::zero();
    let mut negative = OrientMagnitude::zero();
    for (subtract, first, second) in terms {
        let (first_negative, first_mantissa, first_exponent) = binary64_parts(first)?;
        let (second_negative, second_mantissa, second_exponent) = binary64_parts(second)?;
        if first_mantissa == 0 || second_mantissa == 0 {
            continue;
        }
        let offset = usize::try_from(first_exponent + second_exponent - MIN_PAIR_EXPONENT).ok()?;
        let product = u128::from(first_mantissa) * u128::from(second_mantissa);
        let target = if subtract ^ first_negative ^ second_negative {
            &mut negative
        } else {
            &mut positive
        };
        if !target.add_shifted_product(&[product as u64, (product >> 64) as u64, 0], offset) {
            return None;
        }
    }
    Some(match positive.compare(&negative) {
        std::cmp::Ordering::Greater => Sign::Positive,
        std::cmp::Ordering::Less => Sign::Negative,
        std::cmp::Ordering::Equal => Sign::Zero,
    })
}
