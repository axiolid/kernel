//! Fixed-size exact dyadic fallback for `orient3d`, and the accumulator that
//! the `orient2d`, `incircle` and `insphere` fallbacks share.
//!
//! This path is used only when binary64 expansion arithmetic cannot normalize
//! finite inputs without overflow or underflow. It performs bounded stack-only
//! work over the exact binary representation of the original coordinates.

use axiolid_core::Point3;
use axiolid_guarantees::Sign;

const DYADIC_LIMBS: usize = 100;
const MIN_TRIPLE_EXPONENT: i32 = -3222;

/// The accumulator the `orient2d` and `orient3d` fallbacks use.
pub(super) type OrientMagnitude = ExactMagnitude<DYADIC_LIMBS>;

/// A non-negative integer of `LIMBS` 64-bit words, little-endian.
pub(super) struct ExactMagnitude<const LIMBS: usize> {
    limbs: [u64; LIMBS],
}

impl<const LIMBS: usize> ExactMagnitude<LIMBS> {
    pub(super) fn zero() -> Self {
        Self { limbs: [0; LIMBS] }
    }

    /// Adds `product << offset` (little-endian words) whole words at a time,
    /// or returns `false` when a set bit would land beyond the accumulator.
    pub(super) fn add_shifted_product(&mut self, product: &[u64], offset: usize) -> bool {
        let base = offset / 64;
        let shift = offset % 64;
        for (word_index, word) in product.iter().copied().enumerate() {
            if !self.add_limb(base + word_index, word << shift) {
                return false;
            }
            if shift != 0 && !self.add_limb(base + word_index + 1, word >> (64 - shift)) {
                return false;
            }
        }
        true
    }

    /// Adds `value` at limb `index` and propagates the carry upward.
    fn add_limb(&mut self, mut index: usize, mut value: u64) -> bool {
        while value != 0 {
            if index >= LIMBS {
                return false;
            }
            let (sum, carry) = self.limbs[index].overflowing_add(value);
            self.limbs[index] = sum;
            value = u64::from(carry);
            index += 1;
        }
        true
    }

    pub(super) fn compare(&self, other: &Self) -> std::cmp::Ordering {
        for index in (0..LIMBS).rev() {
            match self.limbs[index].cmp(&other.limbs[index]) {
                std::cmp::Ordering::Equal => {}
                ordering => return ordering,
            }
        }
        std::cmp::Ordering::Equal
    }
}

/// Exact sign fallback over the original coordinates, without subtraction.
///
/// Expanding each `(coordinate - d)` factor yields 48 signed products. Every
/// binary64 value is an integer significand times a power of two, so those
/// products accumulate exactly into fixed positive and negative magnitudes.
#[must_use]
pub(super) fn orient3d_exact_dyadic(a: Point3, b: Point3, c: Point3, d: Point3) -> Option<Sign> {
    let points = [
        [a.x, a.y, a.z],
        [b.x, b.y, b.z],
        [c.x, c.y, c.z],
        [d.x, d.y, d.z],
    ];
    let monomials = [
        (false, [(1, 0), (2, 1), (0, 2)]),
        (true, [(2, 0), (1, 1), (0, 2)]),
        (false, [(2, 0), (0, 1), (1, 2)]),
        (true, [(0, 0), (2, 1), (1, 2)]),
        (false, [(0, 0), (1, 1), (2, 2)]),
        (true, [(1, 0), (0, 1), (2, 2)]),
    ];
    let mut positive = OrientMagnitude::zero();
    let mut negative = OrientMagnitude::zero();

    for (subtract, factors) in monomials {
        for choices in 0_u8..8 {
            let mut mantissas = [0_u64; 3];
            let mut exponent = 0_i32;
            let mut is_negative = subtract;
            for index in 0..3 {
                let choose_d = choices & (1 << index) != 0;
                let (point, coordinate) = factors[index];
                let value = points[if choose_d { 3 } else { point }][coordinate];
                let (value_negative, mantissa, value_exponent) = binary64_parts(value)?;
                if mantissa == 0 {
                    mantissas = [0; 3];
                    break;
                }
                mantissas[index] = mantissa;
                exponent += value_exponent;
                is_negative ^= value_negative ^ choose_d;
            }
            if mantissas.contains(&0) {
                continue;
            }

            let offset = usize::try_from(exponent - MIN_TRIPLE_EXPONENT).ok()?;
            let product = multiply_three_significands(mantissas);
            let target = if is_negative {
                &mut negative
            } else {
                &mut positive
            };
            if !target.add_shifted_product(&product, offset) {
                return None;
            }
        }
    }

    Some(match positive.compare(&negative) {
        std::cmp::Ordering::Greater => Sign::Positive,
        std::cmp::Ordering::Less => Sign::Negative,
        std::cmp::Ordering::Equal => Sign::Zero,
    })
}

#[must_use]
pub(super) fn binary64_parts(value: f64) -> Option<(bool, u64, i32)> {
    if !value.is_finite() {
        return None;
    }
    let bits = value.to_bits();
    let negative = bits >> 63 != 0;
    let encoded_exponent = ((bits >> 52) & 0x7ff) as i32;
    let fraction = bits & ((1_u64 << 52) - 1);
    if encoded_exponent == 0 {
        Some((negative, fraction, -1074))
    } else {
        Some((
            negative,
            (1_u64 << 52) | fraction,
            encoded_exponent - 1023 - 52,
        ))
    }
}

#[must_use]
fn multiply_three_significands([a, b, c]: [u64; 3]) -> [u64; 3] {
    let ab = u128::from(a) * u128::from(b);
    let low = u128::from(ab as u64) * u128::from(c);
    let middle = u128::from((ab >> 64) as u64) * u128::from(c) + (low >> 64);
    [low as u64, middle as u64, (middle >> 64) as u64]
}
