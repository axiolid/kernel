//! Fixed-size exact dyadic fallbacks for `incircle` and `insphere`.
//!
//! The expansion paths multiply coordinate differences, a degree-4 (`incircle`)
//! or degree-5 (`insphere`) product, so they overflow beyond about 1e60 and
//! lose error terms to underflow below about 1e-60 (#430). These paths never
//! form a difference: they expand the lifted determinant over the original
//! coordinates and add each product of binary64 significands at its exact bit
//! position, so they are exact for every finite input.
//!
//! Both determinants are the determinant of the matrix whose rows are
//! `(p - o)` and its squared length, for `p` over the points and `o` the last
//! point. Expanding along the permutations, over which row carries the lift
//! and over the axis the lift squares, leaves monomials of `D + 2` differences
//! (`incircle` 4, `insphere` 5); expanding each difference into the point and
//! the reference point leaves signed products of `D + 2` coordinates.

use axiolid_core::{Point2, Point3};
use axiolid_guarantees::Sign;

use crate::orient3_dyadic::{binary64_parts, ExactMagnitude};

/// Smallest exponent of a binary64 significand taken as an integer.
const MIN_FACTOR_EXPONENT: i32 = -1074;

/// Words of the `incircle` accumulator. A product of four coordinates lands
/// below bit `4 * (971 + 53 + 1074) = 8392`, and 192 of them add under 8 more
/// bits, so 133 words (8512 bits) cannot overflow.
const INCIRCLE_LIMBS: usize = 133;

/// Words of the `insphere` accumulator: five coordinates reach bit
/// `5 * 2098 = 10490`, 2304 products add under 12 bits, and 166 words hold
/// 10624.
const INSPHERE_LIMBS: usize = 166;

/// Exact sign of the lifted `incircle` determinant, `None` for a non-finite
/// input.
#[must_use]
pub(super) fn incircle_exact_dyadic(a: Point2, b: Point2, c: Point2, d: Point2) -> Option<Sign> {
    lifted_determinant_sign::<2, 4, INCIRCLE_LIMBS>(
        &[[a.x, a.y], [b.x, b.y], [c.x, c.y]],
        [d.x, d.y],
    )
}

/// Exact sign of the lifted `insphere` determinant, `None` for a non-finite
/// input.
#[must_use]
pub(super) fn insphere_exact_dyadic(
    a: Point3,
    b: Point3,
    c: Point3,
    d: Point3,
    e: Point3,
) -> Option<Sign> {
    lifted_determinant_sign::<3, 5, INSPHERE_LIMBS>(
        &[
            [a.x, a.y, a.z],
            [b.x, b.y, b.z],
            [c.x, c.y, c.z],
            [d.x, d.y, d.z],
        ],
        [e.x, e.y, e.z],
    )
}

/// Sign of the `(D + 1) x (D + 1)` determinant with rows
/// `(p - origin, |p - origin|^2)`, over `D + 1` `rows`.
///
/// `FACTORS` is `D + 2`, the number of coordinates in one monomial.
fn lifted_determinant_sign<const D: usize, const FACTORS: usize, const LIMBS: usize>(
    rows: &[[f64; D]],
    origin: [f64; D],
) -> Option<Sign> {
    debug_assert!(rows.len() == D + 1 && FACTORS == D + 2);
    let mut positive = ExactMagnitude::<LIMBS>::zero();
    let mut negative = ExactMagnitude::<LIMBS>::zero();

    // `columns[row]` is the column of the matrix that row contributes to this
    // permutation; column `D` is the lift.
    let mut columns: [usize; 4] = [0, 1, 2, 3];
    let columns = &mut columns[..=D];
    loop {
        let inversions = (0..columns.len())
            .flat_map(|i| (i + 1..columns.len()).map(move |j| (i, j)))
            .filter(|&(i, j)| columns[i] > columns[j])
            .count();
        let lift_row = columns.iter().position(|&column| column == D)?;
        for lift_axis in 0..D {
            // (row, axis) of every coordinate factor of the monomial.
            let mut factors = [(0_usize, 0_usize); FACTORS];
            let mut count = 0;
            for (row, &column) in columns.iter().enumerate() {
                if row != lift_row {
                    factors[count] = (row, column);
                    count += 1;
                }
            }
            factors[count] = (lift_row, lift_axis);
            factors[count + 1] = (lift_row, lift_axis);

            // The terms with every factor taken from the origin are those of a
            // matrix with equal rows, which sum to zero, so they are left out.
            for choices in 0_u32..(1 << FACTORS) - 1 {
                let mut product = [0_u64; FACTORS];
                product[0] = 1;
                let mut exponent = 0_i32;
                let mut is_negative = inversions % 2 == 1;
                let mut vanishes = false;
                for (index, &(row, axis)) in factors.iter().enumerate() {
                    // Each difference `p - origin` is the point or minus the origin.
                    let from_origin = choices & (1 << index) != 0;
                    let value = if from_origin {
                        origin[axis]
                    } else {
                        rows[row][axis]
                    };
                    let (value_negative, mantissa, value_exponent) = binary64_parts(value)?;
                    if mantissa == 0 {
                        vanishes = true;
                        continue;
                    }
                    multiply_in_place(&mut product, mantissa);
                    exponent += value_exponent;
                    is_negative ^= value_negative ^ from_origin;
                }
                if vanishes {
                    continue;
                }
                let offset =
                    usize::try_from(exponent - MIN_FACTOR_EXPONENT * FACTORS as i32).ok()?;
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
        if !next_permutation(columns) {
            break;
        }
    }

    Some(match positive.compare(&negative) {
        std::cmp::Ordering::Greater => Sign::Positive,
        std::cmp::Ordering::Less => Sign::Negative,
        std::cmp::Ordering::Equal => Sign::Zero,
    })
}

/// `product *= factor` over little-endian words; the caller sizes `product`
/// so the result fits.
fn multiply_in_place(product: &mut [u64], factor: u64) {
    let mut carry = 0_u128;
    for word in product.iter_mut() {
        let wide = u128::from(*word) * u128::from(factor) + carry;
        *word = wide as u64;
        carry = wide >> 64;
    }
    debug_assert!(carry == 0);
}

/// Advances `items` to the next permutation in lexicographic order, or
/// returns `false` after the last one.
fn next_permutation(items: &mut [usize]) -> bool {
    let Some(pivot) = (1..items.len()).rev().find(|&i| items[i - 1] < items[i]) else {
        return false;
    };
    let pivot = pivot - 1;
    let swap = (pivot + 1..items.len())
        .rev()
        .find(|&i| items[i] > items[pivot])
        .expect("a larger item follows the pivot");
    items.swap(pivot, swap);
    items[pivot + 1..].reverse();
    true
}
