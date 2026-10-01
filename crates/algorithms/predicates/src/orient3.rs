//! `orient3d`: which side of a plane a point lies on.
//!
//! Returns the sign of the 3x3 determinant
//!
//! ```text
//! | ax-dx  ay-dy  az-dz |
//! | bx-dx  by-dy  bz-dz |
//! | cx-dx  cy-dy  cz-dz |
//! ```
//!
//! Positive means `d` sees `a, b, c` counter-clockwise, i.e. `d` is below the
//! plane under the right-hand rule. Zero means the four points are exactly
//! coplanar -- the case that decides tetrahedralisation, convex hull facets,
//! and whether a boolean surface passes through a vertex.
//!
//! Same filtered cascade as `orient2d`: cheap f64 with an error bound, then
//! exact expansion arithmetic when the bound cannot exclude zero.

use axiolid_core::Point3;
use axiolid_guarantees::{Certified, Precision, Sign};

use crate::arithmetic::{
    expansion_product, expansion_sign, expansion_sum, grow_expansion, negate_expansion,
    scale_expansion,
};
use crate::expansion::{two_diff, two_product};
use crate::orient3_dyadic::orient3d_exact_dyadic;

/// Machine epsilon for binary64.
const EPSILON: f64 = f64::EPSILON / 2.0;

/// Relative error bound for the 3x3 determinant filter.
///
/// The determinant is a sum of three 2x2 cofactor products. Propagating the
/// `(1 + eps)` model over that expression gives `7*eps`; the second-order term
/// absorbs the remainder so this is a true upper bound.
///
/// The constant is deliberately conservative. A mutation probe lowering it to
/// `3*eps` does not fail the suite -- the filter stays sound at that value for
/// the inputs generated -- while `0.05*eps` is caught immediately by
/// `near_degenerate_cases_recover_a_definite_sign`. The margin between 3 and 7
/// buys nothing measurable in throughput (the escalation-rate gates show the
/// filter settles clean data either way) and costs nothing, so the derivation's
/// value is kept rather than the empirically minimal one: correctness here is
/// argued from the error model, not tuned against a test suite.
const ORIENT3D_ERROR_FACTOR: f64 = (7.0 + 56.0 * EPSILON) * EPSILON;

/// Orientation of `d` relative to the plane through `a`, `b`, `c`.
///
/// Returns an exact certified sign for finite binary64 coordinates. The usual
/// expansion path handles inputs whose intermediates remain representable;
/// other finite exponent ranges fall back to a fixed-size exact dyadic
/// accumulator. Non-finite input returns [`Certified::Uncertain`] rather than
/// guessing.
#[must_use]
pub fn orient3d(a: Point3, b: Point3, c: Point3, d: Point3) -> Certified {
    match orient3d_filter(a, b, c, d) {
        Certified::Certain { sign, .. } => Certified::exact_sign(sign),
        // Non-exhaustive enum: anything we do not recognise must escalate.
        _ => orient3d_exact(a, b, c, d)
            .or_else(|| orient3d_exact_dyadic(a, b, c, d))
            .map_or(
                Certified::Uncertain {
                    attempted: Precision::Exact,
                },
                Certified::exact_sign,
            ),
    }
}

/// The fast filter alone, exposed so escalation can be measured.
#[must_use]
pub fn orient3d_filter(a: Point3, b: Point3, c: Point3, d: Point3) -> Certified {
    let (adx, ady, adz) = (a.x - d.x, a.y - d.y, a.z - d.z);
    let (bdx, bdy, bdz) = (b.x - d.x, b.y - d.y, b.z - d.z);
    let (cdx, cdy, cdz) = (c.x - d.x, c.y - d.y, c.z - d.z);

    let differences = [adx, ady, adz, bdx, bdy, bdz, cdx, cdy, cdz];
    if !differences.iter().all(|value| {
        value.is_finite() && (*value == 0.0 || (value.abs() >= 1.0e-90 && value.abs() <= 1.0e90))
    }) {
        return Certified::Uncertain {
            attempted: Precision::F64,
        };
    }

    let bdxcdy = bdx * cdy;
    let cdxbdy = cdx * bdy;
    let cdxady = cdx * ady;
    let adxcdy = adx * cdy;
    let adxbdy = adx * bdy;
    let bdxady = bdx * ady;

    let determinant = adz * (bdxcdy - cdxbdy) + bdz * (cdxady - adxcdy) + cdz * (adxbdy - bdxady);

    // The bound tracks operand magnitudes, so it scales with the model's units
    // instead of assuming a coordinate range.
    let permanent = (bdxcdy.abs() + cdxbdy.abs()) * adz.abs()
        + (cdxady.abs() + adxcdy.abs()) * bdz.abs()
        + (adxbdy.abs() + bdxady.abs()) * cdz.abs();

    Certified::from_filter(
        determinant,
        ORIENT3D_ERROR_FACTOR * permanent,
        Precision::F64,
    )
}

/// Exact sign of the 3x3 determinant.
///
/// Every coordinate difference includes the tail recovered by [`two_diff`].
/// Cofactors and their remaining products operate on complete expansions, so
/// the result is the determinant of the input binary64 coordinates rather than
/// of rounded coordinate differences.
#[must_use]
fn orient3d_exact(a: Point3, b: Point3, c: Point3, d: Point3) -> Option<Sign> {
    let (adx, adx_error) = two_diff(a.x, d.x);
    let (ady, ady_error) = two_diff(a.y, d.y);
    let (adz, adz_error) = two_diff(a.z, d.z);
    let (bdx, bdx_error) = two_diff(b.x, d.x);
    let (bdy, bdy_error) = two_diff(b.y, d.y);
    let (bdz, bdz_error) = two_diff(b.z, d.z);
    let (cdx, cdx_error) = two_diff(c.x, d.x);
    let (cdy, cdy_error) = two_diff(c.y, d.y);
    let (cdz, cdz_error) = two_diff(c.z, d.z);

    let differences = [[adx, ady, adz], [bdx, bdy, bdz], [cdx, cdy, cdz]];
    let errors = [
        [adx_error, ady_error, adz_error],
        [bdx_error, bdy_error, bdz_error],
        [cdx_error, cdy_error, cdz_error],
    ];
    if !exact_components_are_representable(&differences, &errors) {
        return None;
    }

    if errors.iter().flatten().all(|error| *error == 0.0) {
        return Some(orient3d_exact_differences(
            differences[0],
            differences[1],
            differences[2],
        ));
    }

    let [[adx, ady, adz], [bdx, bdy, bdz], [cdx, cdy, cdz]] = differences;
    let [[adx_error, ady_error, adz_error], [bdx_error, bdy_error, bdz_error], [cdx_error, cdy_error, cdz_error]] =
        errors;

    let adx = difference_expansion(adx, adx_error);
    let ady = difference_expansion(ady, ady_error);
    let adz = difference_expansion(adz, adz_error);
    let bdx = difference_expansion(bdx, bdx_error);
    let bdy = difference_expansion(bdy, bdy_error);
    let bdz = difference_expansion(bdz, bdz_error);
    let cdx = difference_expansion(cdx, cdx_error);
    let cdy = difference_expansion(cdy, cdy_error);
    let cdz = difference_expansion(cdz, cdz_error);

    let bc = orient3d_expansion_cofactor(&bdx, &cdy, &cdx, &bdy);
    let ca = orient3d_expansion_cofactor(&cdx, &ady, &adx, &cdy);
    let ab = orient3d_expansion_cofactor(&adx, &bdy, &bdx, &ady);

    let total = expansion_sum(
        &expansion_sum(&expansion_product(&bc, &adz), &expansion_product(&ca, &bdz)),
        &expansion_product(&ab, &cdz),
    );
    Some(expansion_sign(&total))
}

/// Whether the ordinary expansion evaluator can represent every intermediate.
///
/// Components above exponent 300 could overflow a three-factor term. The
/// monomial check separately proves that each two- and three-factor product
/// stays above binary64's underflow floor.
#[must_use]
fn exact_components_are_representable(differences: &[[f64; 3]; 3], errors: &[[f64; 3]; 3]) -> bool {
    differences
        .iter()
        .flatten()
        .chain(errors.iter().flatten())
        .all(|value| value.is_finite() && (*value == 0.0 || highest_bit_exponent(*value) <= 300))
        && determinant_terms_are_representable(differences, errors)
}

#[must_use]
fn determinant_terms_are_representable(
    differences: &[[f64; 3]; 3],
    errors: &[[f64; 3]; 3],
) -> bool {
    let mut least_bits = [[None; 3]; 3];
    for row in 0..3 {
        for column in 0..3 {
            least_bits[row][column] = [differences[row][column], errors[row][column]]
                .into_iter()
                .filter(|value| *value != 0.0)
                .map(least_significant_bit_exponent)
                .min();
        }
    }

    [
        [(1, 0), (2, 1), (0, 2)],
        [(2, 0), (1, 1), (0, 2)],
        [(2, 0), (0, 1), (1, 2)],
        [(0, 0), (2, 1), (1, 2)],
        [(0, 0), (1, 1), (2, 2)],
        [(1, 0), (0, 1), (2, 2)],
    ]
    .into_iter()
    .all(|term| {
        let exponents = term.map(|(row, column)| least_bits[row][column]);
        let [Some(first), Some(second), Some(third)] = exponents else {
            return true;
        };
        first + second >= -1074 && first + second + third >= -1074
    })
}

#[must_use]
pub(crate) fn highest_bit_exponent(value: f64) -> i32 {
    let bits = value.abs().to_bits();
    let encoded_exponent = ((bits >> 52) & 0x7ff) as i32;
    if encoded_exponent == 0 {
        let fraction = bits & ((1_u64 << 52) - 1);
        -1074 + (63 - fraction.leading_zeros() as i32)
    } else {
        encoded_exponent - 1023
    }
}

#[must_use]
pub(crate) fn least_significant_bit_exponent(value: f64) -> i32 {
    let bits = value.abs().to_bits();
    let encoded_exponent = ((bits >> 52) & 0x7ff) as i32;
    let fraction = bits & ((1_u64 << 52) - 1);
    if encoded_exponent == 0 {
        -1074 + fraction.trailing_zeros() as i32
    } else {
        let significand = (1_u64 << 52) | fraction;
        encoded_exponent - 1023 - 52 + significand.trailing_zeros() as i32
    }
}

/// Evaluate the determinant when all coordinate differences are already exact.
#[must_use]
fn orient3d_exact_differences(a: [f64; 3], b: [f64; 3], c: [f64; 3]) -> Sign {
    let bc = orient3d_cofactor(b[0], c[1], c[0], b[1]);
    let ca = orient3d_cofactor(c[0], a[1], a[0], c[1]);
    let ab = orient3d_cofactor(a[0], b[1], b[0], a[1]);
    let total = expansion_sum(
        &expansion_sum(&scale_expansion(&bc, a[2]), &scale_expansion(&ca, b[2])),
        &scale_expansion(&ab, c[2]),
    );
    expansion_sign(&total)
}

#[must_use]
fn difference_expansion(difference: f64, error: f64) -> Vec<f64> {
    let mut expansion = Vec::new();
    if error != 0.0 {
        expansion.push(error);
    }
    if difference != 0.0 || expansion.is_empty() {
        expansion.push(difference);
    }
    expansion
}

#[must_use]
fn orient3d_expansion_cofactor(p: &[f64], q: &[f64], r: &[f64], s: &[f64]) -> Vec<f64> {
    expansion_sum(
        &expansion_product(p, q),
        &negate_expansion(&expansion_product(r, s)),
    )
}

/// Exact `p*q - r*s` as an expansion.
///
/// Shared with the Delaunay predicates, which need the same 2x2 minors.
#[must_use]
pub(crate) fn orient3d_cofactor(p: f64, q: f64, r: f64, s: f64) -> Vec<f64> {
    let (pq, pq_err) = two_product(p, q);
    let (rs, rs_err) = two_product(r, s);
    // Subtract by adding the negation; the four terms are combined by the
    // carry-propagating grow so the result stays a valid expansion.
    let e = grow_expansion(&[pq_err], -rs_err);
    let e = grow_expansion(&e, pq);
    grow_expansion(&e, -rs)
}
