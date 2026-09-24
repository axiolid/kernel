//! The exact tier: every finite double converts exactly, and arithmetic
//! matches 128-bit integer arithmetic wherever that is exact too.

use axiolid_exact::{Arith, Dyadic};
use axiolid_guarantees::Sign;
use num_bigint::BigInt;

fn d(value: f64) -> Dyadic {
    Dyadic::from_f64(value)
}

#[test]
fn doubles_convert_exactly_across_every_regime() {
    // Smallest subnormal: 1 * 2^-1074.
    let tiny = d(f64::from_bits(1));
    assert_eq!(
        (tiny.mantissa().clone(), tiny.exponent()),
        (BigInt::from(1), -1074)
    );
    // Largest double: (2^53 - 1) * 2^971.
    let max = d(f64::MAX);
    assert_eq!(max.mantissa(), &BigInt::from((1u64 << 53) - 1));
    assert_eq!(max.exponent(), 971);
    // 0.75 = 3 * 2^-2, normalised to an odd mantissa.
    let three_quarters = d(0.75);
    assert_eq!(
        (three_quarters.mantissa().clone(), three_quarters.exponent()),
        (BigInt::from(3), -2)
    );
    assert_eq!(d(-0.0), d(0.0), "both zeros are the same exact value");
    assert!(Dyadic::try_from_f64(f64::NAN).is_none());
    assert!(Dyadic::try_from_f64(f64::NEG_INFINITY).is_none());
}

#[test]
fn equal_values_have_equal_representations() {
    // Normalisation is what lets a zero difference compare equal to zero.
    assert_eq!(d(2.0).mul(&d(0.5)), d(1.0));
    assert_eq!(d(0.25).add(&d(0.75)), d(1.0));
    assert_eq!(d(3.0).sub(&d(3.0)), Dyadic::zero());
    assert_eq!(d(3.0).sub(&d(3.0)).sign(), Some(Sign::Zero));
}

#[test]
fn integer_arithmetic_matches_i128() {
    let mut state = 0x9E37_79B9_7F4A_7C15u64;
    let mut next = || {
        state ^= state << 13;
        state ^= state >> 7;
        state ^= state << 17;
        (state % 2_000_001) as i64 - 1_000_000
    };
    for _ in 0..10_000 {
        let (a, b, c) = (next(), next(), next());
        // a*b - c*c + a, all exact in i128.
        let expected =
            i128::from(a) * i128::from(b) - i128::from(c) * i128::from(c) + i128::from(a);
        let got = d(a as f64)
            .mul(&d(b as f64))
            .sub(&d(c as f64).square())
            .add(&d(a as f64));
        let exact = BigInt::from(expected);
        // Undo the normalisation to compare as integers.
        let value = if got.exponent() >= 0 {
            got.mantissa() << (got.exponent() as u64)
        } else {
            panic!("an integer result has a negative exponent")
        };
        assert_eq!(value, exact);
    }
}

#[test]
fn the_exact_tier_resolves_what_doubles_cannot() {
    // (1 + 2^-60) - 1 is 2^-60 exactly, but 1 + 2^-60 rounds to 1 in f64.
    // Built from exact parts, the dyadic difference keeps it.
    let tiny = d(2f64.powi(-60));
    let sum = d(1.0).add(&tiny);
    assert_eq!(sum.sub(&d(1.0)), tiny);
    assert_eq!(1.0 + 2f64.powi(-60) - 1.0, 0.0, "doubles lose it");
}
