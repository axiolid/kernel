//! `a + b*sqrt(c)` signs against an independent oracle.
//!
//! The code under test decides signs by squaring with case analysis. The
//! oracle does something unrelated: it evaluates the expression in 256-bit
//! fixed point using an integer square root, where the truncation error is
//! bounded by `|b| + |r|` units. For small integer coefficients a non-zero
//! value is astronomically larger than that bound, so the oracle's sign is
//! exact, and anything within the bound is an exact zero.

use axiolid_exact::{
    certify, filter, sign_root, sign_two_roots, Arith, Dyadic, Interval, Root2, SignExpr,
};
use axiolid_guarantees::Sign;
use num_bigint::BigInt;

const FIXED_BITS: u32 = 256;

fn scaled_sqrt(c: i64) -> BigInt {
    // floor(sqrt(c) * 2^256), from the integer square root of c * 2^512.
    (BigInt::from(c) << (2 * FIXED_BITS)).sqrt()
}

/// Sign of `p + q*sqrt(c) + r*sqrt(e)` by fixed-point evaluation.
fn oracle(p: i64, q: i64, c: i64, r: i64, e: i64) -> Sign {
    let approx = (BigInt::from(p) << FIXED_BITS)
        + BigInt::from(q) * scaled_sqrt(c)
        + BigInt::from(r) * scaled_sqrt(e);
    // Each floor loses less than one unit; scaled by |q| and |r|.
    let bound = BigInt::from(q.abs() + r.abs());
    if approx > bound {
        Sign::Positive
    } else if approx < -bound {
        Sign::Negative
    } else {
        Sign::Zero
    }
}

fn exact<T: Arith>(value: i64) -> T {
    T::from_f64(value as f64)
}

#[test]
fn one_root_matches_the_oracle_exhaustively() {
    let (mut zeros, mut cancellations) = (0, 0);
    for a in -7..=7 {
        for b in -7..=7 {
            for c in 0..=50 {
                let want = oracle(a, b, c, 0, 0);
                let got = sign_root::<Dyadic>(&exact(a), &exact(b), &exact(c));
                assert_eq!(got, Some(want), "sign({a} + {b}*sqrt({c}))");
                // A decided interval must agree; undecided is allowed.
                if let Some(fast) = sign_root::<Interval>(&exact(a), &exact(b), &exact(c)) {
                    assert_eq!(fast, want, "interval sign({a} + {b}*sqrt({c}))");
                }
                if want == Sign::Zero {
                    zeros += 1;
                    cancellations += usize::from(a != 0);
                }
            }
        }
    }
    // Counted independently: a + b*sqrt(c) = 0 needs a = 0 with b*sqrt(c) =
    // 0 (65 cases), or c a perfect square s^2 with a = -b*s (32 more).
    // Exact counts make a vacuous pass impossible.
    assert_eq!((zeros, cancellations), (97, 32));
}

#[test]
fn two_roots_match_the_oracle_exhaustively() {
    let (mut zeros, mut cancellations) = (0, 0);
    for p in -3..=3 {
        for q in -3..=3 {
            for r in -3..=3 {
                for c in 0..=9 {
                    for e in 0..=9 {
                        let want = oracle(p, q, c, r, e);
                        let got = sign_two_roots::<Dyadic>(
                            &exact(p),
                            &exact(q),
                            &exact(c),
                            &exact(r),
                            &exact(e),
                        );
                        assert_eq!(got, Some(want), "{p} + {q}√{c} + {r}√{e}");
                        if want == Sign::Zero {
                            zeros += 1;
                            let both_live = q != 0 && c != 0 && r != 0 && e != 0;
                            cancellations += usize::from(p != 0 || both_live);
                        }
                    }
                }
            }
        }
    }
    // Counted independently (80-digit decimal evaluation). The 484 include
    // cross-radicand cancellations such as sqrt(8) - 2*sqrt(2).
    assert_eq!((zeros, cancellations), (740, 484));
}

#[test]
fn a_negative_radicand_or_zero_denominator_is_undefined() {
    assert_eq!(sign_root::<Dyadic>(&exact(1), &exact(1), &exact(-1)), None);
    let bad = Root2 {
        a: exact::<Dyadic>(1),
        b: exact(1),
        c: exact(2),
        d: exact(0),
    };
    assert_eq!(bad.sign(), None);
}

#[test]
fn comparison_orders_values_with_different_radicands() {
    // (1 + sqrt(2)) / 2 = 1.207..., (3 - sqrt(3)) / 1 = 1.267...
    let x = Root2 {
        a: exact::<Dyadic>(1),
        b: exact(1),
        c: exact(2),
        d: exact(2),
    };
    let y = Root2 {
        a: exact::<Dyadic>(3),
        b: exact(-1),
        c: exact(3),
        d: exact(1),
    };
    assert_eq!(x.cmp_sign(&y), Some(Sign::Negative));
    assert_eq!(y.cmp_sign(&x), Some(Sign::Positive));
    assert_eq!(x.cmp_sign(&x), Some(Sign::Zero));
    // Negative denominators flip the value, not the logic: (-1 - sqrt 2)/-2 == x.
    let x_flipped = Root2 {
        a: exact::<Dyadic>(-1),
        b: exact(-1),
        c: exact(2),
        d: exact(-2),
    };
    assert_eq!(x.cmp_sign(&x_flipped), Some(Sign::Zero));
}

/// `p - q*sqrt(2)` for the sign question machinery.
struct PellGap {
    p: f64,
    q: f64,
}

impl SignExpr for PellGap {
    fn sign_in<T: Arith>(&self) -> Option<Sign> {
        sign_root(
            &T::from_f64(self.p),
            &T::from_f64(-self.q),
            &T::from_f64(2.0),
        )
    }
}

#[test]
fn pell_convergents_need_the_exact_tier_and_get_it_right() {
    // p/q -> sqrt(2) with p^2 - 2 q^2 = +-1, so p - q*sqrt(2) is about
    // 1 / (2 q sqrt 2): near 1e-17 for the largest pair doubles can hold.
    // f64 cannot resolve that against magnitudes near 1e16; the filter
    // must admit it cannot, and the exact tier must get it right.
    let (mut p, mut q) = (1i128, 1i128);
    let mut largest = Vec::new();
    while p < (1i128 << 53) {
        largest.push((p, q));
        (p, q) = (p + 2 * q, p + q);
    }
    let mut escalated = 0;
    for &(p, q) in largest.iter().rev().take(4) {
        let norm = p * p - 2 * q * q;
        assert!(norm == 1 || norm == -1, "Pell identity for {p}, {q}");
        let want = if norm > 0 {
            Sign::Positive
        } else {
            Sign::Negative
        };
        let expr = PellGap {
            p: p as f64,
            q: q as f64,
        };
        assert_eq!(certify(&expr), Ok(want), "sign of {p} - {q}*sqrt(2)");
        escalated += usize::from(!filter(&expr).is_certain());
    }
    assert!(escalated >= 2, "the largest convergents must escalate");
}
