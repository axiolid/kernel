//! Exact real-root isolation, checked against polynomials built from
//! known roots. Every expected answer is computed without floating point:
//! dyadic roots compare exactly, and sqrt(k) against a dyadic d compares
//! by squaring (sign(sqrt(k) - d) = sign(k - d^2) for d >= 0).

use axiolid_exact::{Arith, Dyadic, IntPoly, RealRoot};
use axiolid_guarantees::Sign;
use num_bigint::BigInt;

fn big(v: i64) -> BigInt {
    BigInt::from(v)
}

fn mul(a: &[BigInt], b: &[BigInt]) -> Vec<BigInt> {
    let mut out = vec![big(0); a.len() + b.len() - 1];
    for (i, x) in a.iter().enumerate() {
        for (j, y) in b.iter().enumerate() {
            out[i + j] += x * y;
        }
    }
    out
}

/// A known root: dyadic m * 2^-e, or +-sqrt(k).
#[derive(Clone, Debug)]
enum Known {
    Dyadic(i64, i64),
    Sqrt(i64, bool),
}

impl Known {
    fn value(&self) -> f64 {
        match *self {
            Known::Dyadic(m, e) => m as f64 / (1i64 << e) as f64,
            Known::Sqrt(k, neg) => {
                let r = (k as f64).sqrt();
                if neg {
                    -r
                } else {
                    r
                }
            }
        }
    }

    /// Exact sign of `self - d`.
    fn cmp_dyadic(&self, d: &Dyadic) -> Sign {
        match *self {
            Known::Dyadic(m, e) => Dyadic::from_parts(big(m), -e).sub(d).sign().unwrap(),
            Known::Sqrt(k, neg) => {
                let dneg = d.sign() == Some(Sign::Negative);
                let dd = if dneg { d.neg() } else { d.clone() };
                // |root| vs |d|
                let mag = Dyadic::from_f64(k as f64).sub(&dd.square()).sign().unwrap();
                match (neg, dneg) {
                    (false, false) => mag,
                    (true, true) => mag.flip(),
                    (false, true) => Sign::Positive,
                    (true, false) => Sign::Negative,
                }
            }
        }
    }
}

fn factor(root: &Known) -> Vec<BigInt> {
    match *root {
        // 2^e x - m
        Known::Dyadic(m, e) => vec![big(-m), big(1) << e as u64],
        // x^2 - k (both signs share it)
        Known::Sqrt(k, _) => vec![big(-k), big(0), big(1)],
    }
}

/// Exact check that `found` and the distinct `known` roots agree one to one,
/// in order.
fn check(roots: &[Known], context: &str) {
    let mut coeffs = vec![big(1)];
    for r in roots {
        coeffs = mul(&coeffs, &factor(r));
    }
    let poly = IntPoly::new(coeffs);
    // Distinct known roots, sorted by value (values here are far apart
    // enough, or equal, that f64 sorting of the KNOWN list is safe).
    let mut distinct: Vec<Known> = Vec::new();
    for r in roots {
        let rs: Vec<Known> = match *r {
            Known::Sqrt(k, _) => vec![Known::Sqrt(k, false), Known::Sqrt(k, true)],
            ref d => vec![d.clone()],
        };
        for x in rs {
            if !distinct
                .iter()
                .any(|y| (y.value() - x.value()).abs() < 1e-12)
            {
                distinct.push(x);
            }
        }
    }
    distinct.sort_by(|a, b| a.value().partial_cmp(&b.value()).unwrap());
    let found = poly.real_roots();
    assert_eq!(found.len(), distinct.len(), "{context}: root count");
    for (f, k) in found.iter().zip(&distinct) {
        let (lo, hi) = f.bounds();
        // The known root lies in the isolating interval (exactly at it if
        // the interval is a point).
        if f.is_exact() {
            assert_eq!(k.cmp_dyadic(lo), Sign::Zero, "{context}: exact root {k:?}");
        } else {
            assert_eq!(
                k.cmp_dyadic(lo),
                Sign::Positive,
                "{context}: {k:?} above lo"
            );
            assert_eq!(
                k.cmp_dyadic(hi),
                Sign::Negative,
                "{context}: {k:?} below hi"
            );
            assert_ne!(
                f.poly().sign_at(lo),
                Sign::Zero,
                "{context}: open lo is no root"
            );
            assert_ne!(
                f.poly().sign_at(hi),
                Sign::Zero,
                "{context}: open hi is no root"
            );
        }
    }
    for pair in found.windows(2) {
        assert_eq!(
            pair[0].cmp_root(&pair[1]),
            Sign::Negative,
            "{context}: order"
        );
    }
}

#[test]
fn simple_and_repeated_roots() {
    use Known::*;
    check(&[Dyadic(1, 0), Dyadic(-3, 1), Dyadic(5, 2)], "three dyadic");
    // Repeated roots: reported once, as a square-free polynomial would.
    check(
        &[Dyadic(1, 0), Dyadic(1, 0), Dyadic(1, 0), Dyadic(-2, 0)],
        "triple root",
    );
    check(&[Sqrt(2, false), Sqrt(2, false)], "(x^2-2)^2");
    // Zero, and roots exactly on the first bisection points.
    check(
        &[Dyadic(0, 0), Dyadic(4, 0), Dyadic(-4, 0), Dyadic(2, 0)],
        "split points",
    );
    check(
        &[Sqrt(2, false), Sqrt(3, false), Dyadic(3, 1), Dyadic(-7, 2)],
        "mixed",
    );
}

#[test]
fn roots_closer_than_f64_can_separate() {
    // sqrt(2) against the double nearest it, as an exact dyadic root.
    let near = 2f64.sqrt();
    let d = Dyadic::from_f64(near);
    let linear = IntPoly::from_dyadic(&[d.neg(), Dyadic::from_f64(1.0)]);
    let quad = IntPoly::new(vec![big(-2), big(0), big(1)]);
    let root2 = quad.real_roots().pop().unwrap();
    let dbl = linear.real_roots().pop().unwrap();
    // The double 1.4142135623730951 is above sqrt(2).
    assert_eq!(root2.cmp_root(&dbl), Sign::Negative);
    assert_eq!(dbl.cmp_root(&root2), Sign::Positive);
    assert_eq!(root2.cmp_dyadic(&d), Sign::Negative);
}

#[test]
fn equal_roots_of_different_polynomials_compare_equal() {
    // sqrt(2) as a root of three different polynomials.
    let a = IntPoly::new(mul(&[big(-2), big(0), big(1)], &[big(-1), big(1)]));
    let b = IntPoly::new(mul(&[big(-2), big(0), big(1)], &[big(3), big(1)]));
    let c = IntPoly::new(vec![big(-4), big(0), big(0), big(0), big(1)]); // x^4 - 4
    let pick = |p: &IntPoly| {
        p.real_roots()
            .into_iter()
            .find(|r| {
                r.cmp_dyadic(&Dyadic::from_f64(1.4)) == Sign::Positive
                    && r.cmp_dyadic(&Dyadic::from_f64(1.5)) == Sign::Negative
            })
            .expect("sqrt 2")
    };
    let (ra, rb, rc) = (pick(&a), pick(&b), pick(&c));
    assert_eq!(ra.cmp_root(&rb), Sign::Zero);
    assert_eq!(rb.cmp_root(&rc), Sign::Zero);
    // And not equal to sqrt(3), which shares no factor.
    let s3 = IntPoly::new(vec![big(-3), big(0), big(1)])
        .real_roots()
        .pop()
        .unwrap();
    assert_eq!(ra.cmp_root(&s3), Sign::Negative);
}

#[test]
fn mignotte_close_roots_are_separated() {
    // x^4 - 2*(1000x - 1)^2 has two real roots within about 1e-9 of 1e-3.
    let sq = mul(&[big(-1), big(1000)], &[big(-1), big(1000)]);
    let mut coeffs = vec![big(0); 5];
    for (i, c) in sq.iter().enumerate() {
        coeffs[i] -= c * big(2);
    }
    coeffs[4] += big(1);
    let poly = IntPoly::new(coeffs);
    let roots = poly.real_roots();
    let near: Vec<&RealRoot> = roots
        .iter()
        .filter(|r| (r.approx() - 1e-3).abs() < 1e-6)
        .collect();
    assert_eq!(near.len(), 2, "both close roots isolated");
    assert_eq!(near[0].cmp_root(near[1]), Sign::Negative);
    for r in &roots {
        assert!(!r.is_exact() || r.poly().sign_at(r.bounds().0) == Sign::Zero);
    }
}

#[test]
fn random_products_of_known_factors() {
    let mut s: u64 = 0xC0FF_EE00_1234_5678;
    let mut next = || {
        s ^= s << 13;
        s ^= s >> 7;
        s ^= s << 17;
        s
    };
    let squares = [4i64, 9, 16, 25];
    let nonsquares = [2i64, 3, 5, 6, 7, 10, 11];
    for case in 0..300 {
        let n = 1 + (next() % 5) as usize;
        let mut roots = Vec::new();
        for _ in 0..n {
            if next() % 3 == 0 {
                let k = nonsquares[(next() % nonsquares.len() as u64) as usize];
                roots.push(Known::Sqrt(k, false));
            } else {
                // Dyadic roots in [-8, 8], often on bisection points; a
                // perfect square's root duplicates an integer root.
                let e = (next() % 4) as i64;
                let m = (next() % 33) as i64 - 16;
                roots.push(Known::Dyadic(m << (e.min(1)), e));
                if next() % 7 == 0 {
                    let k = squares[(next() % 4) as usize];
                    let r = (k as f64).sqrt() as i64;
                    roots.push(Known::Dyadic(r, 0));
                }
            }
        }
        check(&roots, &format!("case {case}"));
    }
}

#[test]
fn approximations_are_close() {
    let quad = IntPoly::new(vec![big(-2), big(0), big(1)]);
    let roots = quad.real_roots();
    assert_eq!(roots.len(), 2);
    assert!((roots[1].approx() - 2f64.sqrt()).abs() < 1e-15);
    assert!((roots[0].approx() + 2f64.sqrt()).abs() < 1e-15);
    assert!(IntPoly::new(vec![big(1), big(0), big(1)])
        .real_roots()
        .is_empty());
}
