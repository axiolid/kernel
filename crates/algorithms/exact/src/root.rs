//! Numbers of the form `(a + b*sqrt(c)) / d`.
//!
//! A line meets a circle at such a parameter; comparing two such hits, or
//! asking which side of a line one lies on, is a sign question about them.
//! Signs are decided by squaring with case analysis (the approach of CGAL's
//! `Root_of_2`), never by evaluating a root, so in exact arithmetic the
//! answer is exact, including an exact zero.
//!
//! Degree note: `sign_root` squares once (degree 2 in `a`, `b`, `c`);
//! `sign_two_roots` squares twice. Mantissa length, and so exact-tier cost,
//! grows with that degree. Both run the interval filter first.

use axiolid_guarantees::Sign;

use crate::arith::{sign_product, Arith};

/// The value `(a + b*sqrt(c)) / d`, with `c >= 0` and `d != 0`.
#[derive(Debug, Clone, PartialEq)]
pub struct Root2<T> {
    /// Rational part of the numerator.
    pub a: T,
    /// Coefficient of the root.
    pub b: T,
    /// Radicand; must not be negative.
    pub c: T,
    /// Denominator; must not be zero.
    pub d: T,
}

impl<T: Arith> Root2<T> {
    /// Sign of the value, or `None` when `T` cannot decide (or the
    /// preconditions fail: negative radicand, zero denominator).
    #[must_use]
    pub fn sign(&self) -> Option<Sign> {
        let denominator = self.d.sign()?;
        if denominator == Sign::Zero {
            return None;
        }
        Some(sign_product(
            sign_root(&self.a, &self.b, &self.c)?,
            denominator,
        ))
    }

    /// Sign of `self - other`: which of the two values is larger.
    ///
    /// Works for different radicands. With `x = (a1 + b1*sqrt(c1)) / d1`
    /// and `y` likewise, `x - y` has the sign of `d1 * d2` times
    /// `(a1*d2 - a2*d1) + b1*d2*sqrt(c1) - b2*d1*sqrt(c2)`.
    #[must_use]
    pub fn cmp_sign(&self, other: &Self) -> Option<Sign> {
        let d1 = self.d.sign()?;
        let d2 = other.d.sign()?;
        if d1 == Sign::Zero || d2 == Sign::Zero {
            return None;
        }
        let p = self.a.mul(&other.d).sub(&other.a.mul(&self.d));
        let q = self.b.mul(&other.d);
        let r = other.b.mul(&self.d).neg();
        let inner = sign_two_roots(&p, &q, &self.c, &r, &other.c)?;
        Some(sign_product(sign_product(inner, d1), d2))
    }
}

/// Sign of `a + b*sqrt(c)` for `c >= 0`.
///
/// If `a` and `b*sqrt(c)` agree in sign (or one is zero), that is the
/// answer. If they disagree, the larger magnitude wins, and comparing
/// magnitudes is comparing squares: the result is `sign(a) *
/// sign(a^2 - b^2*c)`. An exact zero comes out when they cancel exactly.
///
/// Returns `None` when `T` cannot decide, or when `c` is negative (the
/// value is not real).
#[must_use]
pub fn sign_root<T: Arith>(a: &T, b: &T, c: &T) -> Option<Sign> {
    let radicand = c.sign()?;
    if radicand == Sign::Negative {
        return None;
    }
    let sa = a.sign()?;
    // b*sqrt(c) has the sign of b, unless c is zero.
    let sb = if radicand == Sign::Zero {
        Sign::Zero
    } else {
        b.sign()?
    };
    if sb == Sign::Zero {
        return Some(sa);
    }
    if sa == Sign::Zero || sa == sb {
        return Some(sb);
    }
    let dominance = a.square().sub(&b.square().mul(c)).sign()?;
    Some(sign_product(sa, dominance))
}

/// Sign of `p + q*sqrt(c) + r*sqrt(e)` for `c, e >= 0`.
///
/// Split as `u + v` with `u = p + q*sqrt(c)` and `v = r*sqrt(e)`. If they
/// agree in sign (or one is zero) that is the answer. Otherwise the result
/// is `sign(u) * sign(u^2 - v^2)`, and `u^2 - v^2 = (p^2 + q^2*c - r^2*e) +
/// 2*p*q*sqrt(c)` is again one root, decided by [`sign_root`].
#[must_use]
pub fn sign_two_roots<T: Arith>(p: &T, q: &T, c: &T, r: &T, e: &T) -> Option<Sign> {
    let su = sign_root(p, q, c)?;
    let radicand = e.sign()?;
    if radicand == Sign::Negative {
        return None;
    }
    let sv = if radicand == Sign::Zero {
        Sign::Zero
    } else {
        r.sign()?
    };
    if sv == Sign::Zero {
        return Some(su);
    }
    if su == Sign::Zero || su == sv {
        return Some(sv);
    }
    let rational = p.square().add(&q.square().mul(c)).sub(&r.square().mul(e));
    let two = T::from_f64(2.0);
    let irrational = two.mul(p).mul(q);
    let dominance = sign_root(&rational, &irrational, c)?;
    Some(sign_product(su, dominance))
}
