//! Integer polynomials and exact real-root isolation.
//!
//! Conic intersections lead to a quartic whose roots are nested radicals
//! at best and not radicals at all in general (casus irreducibilis). So an
//! exact root is represented the way CGAL's `Algebraic_kernel_d` does it:
//! a square-free integer polynomial plus a dyadic interval that contains
//! exactly one of its roots. Every comparison is then decided exactly:
//! by Sturm counts, by refining the intervals, or by a gcd when two roots
//! might be equal.
//!
//! Coefficients are `BigInt`. A polynomial with dyadic coefficients is
//! scaled by a power of two first ([`IntPoly::from_dyadic`]), which does
//! not move its roots.

use std::cmp::Ordering;

use num_bigint::BigInt;
use num_integer::Integer;

use axiolid_guarantees::Sign;

use crate::arith::Arith;
use crate::dyadic::Dyadic;

/// A polynomial with integer coefficients, lowest degree first, no
/// trailing zeros (the zero polynomial is empty).
#[derive(Debug, Clone, PartialEq, Eq, Hash)]
pub struct IntPoly {
    coeffs: Vec<BigInt>,
}

impl IntPoly {
    /// From coefficients, lowest degree first. Trailing zeros are dropped.
    #[must_use]
    pub fn new(mut coeffs: Vec<BigInt>) -> Self {
        while coeffs
            .last()
            .is_some_and(|c| c.sign() == num_bigint::Sign::NoSign)
        {
            coeffs.pop();
        }
        Self { coeffs }
    }

    /// From dyadic coefficients, lowest degree first, scaled by a power of
    /// two to integers. The roots are those of the dyadic polynomial.
    #[must_use]
    pub fn from_dyadic(coeffs: &[Dyadic]) -> Self {
        let lowest = coeffs
            .iter()
            .filter(|c| c.sign() != Some(Sign::Zero))
            .map(Dyadic::exponent)
            .min()
            .unwrap_or(0);
        Self::new(
            coeffs
                .iter()
                .map(|c| {
                    if c.sign() == Some(Sign::Zero) {
                        BigInt::from(0)
                    } else {
                        c.mantissa() << (c.exponent() - lowest) as u64
                    }
                })
                .collect(),
        )
    }

    /// Coefficients, lowest degree first.
    #[must_use]
    pub fn coeffs(&self) -> &[BigInt] {
        &self.coeffs
    }

    /// Degree, or `None` for the zero polynomial.
    #[must_use]
    pub fn degree(&self) -> Option<usize> {
        self.coeffs.len().checked_sub(1)
    }

    /// Whether this is the zero polynomial.
    #[must_use]
    pub fn is_zero(&self) -> bool {
        self.coeffs.is_empty()
    }

    fn lead(&self) -> &BigInt {
        self.coeffs.last().expect("non-zero polynomial")
    }

    /// Exact value at a dyadic point.
    #[must_use]
    pub fn eval(&self, x: &Dyadic) -> Dyadic {
        // Horner, exact.
        let mut acc = Dyadic::zero();
        for c in self.coeffs.iter().rev() {
            acc = acc.mul(x).add(&Dyadic::from_parts(c.clone(), 0));
        }
        acc
    }

    /// Exact sign at a dyadic point.
    #[must_use]
    pub fn sign_at(&self, x: &Dyadic) -> Sign {
        self.eval(x).sign().expect("exact")
    }

    /// The derivative.
    #[must_use]
    pub fn derivative(&self) -> Self {
        Self::new(
            self.coeffs
                .iter()
                .enumerate()
                .skip(1)
                .map(|(i, c)| c * BigInt::from(i))
                .collect(),
        )
    }

    fn content(&self) -> BigInt {
        self.coeffs.iter().fold(BigInt::from(0), |g, c| g.gcd(c))
    }

    /// Divided by the gcd of its coefficients, leading coefficient
    /// positive. Keeps coefficient growth in Sturm chains in check.
    #[must_use]
    pub fn primitive(&self) -> Self {
        if self.is_zero() {
            return self.clone();
        }
        let mut g = self.content();
        if self.lead().sign() == num_bigint::Sign::Minus {
            g = -g;
        }
        Self::new(self.coeffs.iter().map(|c| c / &g).collect())
    }

    /// Divided by the absolute gcd of its coefficients: same signs, same
    /// roots, smaller numbers. Sturm chains need this rather than
    /// [`Self::primitive`], which may flip the sign.
    fn scaled_down(&self) -> Self {
        if self.is_zero() {
            return self.clone();
        }
        let g = self.content();
        let g = if g.sign() == num_bigint::Sign::Minus {
            -g
        } else {
            g
        };
        Self::new(self.coeffs.iter().map(|c| c / &g).collect())
    }

    /// Pseudo-remainder of `self` by `divisor`: the remainder of
    /// `lead(divisor)^k * self` by `divisor`, with `k` just large enough to
    /// keep all arithmetic in the integers. Positive leading factor, so the
    /// sign structure a Sturm chain needs is preserved.
    fn pseudo_rem(&self, divisor: &Self) -> Self {
        let ld = divisor.lead().clone();
        let dd = divisor.degree().expect("non-zero divisor");
        let mut r = self.coeffs.clone();
        // lead(divisor) squared is positive, so multiply by |lead| as needed
        // while keeping the sign right: using ld^2 per step when ld < 0
        // would also work, but scaling by |ld| suffices because we only
        // eliminate the top term and the factor is positive.
        let scale = if ld.sign() == num_bigint::Sign::Minus {
            -ld.clone()
        } else {
            ld.clone()
        };
        while r.len() > dd && !r.is_empty() {
            let top = r.last().cloned().unwrap_or_default();
            if top.sign() == num_bigint::Sign::NoSign {
                r.pop();
                continue;
            }
            let shift = r.len() - 1 - dd;
            // r := scale * r - (top * sign(ld)) * x^shift * divisor
            let factor = if ld.sign() == num_bigint::Sign::Minus {
                -top
            } else {
                top
            };
            for c in r.iter_mut() {
                *c *= &scale;
            }
            for (i, dc) in divisor.coeffs.iter().enumerate() {
                r[shift + i] -= &factor * dc;
            }
            r.pop();
        }
        Self::new(r)
    }

    fn gcd_poly(&self, other: &Self) -> Self {
        let (mut a, mut b) = (self.primitive(), other.primitive());
        while !b.is_zero() {
            let r = a.pseudo_rem(&b).primitive();
            a = b;
            b = r;
        }
        a.primitive()
    }

    fn exact_div(&self, divisor: &Self) -> Self {
        // Polynomial long division known to be exact over the rationals;
        // done over the integers after scaling, then made primitive.
        let dd = divisor.degree().expect("non-zero divisor");
        let Some(n) = self.degree() else {
            return self.clone();
        };
        if n < dd {
            return Self::new(vec![]);
        }
        let ld = divisor.lead().clone();
        let mut r = self.coeffs.clone();
        let mut q = vec![BigInt::from(0); n - dd + 1];
        // Scale the dividend so each step divides exactly.
        let steps = (n - dd + 1) as u32;
        let scale = num_traits_pow(&ld, steps);
        for c in r.iter_mut() {
            *c *= &scale;
        }
        for k in (0..=n - dd).rev() {
            let top = r[k + dd].clone();
            let coef = &top / &ld;
            debug_assert_eq!(&coef * &ld, top, "exact division");
            for (i, dc) in divisor.coeffs.iter().enumerate() {
                r[k + i] -= &coef * dc;
            }
            q[k] = coef;
        }
        Self::new(q).primitive()
    }

    /// The square-free part: same real roots, each with multiplicity one.
    #[must_use]
    pub fn square_free(&self) -> Self {
        if self.degree().unwrap_or(0) < 1 {
            return self.primitive();
        }
        let g = self.gcd_poly(&self.derivative());
        if g.degree() == Some(0) {
            self.primitive()
        } else {
            self.exact_div(&g)
        }
    }

    /// Greatest common divisor, primitive with positive leading
    /// coefficient.
    #[must_use]
    pub fn gcd(&self, other: &Self) -> Self {
        if self.is_zero() {
            return other.primitive();
        }
        if other.is_zero() {
            return self.primitive();
        }
        self.gcd_poly(other)
    }

    /// The Sturm chain of a square-free polynomial.
    fn sturm(&self) -> Vec<Self> {
        let mut chain = vec![self.clone(), self.derivative()];
        loop {
            let n = chain.len();
            if chain[n - 1].is_zero() {
                chain.pop();
                break;
            }
            // -rem(p_{k-2}, p_{k-1}); a positive scaling keeps signs.
            let r = chain[n - 2].pseudo_rem(&chain[n - 1]);
            if r.is_zero() {
                break;
            }
            // Sturm needs exactly -rem up to a POSITIVE factor.
            let neg = Self::new(r.scaled_down().coeffs.iter().map(|c| -c).collect());
            chain.push(neg);
        }
        chain
    }

    /// A bound `B` (a power of two) with every real root in `(-B, B)`.
    fn root_bound(&self) -> Dyadic {
        // Cauchy: 1 + max |c_i / c_n|, rounded up to a power of two.
        let lead_bits = self.lead().bits();
        let max_bits = self.coeffs.iter().map(BigInt::bits).max().unwrap_or(0);
        let shift = max_bits.saturating_sub(lead_bits) + 2;
        Dyadic::from_parts(BigInt::from(1), shift as i64)
    }

    /// Isolate every real root of this polynomial.
    ///
    /// Returns roots in increasing order, each as a [`RealRoot`] whose
    /// interval contains exactly that root. Multiple roots are reported
    /// once. The zero polynomial has no isolated roots (it vanishes
    /// everywhere); callers must handle it before asking.
    #[must_use]
    pub fn real_roots(&self) -> Vec<RealRoot> {
        if self.degree().unwrap_or(0) < 1 {
            return Vec::new();
        }
        let sf = self.square_free();
        let chain = sf.sturm();
        let bound = sf.root_bound();
        let mut out = Vec::new();
        isolate(&sf, &chain, bound.neg(), bound, &mut out);
        out
    }
}

fn num_traits_pow(base: &BigInt, exp: u32) -> BigInt {
    let mut out = BigInt::from(1);
    for _ in 0..exp {
        out *= base;
    }
    out
}

/// Sign changes of a Sturm chain at `x`, zeros skipped.
fn variations(chain: &[IntPoly], x: &Dyadic) -> usize {
    let mut count = 0;
    let mut last = Sign::Zero;
    for p in chain {
        let s = p.sign_at(x);
        if s == Sign::Zero {
            continue;
        }
        if last != Sign::Zero && s != last {
            count += 1;
        }
        last = s;
    }
    count
}

/// Roots in `(lo, hi]` by Sturm's theorem.
fn roots_in(chain: &[IntPoly], lo: &Dyadic, hi: &Dyadic) -> usize {
    variations(chain, lo) - variations(chain, hi)
}

fn midpoint(lo: &Dyadic, hi: &Dyadic) -> Dyadic {
    lo.add(hi).mul(&Dyadic::from_parts(BigInt::from(1), -1))
}

fn isolate(poly: &IntPoly, chain: &[IntPoly], lo: Dyadic, hi: Dyadic, out: &mut Vec<RealRoot>) {
    let count = roots_in(chain, &lo, &hi);
    // An open interval must not start at a root: `lo` can be the previous
    // interval's right end, and a root there. Bisect until it does not.
    let lo_is_root = poly.sign_at(&lo) == Sign::Zero;
    match count {
        0 => {}
        1 if !lo_is_root || poly.sign_at(&hi) == Sign::Zero => {
            // Exactly one root in (lo, hi]. If it is hi itself, record it as
            // an exact point, so intervals stay open at both ends otherwise.
            if poly.sign_at(&hi) == Sign::Zero {
                out.push(RealRoot::exact(poly.clone(), hi));
            } else {
                out.push(RealRoot {
                    poly: poly.clone(),
                    lo,
                    hi,
                });
            }
        }
        _ => {
            let mid = midpoint(&lo, &hi);
            isolate(poly, chain, lo, mid.clone(), out);
            isolate(poly, chain, mid, hi, out);
        }
    }
}

/// One real root of a square-free integer polynomial.
///
/// Either an exact dyadic point (`lo == hi`), or the unique root in the
/// open interval `(lo, hi)` with the polynomial non-zero at both ends.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct RealRoot {
    poly: IntPoly,
    lo: Dyadic,
    hi: Dyadic,
}

impl RealRoot {
    fn exact(poly: IntPoly, at: Dyadic) -> Self {
        Self {
            poly,
            lo: at.clone(),
            hi: at,
        }
    }

    /// The defining square-free polynomial.
    #[must_use]
    pub fn poly(&self) -> &IntPoly {
        &self.poly
    }

    /// Current isolating interval `[lo, hi]`.
    #[must_use]
    pub fn bounds(&self) -> (&Dyadic, &Dyadic) {
        (&self.lo, &self.hi)
    }

    /// Whether the root is a known dyadic value.
    #[must_use]
    pub fn is_exact(&self) -> bool {
        self.lo == self.hi
    }

    /// Halve the isolating interval once (or pin the root exactly).
    pub fn refine(&mut self) {
        if self.is_exact() {
            return;
        }
        let mid = midpoint(&self.lo, &self.hi);
        let s_mid = self.poly.sign_at(&mid);
        if s_mid == Sign::Zero {
            self.lo = mid.clone();
            self.hi = mid;
            return;
        }
        if s_mid == self.poly.sign_at(&self.lo) {
            self.lo = mid;
        } else {
            self.hi = mid;
        }
    }

    /// Refine until the interval is no wider than `width`, or the root is
    /// exact.
    pub fn refine_to_width(&mut self, width: &Dyadic) {
        while !self.is_exact() && self.hi.sub(&self.lo).sub(width).sign() == Some(Sign::Positive) {
            self.refine();
        }
    }

    /// A double near the root (for output only; never for decisions).
    #[must_use]
    pub fn approx(&self) -> f64 {
        let mut copy = self.clone();
        // 2^-60 relative is well below f64 resolution for any root the
        // bound admits; exact roots return immediately.
        for _ in 0..200 {
            if copy.is_exact() {
                break;
            }
            let lo = copy.lo.to_f64();
            let hi = copy.hi.to_f64();
            if lo == hi || hi.next_down() <= lo {
                break;
            }
            copy.refine();
        }
        let (lo, hi) = (copy.lo.to_f64(), copy.hi.to_f64());
        lo + (hi - lo) / 2.0
    }

    /// Exact sign of `root - x` for a dyadic `x`.
    #[must_use]
    pub fn cmp_dyadic(&self, x: &Dyadic) -> Sign {
        if self.is_exact() {
            return self.lo.sub(x).sign().expect("exact");
        }
        // x outside the open interval decides immediately.
        if x.sub(&self.lo).sign() != Some(Sign::Positive) {
            return Sign::Positive;
        }
        if x.sub(&self.hi).sign() != Some(Sign::Negative) {
            return Sign::Negative;
        }
        // x strictly inside: the polynomial's sign at x relative to its
        // sign at lo says which side the root is on.
        let s = self.poly.sign_at(x);
        if s == Sign::Zero {
            return Sign::Zero;
        }
        if s == self.poly.sign_at(&self.lo) {
            // No sign change between lo and x: root is above x.
            Sign::Positive
        } else {
            Sign::Negative
        }
    }

    /// Exact sign of `self - other`.
    ///
    /// Distinct roots are separated by refinement, which terminates because
    /// they differ. Equality is decided exactly first: `g = gcd` of the two
    /// polynomials. A root of `g` in the overlap of the two isolating
    /// intervals is a root of each polynomial in each interval, and each
    /// interval holds only one, so it is both roots. No such root means
    /// they differ.
    #[must_use]
    pub fn cmp_root(&self, other: &Self) -> Sign {
        if self.is_exact() {
            return other.cmp_dyadic(&self.lo).flip();
        }
        if other.is_exact() {
            return self.cmp_dyadic(&other.lo);
        }
        if let Some(sign) = disjoint(self, other) {
            return sign;
        }
        let g = self.poly.gcd(&other.poly);
        if g.degree().unwrap_or(0) >= 1 {
            // Overlap is open at both ends, and neither end is a root of
            // g: each end is an end of one interval, where that
            // interval's polynomial -- a multiple of g -- is non-zero.
            let lo = max_dyadic(&self.lo, &other.lo);
            let hi = min_dyadic(&self.hi, &other.hi);
            if roots_in(&g.sturm(), &lo, &hi) > 0 {
                return Sign::Zero;
            }
        }
        let (mut a, mut b) = (self.clone(), other.clone());
        loop {
            a.refine();
            b.refine();
            if a.is_exact() || b.is_exact() {
                return a.cmp_root(&b);
            }
            if let Some(sign) = disjoint(&a, &b) {
                return sign;
            }
        }
    }
}

fn disjoint(a: &RealRoot, b: &RealRoot) -> Option<Sign> {
    if a.hi.sub(&b.lo).sign() != Some(Sign::Positive) {
        return Some(Sign::Negative);
    }
    if b.hi.sub(&a.lo).sign() != Some(Sign::Positive) {
        return Some(Sign::Positive);
    }
    None
}

fn max_dyadic(a: &Dyadic, b: &Dyadic) -> Dyadic {
    if a.sub(b).sign() == Some(Sign::Negative) {
        b.clone()
    } else {
        a.clone()
    }
}

fn min_dyadic(a: &Dyadic, b: &Dyadic) -> Dyadic {
    if a.sub(b).sign() == Some(Sign::Positive) {
        b.clone()
    } else {
        a.clone()
    }
}

impl PartialOrd for RealRoot {
    fn partial_cmp(&self, other: &Self) -> Option<Ordering> {
        Some(match self.cmp_root(other) {
            Sign::Negative => Ordering::Less,
            Sign::Positive => Ordering::Greater,
            _ => Ordering::Equal,
        })
    }
}
