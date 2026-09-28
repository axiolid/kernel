//! Fixed-point intervals at a chosen precision, with certified `sin` and
//! `cos`: the tier for values no finite arithmetic holds exactly.
//!
//! A [`FixedInterval`] is `[lo, hi] * 2^-bits` with big-integer bounds.
//! Every operation rounds its bounds outward, so the true value stays
//! inside. Unlike [`Dyadic`] the result is not exact, but
//! its width shrinks as `bits` grows: a sign the interval cannot decide
//! at one precision is asked again at a higher one, and any nonzero value
//! is decided at some precision. Deciding that a value is exactly zero is
//! the caller's business (by an identity, as for harmonics of dyadic
//! angles).
//!
//! `sin` and `cos` of a dyadic angle are summed from their Taylor series
//! after halving the angle below `1/16`, with the Lagrange remainder added
//! to the bounds, and doubled back up in interval arithmetic.

use num_bigint::BigInt;
use num_integer::Integer;

use axiolid_guarantees::Sign;

use crate::dyadic::Dyadic;
use crate::interval::Interval;

/// A real number between `lo * 2^-bits` and `hi * 2^-bits`.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct FixedInterval {
    lo: BigInt,
    hi: BigInt,
    bits: u32,
}

/// `x / 2^k`, rounded up.
fn ceil_shr(x: &BigInt, k: u32) -> BigInt {
    -((-x) >> k)
}

impl FixedInterval {
    /// Exactly `value` (an integer), at `bits` fractional bits.
    #[must_use]
    pub fn integer(value: i64, bits: u32) -> Self {
        let v = BigInt::from(value) << bits;
        Self {
            lo: v.clone(),
            hi: v,
            bits,
        }
    }

    /// The narrowest interval at `bits` holding `value`.
    #[must_use]
    pub fn from_dyadic(value: &Dyadic, bits: u32) -> Self {
        let shift = value.exponent() + i64::from(bits);
        let m = value.mantissa();
        if shift >= 0 {
            let v = m << (shift as u64);
            return Self {
                lo: v.clone(),
                hi: v,
                bits,
            };
        }
        let k = u32::try_from(-shift).unwrap_or(u32::MAX);
        Self {
            lo: m >> k,
            hi: ceil_shr(m, k),
            bits,
        }
    }

    /// The narrowest interval at `bits` holding a finite `value`.
    #[must_use]
    pub fn from_f64(value: f64, bits: u32) -> Option<Self> {
        Dyadic::try_from_f64(value).map(|d| Self::from_dyadic(&d, bits))
    }

    /// The fractional bits of the bounds.
    #[must_use]
    pub fn bits(&self) -> u32 {
        self.bits
    }

    fn same(&self, other: &Self) {
        debug_assert_eq!(self.bits, other.bits, "fixed intervals of one precision");
    }

    /// The sum.
    #[must_use]
    pub fn add(&self, other: &Self) -> Self {
        self.same(other);
        Self {
            lo: &self.lo + &other.lo,
            hi: &self.hi + &other.hi,
            bits: self.bits,
        }
    }

    /// The difference.
    #[must_use]
    pub fn sub(&self, other: &Self) -> Self {
        self.add(&other.neg())
    }

    /// The negation.
    #[must_use]
    pub fn neg(&self) -> Self {
        Self {
            lo: -&self.hi,
            hi: -&self.lo,
            bits: self.bits,
        }
    }

    /// The product, rounded outward.
    #[must_use]
    pub fn mul(&self, other: &Self) -> Self {
        self.same(other);
        let products = [
            &self.lo * &other.lo,
            &self.lo * &other.hi,
            &self.hi * &other.lo,
            &self.hi * &other.hi,
        ];
        let lo = products.iter().min().expect("four products");
        let hi = products.iter().max().expect("four products");
        Self {
            lo: lo >> self.bits,
            hi: ceil_shr(hi, self.bits),
            bits: self.bits,
        }
    }

    /// The product with an integer, exactly.
    #[must_use]
    pub fn mul_int(&self, factor: i128) -> Self {
        let factor = BigInt::from(factor);
        let (a, b) = (&self.lo * &factor, &self.hi * &factor);
        if a <= b {
            Self {
                lo: a,
                hi: b,
                bits: self.bits,
            }
        } else {
            Self {
                lo: b,
                hi: a,
                bits: self.bits,
            }
        }
    }

    /// The quotient by a positive integer, rounded outward.
    #[must_use]
    pub fn div_int(&self, divisor: u64) -> Self {
        let d = BigInt::from(divisor.max(1));
        Self {
            lo: self.lo.div_floor(&d),
            hi: -((-&self.hi).div_floor(&d)),
            bits: self.bits,
        }
    }

    /// The same interval at `bits` fractional bits, rounded outward.
    #[must_use]
    pub fn with_bits(&self, bits: u32) -> Self {
        if bits >= self.bits {
            let k = bits - self.bits;
            return Self {
                lo: &self.lo << k,
                hi: &self.hi << k,
                bits,
            };
        }
        let k = self.bits - bits;
        Self {
            lo: &self.lo >> k,
            hi: ceil_shr(&self.hi, k),
            bits,
        }
    }

    /// The sign every value in the interval shares: `Zero` only for the
    /// single point zero, `None` where the interval straddles zero.
    #[must_use]
    pub fn sign(&self) -> Option<Sign> {
        let zero = BigInt::from(0);
        if self.lo > zero {
            Some(Sign::Positive)
        } else if self.hi < zero {
            Some(Sign::Negative)
        } else if self.lo == zero && self.hi == zero {
            Some(Sign::Zero)
        } else {
            None
        }
    }

    /// A sound `f64` interval holding this one.
    #[must_use]
    pub fn enclosure(&self) -> Interval {
        let exponent = -i64::from(self.bits);
        let lo = Dyadic::from_parts(self.lo.clone(), exponent).enclosure();
        let hi = Dyadic::from_parts(self.hi.clone(), exponent).enclosure();
        Interval::from_bounds(lo.lo(), hi.hi())
    }

    /// Bounds on the magnitude of every value in the interval, as sound
    /// `f64`s: `(0, _)` where the interval holds zero.
    #[must_use]
    pub fn magnitude(&self) -> (f64, f64) {
        let e = self.enclosure();
        let upper = e.lo().abs().max(e.hi().abs());
        let lower = if e.lo() > 0.0 {
            e.lo()
        } else if e.hi() < 0.0 {
            -e.hi()
        } else {
            0.0
        };
        (lower, upper)
    }

    /// The interval cut to `[-1, 1]`: sound for a value known to lie
    /// there, such as a sine.
    fn clamp_unit(mut self) -> Self {
        let one = BigInt::from(1) << self.bits;
        if self.hi > one {
            self.hi = one.clone();
        }
        if self.lo < -&one {
            self.lo = -one;
        }
        self
    }

    /// Certified `sin(x)` and `cos(x)` at `bits` fractional bits; `None`
    /// for an angle beyond `2^40` in magnitude.
    #[must_use]
    pub fn sin_cos(x: &Dyadic, bits: u32) -> Option<(Self, Self)> {
        let magnitude = x.to_f64().abs();
        if magnitude > 2f64.powi(40) {
            return None;
        }
        // Halvings that bring the angle below 1/16.
        let mut halvings = 0u32;
        while magnitude * 16.0 > 2f64.powi(halvings as i32) {
            halvings += 1;
        }
        // Each doubling may double the error: guard bits for them and for
        // the series' own rounding.
        let w = bits + halvings + 32;
        let y = if x.mantissa() == &BigInt::from(0) {
            Self::integer(0, w)
        } else {
            let halved =
                Dyadic::from_parts(x.mantissa().clone(), x.exponent() - i64::from(halvings));
            Self::from_dyadic(&halved, w)
        };
        let y2 = y.mul(&y);
        let one_ulp = BigInt::from(1);
        let series = |mut term: Self, first: u64| {
            // Terms `y^(n) / n!` with alternating signs; `first` is the
            // power of the first term (0 for cos, 1 for sin).
            let mut sum = term.clone();
            let mut n = first;
            loop {
                term = term.mul(&y2).div_int((n + 1) * (n + 2)).neg();
                n += 2;
                let (_, size) = term.magnitude_ulps();
                if size <= one_ulp {
                    // The Lagrange remainder is at most this next term's
                    // magnitude: widen by it instead of adding it.
                    let pad = size + &one_ulp;
                    sum.lo -= &pad;
                    sum.hi += &pad;
                    return sum;
                }
                sum = sum.add(&term);
            }
        };
        let mut s = series(y.clone(), 1);
        let mut c = series(Self::integer(1, w), 0);
        for _ in 0..halvings {
            let s2 = s.mul(&c).mul_int(2);
            let c2 = Self::integer(1, w).sub(&s.mul(&s).mul_int(2));
            s = s2.clamp_unit();
            c = c2.clamp_unit();
        }
        Some((s.with_bits(bits), c.with_bits(bits)))
    }

    /// `(lower, upper)` bounds on the magnitude in units of `2^-bits`.
    fn magnitude_ulps(&self) -> (BigInt, BigInt) {
        let zero = BigInt::from(0);
        let (a, b) = (self.lo.clone(), self.hi.clone());
        let upper = if -&a > b { -&a } else { b.clone() };
        let lower = if a > zero {
            a
        } else if b < zero {
            -b
        } else {
            zero
        };
        (lower, upper)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn d(x: f64) -> Dyadic {
        Dyadic::try_from_f64(x).unwrap()
    }

    #[test]
    fn floors_and_ceilings_round_outward() {
        let x = FixedInterval::from_f64(-0.1, 8).unwrap();
        // -0.1 * 256 = -25.6
        assert_eq!(
            (x.lo.clone(), x.hi.clone()),
            (BigInt::from(-26), BigInt::from(-25))
        );
        let p = x.mul(&x);
        // 0.01 * 256 = 2.56; products 625..676 over 256.
        assert!(p.lo <= BigInt::from(2) && p.hi >= BigInt::from(3));
        let q = FixedInterval::integer(7, 4).div_int(3);
        // 7/3 * 16 = 37.33
        assert_eq!((q.lo, q.hi), (BigInt::from(37), BigInt::from(38)));
    }

    #[test]
    fn sine_and_cosine_of_one_hold_their_digits() {
        // floor(sin(1) * 2^200) and floor(cos(1) * 2^200), from exact
        // rational Taylor sums.
        let sin1: BigInt = "1352191738627887730370568382426388711294672994500390224088559"
            .parse()
            .unwrap();
        let cos1: BigInt = "868232330700371202471720823065555340568207867417093696894962"
            .parse()
            .unwrap();
        let (s, c) = FixedInterval::sin_cos(&d(1.0), 200).unwrap();
        assert!(s.lo <= sin1 && sin1 < s.hi, "{s:?}");
        assert!(c.lo <= cos1 && cos1 < c.hi, "{c:?}");
        assert!(&s.hi - &s.lo < BigInt::from(8), "{s:?}");
        assert!(&c.hi - &c.lo < BigInt::from(8), "{c:?}");
    }

    #[test]
    fn large_and_tiny_angles_agree_with_f64() {
        for x in [0.0, 1e-300, -3e-9, 0.75, -2.5, 7.0, 100.25, -12345.678] {
            let (s, c) = FixedInterval::sin_cos(&d(x), 120).unwrap();
            let (es, ec) = (s.enclosure(), c.enclosure());
            assert!(
                es.lo() <= x.sin() + 1e-15 && x.sin() - 1e-15 <= es.hi(),
                "{x}"
            );
            assert!(
                ec.lo() <= x.cos() + 1e-15 && x.cos() - 1e-15 <= ec.hi(),
                "{x}"
            );
            // sin^2 + cos^2 = 1, to a few units of the last place.
            let one = s.mul(&s).add(&c.mul(&c));
            let unit = BigInt::from(1) << 120u32;
            assert!(one.lo <= unit && unit <= one.hi, "{x}: {one:?}");
            assert!(&one.hi - &one.lo < BigInt::from(64), "{x}: {one:?}");
        }
        assert!(FixedInterval::sin_cos(&d(1e13), 64).is_none());
    }

    #[test]
    fn signs_are_certain_or_withheld() {
        assert_eq!(FixedInterval::integer(0, 10).sign(), Some(Sign::Zero));
        let (s, _) = FixedInterval::sin_cos(&d(1e-40), 64).unwrap();
        assert_eq!(s.sign(), None);
        let (s, _) = FixedInterval::sin_cos(&d(1e-40), 200).unwrap();
        assert_eq!(s.sign(), Some(Sign::Positive));
    }
}
