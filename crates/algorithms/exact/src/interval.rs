//! Outward-rounded interval arithmetic: the fast tier.
//!
//! Every result bound is the round-to-nearest `f64` result moved one step
//! outward. Round-to-nearest errs by at most half the gap to the adjacent
//! float in the direction of the error, so one step (`next_down` /
//! `next_up`) always covers the true value. This holds in the subnormal
//! range, at binade boundaries (where the gap below is half the gap above,
//! and so is the error), and on overflow (`+inf` rounded from a finite true
//! value steps down to `f64::MAX`, which is still below it).
//!
//! NaN never escapes as a bound: `f64::min`/`max` silently discard NaN,
//! which would turn "unknown" into a confident wrong bound. Any NaN
//! collapses the interval to the whole line, which is always sound and
//! never decides a sign.

use axiolid_guarantees::Sign;

use crate::arith::Arith;

/// A closed interval `[lo, hi]` known to contain the true value.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Interval {
    lo: f64,
    hi: f64,
}

impl Interval {
    /// The whole extended real line: contains every value, decides nothing.
    pub const WHOLE: Self = Self {
        lo: f64::NEG_INFINITY,
        hi: f64::INFINITY,
    };

    /// The degenerate interval `[value, value]`, exact for finite input.
    #[must_use]
    pub fn point(value: f64) -> Self {
        if value.is_finite() {
            Self {
                lo: value,
                hi: value,
            }
        } else {
            Self::WHOLE
        }
    }

    /// Lower bound.
    #[must_use]
    pub const fn lo(self) -> f64 {
        self.lo
    }

    /// Upper bound.
    #[must_use]
    pub const fn hi(self) -> f64 {
        self.hi
    }

    /// Whether `value` lies in `[lo, hi]`.
    #[must_use]
    pub fn contains(self, value: f64) -> bool {
        self.lo <= value && value <= self.hi
    }

    /// Bounds already known to be sound, widened one step outward.
    fn outward(lo: f64, hi: f64) -> Self {
        if lo.is_nan() || hi.is_nan() {
            return Self::WHOLE;
        }
        Self {
            lo: lo.next_down(),
            hi: hi.next_up(),
        }
    }
}

impl Arith for Interval {
    fn from_f64(value: f64) -> Self {
        Self::point(value)
    }

    fn add(&self, other: &Self) -> Self {
        Self::outward(self.lo + other.lo, self.hi + other.hi)
    }

    fn sub(&self, other: &Self) -> Self {
        Self::outward(self.lo - other.hi, self.hi - other.lo)
    }

    fn mul(&self, other: &Self) -> Self {
        let products = [
            self.lo * other.lo,
            self.lo * other.hi,
            self.hi * other.lo,
            self.hi * other.hi,
        ];
        // 0 * inf is NaN; see the module note on why NaN must not reach
        // min/max.
        if products.iter().any(|p| p.is_nan()) {
            return Self::WHOLE;
        }
        let lo = products.iter().copied().fold(f64::INFINITY, f64::min);
        let hi = products.iter().copied().fold(f64::NEG_INFINITY, f64::max);
        Self::outward(lo, hi)
    }

    fn neg(&self) -> Self {
        // Negation is exact in IEEE arithmetic: no widening needed.
        Self {
            lo: -self.hi,
            hi: -self.lo,
        }
    }

    fn sign(&self) -> Option<Sign> {
        if self.lo > 0.0 {
            Some(Sign::Positive)
        } else if self.hi < 0.0 {
            Some(Sign::Negative)
        } else if self.lo == 0.0 && self.hi == 0.0 {
            // [0, 0] contains only zero. It arises from literal zero inputs
            // (and exact negation), never from a widened operation.
            Some(Sign::Zero)
        } else {
            None
        }
    }

    fn square(&self) -> Self {
        // Tighter than mul(self, self): a square is never negative, so an
        // interval straddling zero squares to [0, max].
        if self.lo >= 0.0 || self.hi <= 0.0 {
            return self.mul(self);
        }
        let top = (self.lo * self.lo).max(self.hi * self.hi);
        if top.is_nan() {
            return Self::WHOLE;
        }
        Self {
            lo: 0.0,
            hi: top.next_up(),
        }
    }
}
