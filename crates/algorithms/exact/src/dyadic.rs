//! Dyadic rationals: the exact tier.
//!
//! A dyadic number is `mantissa * 2^exponent` with a big-integer mantissa.
//! Every finite `f64` is one exactly, and the set is closed under `+`, `-`
//! and `*`, so any polynomial in `f64` inputs has an exact dyadic value and
//! an exact sign. Division is deliberately absent (ADR 0068): callers clear
//! denominators instead.
//!
//! Values are kept normalised -- the mantissa odd, or zero with exponent 0 --
//! so equal values have equal representations and mantissas stay as short
//! as the value allows.

use num_bigint::{BigInt, Sign as BigSign};

use axiolid_guarantees::Sign;

use crate::arith::Arith;

/// An exact value `mantissa * 2^exponent`.
#[derive(Debug, Clone, PartialEq, Eq, Hash)]
pub struct Dyadic {
    mantissa: BigInt,
    exponent: i64,
}

impl Dyadic {
    /// Zero.
    #[must_use]
    pub fn zero() -> Self {
        Self {
            mantissa: BigInt::from(0),
            exponent: 0,
        }
    }

    /// `mantissa * 2^exponent`, normalised.
    #[must_use]
    pub fn from_parts(mantissa: BigInt, exponent: i64) -> Self {
        Self { mantissa, exponent }.normalised()
    }

    /// The exact value of a finite `f64`, or `None` for NaN and infinities.
    #[must_use]
    pub fn try_from_f64(value: f64) -> Option<Self> {
        if !value.is_finite() {
            return None;
        }
        let bits = value.to_bits();
        let negative = bits >> 63 == 1;
        let biased = ((bits >> 52) & 0x7ff) as i64;
        let fraction = bits & ((1u64 << 52) - 1);
        // Subnormals have no implicit leading bit and a fixed exponent.
        let (magnitude, exponent) = if biased == 0 {
            (fraction, -1074)
        } else {
            (fraction | (1u64 << 52), biased - 1075)
        };
        let mut mantissa = BigInt::from(magnitude);
        if negative {
            mantissa = -mantissa;
        }
        Some(Self::from_parts(mantissa, exponent))
    }

    /// The mantissa of the normalised form.
    #[must_use]
    pub fn mantissa(&self) -> &BigInt {
        &self.mantissa
    }

    /// The power-of-two exponent of the normalised form.
    #[must_use]
    pub fn exponent(&self) -> i64 {
        self.exponent
    }

    /// Bits in the mantissa: the cost driver of every operation.
    #[must_use]
    pub fn bits(&self) -> u64 {
        self.mantissa.bits()
    }

    fn normalised(mut self) -> Self {
        match self.mantissa.trailing_zeros() {
            // Only zero has no trailing-zero count.
            None => Self::zero(),
            Some(0) => self,
            Some(shift) => {
                self.mantissa >>= shift;
                self.exponent += shift as i64;
                self
            }
        }
    }

    /// Both mantissas on the smaller exponent.
    fn aligned(&self, other: &Self) -> (BigInt, BigInt, i64) {
        if self.exponent <= other.exponent {
            let shift = (other.exponent - self.exponent) as u64;
            (
                self.mantissa.clone(),
                &other.mantissa << shift,
                self.exponent,
            )
        } else {
            let shift = (self.exponent - other.exponent) as u64;
            (
                &self.mantissa << shift,
                other.mantissa.clone(),
                other.exponent,
            )
        }
    }
}

impl Arith for Dyadic {
    /// # Panics
    ///
    /// On NaN or an infinity, which have no exact value. The public entry
    /// points reject them first ([`crate::require_finite`]).
    fn from_f64(value: f64) -> Self {
        Self::try_from_f64(value).expect("exact arithmetic needs finite input")
    }

    fn add(&self, other: &Self) -> Self {
        let (left, right, exponent) = self.aligned(other);
        Self::from_parts(left + right, exponent)
    }

    fn sub(&self, other: &Self) -> Self {
        let (left, right, exponent) = self.aligned(other);
        Self::from_parts(left - right, exponent)
    }

    fn mul(&self, other: &Self) -> Self {
        // Odd times odd is odd, so a product of normalised values is already
        // normalised unless it is zero. `normalised` finds the lowest set
        // bit in the first limb for an odd mantissa, so this costs nothing
        // measurable and keeps zero canonical.
        Self {
            mantissa: &self.mantissa * &other.mantissa,
            exponent: self.exponent + other.exponent,
        }
        .normalised()
    }

    fn neg(&self) -> Self {
        Self {
            mantissa: -&self.mantissa,
            exponent: self.exponent,
        }
    }

    fn sign(&self) -> Option<Sign> {
        Some(match self.mantissa.sign() {
            BigSign::Plus => Sign::Positive,
            BigSign::Minus => Sign::Negative,
            BigSign::NoSign => Sign::Zero,
        })
    }
}
