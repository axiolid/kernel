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
use crate::interval::Interval;

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

    /// A sound `f64` interval containing the value.
    ///
    /// The top 64 mantissa bits are taken with a floor shift (the value
    /// lies within one unit of them), converted with one rounding, then
    /// scaled by an exact power of two and widened two steps each way,
    /// which covers both errors even at a binade boundary. Values whose
    /// scaled result would leave the normal `f64` range get the whole line:
    /// sound, and never decisive.
    #[must_use]
    pub fn enclosure(&self) -> Interval {
        if self.mantissa.sign() == BigSign::NoSign {
            return Interval::point(0.0);
        }
        let bits = self.mantissa.bits();
        let drop = bits.saturating_sub(64);
        let exponent = self.exponent + drop as i64;
        // top < 2^64, so the result is normal iff 2^(exponent+63) is.
        if !(-1000..=900).contains(&exponent) {
            return Interval::WHOLE;
        }
        let scale = 2f64.powi(exponent as i32);
        if drop == 0 && bits <= 53 {
            // Nothing dropped and the mantissa fits a double: the value is
            // exactly representable (the scale is a power of two within
            // the normal range), so the enclosure is a point. The filter's
            // strength depends on this: input coordinates arrive here.
            let exact = i64::try_from(&self.mantissa).expect("at most 53 bits") as f64 * scale;
            return Interval::from_bounds(exact, exact);
        }
        // BigInt >> floors (towards -infinity), so value is in
        // [top, top + 1] * 2^exponent; with nothing dropped it is exactly
        // top, and only the conversion to f64 rounds.
        let top = &self.mantissa >> drop;
        let top = i128::try_from(&top).expect("at most 65 bits");
        let low = (top as f64) * scale;
        let high = if drop == 0 {
            low
        } else {
            ((top + 1) as f64) * scale
        };
        Interval::from_bounds(low.next_down().next_down(), high.next_up().next_up())
    }

    /// `(m, e)` with the value within a relative `2^-52` of `m * 2^e`,
    /// where `m` is a double of magnitude in `[2^52, 2^53)`, or `(0, 0)`.
    ///
    /// Unlike [`Dyadic::to_f64`] this never overflows or underflows, so
    /// ratios and roots of huge or tiny values stay accurate: combine the
    /// parts first, apply the exponent last. For output only.
    #[must_use]
    pub fn approx_parts(&self) -> (f64, i64) {
        let bits = self.mantissa.bits();
        if bits == 0 {
            return (0.0, 0);
        }
        let drop = bits.saturating_sub(53);
        let top = &self.mantissa >> drop;
        let top = i64::try_from(&top).expect("at most 53 bits");
        let mut m = top as f64;
        let mut e = self.exponent + drop as i64;
        // Normalise short mantissas up to 53 bits so every caller gets the
        // same magnitude range.
        let lift = 53 - bits.min(53) as i64;
        m *= 2f64.powi(lift as i32);
        e -= lift;
        (m, e)
    }

    /// Bits in the mantissa: the cost driver of every operation.
    #[must_use]
    pub fn bits(&self) -> u64 {
        self.mantissa.bits()
    }

    /// A nearby double, for output only; never for decisions.
    ///
    /// Takes the top 64 mantissa bits, so the result is within a few ulps
    /// of the value. Overflows to an infinity and underflows to zero like
    /// any `f64` conversion.
    #[must_use]
    pub fn to_f64(&self) -> f64 {
        let bits = self.mantissa.bits();
        let drop = bits.saturating_sub(64);
        let top = &self.mantissa >> drop;
        // `top` fits in an i128 comfortably (at most 64 magnitude bits).
        let top = i128::try_from(&top).expect("at most 64 bits");
        let exponent = self.exponent + drop as i64;
        let exponent = exponent.clamp(-2000, 2000) as i32;
        (top as f64) * 2f64.powi(exponent / 2) * 2f64.powi(exponent - exponent / 2)
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

    fn from_dyadic(value: &Dyadic) -> Self {
        value.clone()
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
