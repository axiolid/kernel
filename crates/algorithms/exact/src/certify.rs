//! The two-tier driver: interval filter first, exact fallback second.

use axiolid_guarantees::{Certified, Precision, Sign};

use crate::arith::Arith;
use crate::dyadic::Dyadic;
use crate::interval::Interval;

/// Why an exact computation was refused.
#[non_exhaustive]
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum ExactError {
    /// An input was NaN or infinite. Exact values exist only for finite
    /// numbers.
    NonFinite,
    /// A line's two defining points coincide, so it has no direction.
    DegenerateLine,
    /// A circle's radius was negative.
    NegativeRadius,
    /// Two parameters were compared along different lines, where
    /// parameters mean different things.
    DifferentLines,
    /// The expression has no real value in exact arithmetic (for example a
    /// square root of a negative number). A construction that validated its
    /// inputs never reports this.
    Undefined,
    /// More nested square roots than [`crate::tower::MAX_DEPTH`]: cost
    /// grows exponentially with depth, so the tower refuses rather than
    /// run unboundedly.
    TooDeep,
}

/// Refuse NaN and infinities before any arithmetic runs.
pub fn require_finite(values: &[f64]) -> Result<(), ExactError> {
    if values.iter().all(|value| value.is_finite()) {
        Ok(())
    } else {
        Err(ExactError::NonFinite)
    }
}

/// A sign question, written once and evaluable in any [`Arith`].
///
/// Implementors build their polynomial from `T::from_f64` of the inputs and
/// return `None` only when `T` cannot decide. The same code then runs as
/// the interval filter and as the exact fallback.
pub trait SignExpr {
    /// The sign in arithmetic `T`, or `None` when `T` cannot decide.
    fn sign_in<T: Arith>(&self) -> Option<Sign>;
}

/// The interval filter alone, so escalation can be observed and measured.
///
/// [`Certified::Uncertain`] means the exact tier is needed.
pub fn filter<E: SignExpr>(expr: &E) -> Certified {
    match expr.sign_in::<Interval>() {
        Some(sign) => Certified::Certain {
            sign,
            precision: Precision::F64,
        },
        None => Certified::Uncertain {
            attempted: Precision::F64,
        },
    }
}

/// The proven sign: filter, then exact arithmetic if the filter was
/// undecided.
///
/// Callers must have checked their inputs are finite
/// ([`require_finite`]); the exact tier has no value for NaN or infinity.
pub fn certify<E: SignExpr>(expr: &E) -> Result<Sign, ExactError> {
    if let Some(sign) = expr.sign_in::<Interval>() {
        return Ok(sign);
    }
    expr.sign_in::<Dyadic>().ok_or(ExactError::Undefined)
}
