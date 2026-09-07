//! Validation that a timed result was actually correct.
//!
//! # Why timing alone is a lie
//!
//! A benchmark that measures only elapsed time rewards the fastest way to be
//! wrong. Three failure modes all look like wins on a stopwatch:
//!
//! - **Declining.** An implementation that returns "unsupported" for a hard
//!   input finishes instantly. It did not compute the answer faster; it did not
//!   compute it at all.
//! - **Failing silently.** A boolean that produces a zero-volume result has
//!   failed, but it failed *quickly*.
//! - **Being optimised away.** A pure computation whose result is discarded can
//!   be deleted wholesale by the optimiser, timing an empty loop.
//!
//! Every benchmark in this crate therefore states what the answer must be, in
//! terms derived from how the input was built, and fails on mismatch.

use axiolid_core::Scalar;

/// How close a measured value must be to its expected value.
///
/// Geometry accumulates floating-point error proportional to the magnitudes
/// involved, so a fixed absolute epsilon is wrong for large inputs and a fixed
/// relative one is wrong near zero. This carries both and accepts either.
#[derive(Debug, Clone, Copy)]
pub struct Tolerance {
    /// Absolute allowance, for values near zero.
    pub absolute: Scalar,
    /// Relative allowance, scaled by the larger magnitude.
    pub relative: Scalar,
}

impl Tolerance {
    /// A tolerance suitable for exact-arithmetic results, where agreement
    /// should be near machine precision.
    pub const EXACT: Self = Self {
        absolute: 1e-12,
        relative: 1e-12,
    };

    /// A tolerance suitable for accumulated floating-point geometry.
    pub const GEOMETRIC: Self = Self {
        absolute: 1e-9,
        relative: 1e-9,
    };

    /// Whether two values agree within this tolerance.
    #[must_use]
    pub fn accepts(self, measured: Scalar, expected: Scalar) -> bool {
        if !measured.is_finite() || !expected.is_finite() {
            // A non-finite result is never "close enough"; it is a failure that
            // a naive difference test would silently pass as NaN comparisons
            // are false in both directions.
            return false;
        }
        let difference = (measured - expected).abs();
        if difference <= self.absolute {
            return true;
        }
        let magnitude = measured.abs().max(expected.abs());
        difference <= self.relative * magnitude
    }
}

/// Why a benchmark's result was rejected.
#[derive(Debug, Clone, PartialEq)]
pub enum ValidationError {
    /// The implementation declined to produce a result.
    ///
    /// Kept distinct from a wrong answer: declining is a legitimate contract
    /// outcome, but it must never be reported as a fast one.
    Declined {
        /// Which measurement declined.
        what: &'static str,
    },
    /// A value disagreed with the derived ground truth.
    Mismatch {
        /// Which measurement disagreed.
        what: &'static str,
        /// What the implementation produced.
        measured: Scalar,
        /// What the construction guarantees.
        expected: Scalar,
    },
    /// An integer invariant disagreed.
    CountMismatch {
        /// Which measurement disagreed.
        what: &'static str,
        /// What the implementation produced.
        measured: i64,
        /// What the construction guarantees.
        expected: i64,
    },
}

impl core::fmt::Display for ValidationError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        match self {
            Self::Declined { what } => {
                write!(
                    f,
                    "{what}: implementation declined; a refusal is not a faster answer"
                )
            }
            Self::Mismatch {
                what,
                measured,
                expected,
            } => write!(f, "{what}: measured {measured}, expected {expected}"),
            Self::CountMismatch {
                what,
                measured,
                expected,
            } => write!(f, "{what}: measured {measured}, expected {expected}"),
        }
    }
}

impl std::error::Error for ValidationError {}

/// A checked benchmark result.
///
/// Constructing one is the only way to assert a measurement was correct, so a
/// benchmark cannot report a number without having stated what it should be.
#[derive(Debug, Clone, Copy)]
pub struct Validated;

impl Validated {
    /// Check a scalar against a derived expectation.
    ///
    /// # Errors
    ///
    /// Returns [`ValidationError::Mismatch`] when the values disagree beyond
    /// `tolerance`, including when either is non-finite.
    pub fn scalar(
        what: &'static str,
        measured: Scalar,
        expected: Scalar,
        tolerance: Tolerance,
    ) -> Result<Self, ValidationError> {
        if tolerance.accepts(measured, expected) {
            Ok(Self)
        } else {
            Err(ValidationError::Mismatch {
                what,
                measured,
                expected,
            })
        }
    }

    /// Check an integer invariant, such as an Euler characteristic.
    ///
    /// # Errors
    ///
    /// Returns [`ValidationError::CountMismatch`] when the values differ.
    pub fn count(
        what: &'static str,
        measured: i64,
        expected: i64,
    ) -> Result<Self, ValidationError> {
        if measured == expected {
            Ok(Self)
        } else {
            Err(ValidationError::CountMismatch {
                what,
                measured,
                expected,
            })
        }
    }

    /// Check that an implementation produced a result at all.
    ///
    /// # Errors
    ///
    /// Returns [`ValidationError::Declined`] when the option is empty.
    pub fn present<T>(what: &'static str, value: Option<T>) -> Result<T, ValidationError> {
        value.ok_or(ValidationError::Declined { what })
    }
}

/// Validate or panic, for use inside a benchmark's setup.
///
/// A benchmark that measures an incorrect implementation is worse than no
/// benchmark, so the failure is loud and immediate rather than a logged warning
/// that a later reader mistakes for a passing run.
///
/// # Panics
///
/// Panics with the validation error when `result` is `Err`.
pub fn expect_valid(result: Result<Validated, ValidationError>) {
    if let Err(error) = result {
        panic!("benchmark validation failed: {error}");
    }
}
