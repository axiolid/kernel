//! The arithmetic both evaluation tiers share.

use axiolid_guarantees::Sign;

/// A number system an expression can be evaluated in.
///
/// Implemented by [`crate::Interval`] (fast, may be undecided) and
/// [`crate::Dyadic`] (exact, always decided). Writing an expression once,
/// generically, is what guarantees the filter and the exact fallback
/// evaluate the same polynomial.
pub trait Arith: Clone {
    /// The value of a finite `f64`.
    ///
    /// Callers must pass finite values; check with [`crate::require_finite`].
    /// Non-finite input is never silently accepted: the interval tier widens
    /// it to the whole line and the exact tier panics.
    fn from_f64(value: f64) -> Self;

    /// `self + other`.
    #[must_use]
    fn add(&self, other: &Self) -> Self;

    /// `self - other`.
    #[must_use]
    fn sub(&self, other: &Self) -> Self;

    /// `self * other`.
    #[must_use]
    fn mul(&self, other: &Self) -> Self;

    /// `-self`.
    #[must_use]
    fn neg(&self) -> Self;

    /// The sign, or `None` when this arithmetic cannot decide it.
    ///
    /// Exact arithmetics always return `Some`.
    fn sign(&self) -> Option<Sign>;

    /// An enclosure of `sqrt(self)`, for approximate arithmetics only.
    ///
    /// The interval tier uses it to evaluate nested radicals numerically,
    /// which decides far more signs than case analysis on coefficients
    /// (a coefficient that is exactly zero becomes an interval straddling
    /// zero, and case analysis stops there). Exact arithmetics return
    /// `None`: they decide by case analysis instead, never by a rounded
    /// root. Also `None` whenever `self` might be negative, so a filter
    /// never assigns a sign to a value that is not real.
    fn sqrt_enclosure(&self) -> Option<Self> {
        None
    }

    /// `self * self`.
    #[must_use]
    fn square(&self) -> Self {
        self.mul(self)
    }
}

/// Product of two signs.
#[must_use]
pub(crate) fn sign_product(left: Sign, right: Sign) -> Sign {
    match right {
        Sign::Positive => left,
        Sign::Negative => left.flip(),
        // Zero, and any future variant: a product with zero is zero, and an
        // unknown variant must not be read as decisive.
        _ => Sign::Zero,
    }
}
