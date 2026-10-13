//! A certified enclosure of a distance (#294).

use axiolid_core::Scalar;

/// A distance certified to lie in `[lower, upper]`: what
/// [`CurveEvaluator::path_joint_distances`](crate::CurveEvaluator::path_joint_distances)
/// returns for each joint of a curve path.
///
/// The bounds are a proof, not an estimate: the true value lies between
/// them. A distance read from stored data alone (a line's, a circle's, a
/// sum of stored lengths) is a point interval, `lower == upper`.
#[non_exhaustive]
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct DistanceInterval {
    /// The largest value the distance is certified not to be below.
    pub lower: Scalar,
    /// The smallest value the distance is certified not to exceed.
    pub upper: Scalar,
}

impl DistanceInterval {
    /// The interval `[lower, upper]`. Nothing is checked here; a provider
    /// returns `lower <= upper`, both finite.
    #[must_use]
    pub const fn new(lower: Scalar, upper: Scalar) -> Self {
        Self { lower, upper }
    }

    /// The point interval `[value, value]`: a distance known exactly.
    #[must_use]
    pub const fn point(value: Scalar) -> Self {
        Self::new(value, value)
    }

    /// `upper - lower`.
    #[must_use]
    pub fn width(&self) -> Scalar {
        self.upper - self.lower
    }

    /// The midpoint, `(lower + upper) / 2`.
    #[must_use]
    pub fn midpoint(&self) -> Scalar {
        0.5 * (self.lower + self.upper)
    }

    /// Whether `value` lies in `[lower, upper]`.
    #[must_use]
    pub fn contains(&self, value: Scalar) -> bool {
        self.lower <= value && value <= self.upper
    }

    /// Whether the two intervals share a point, their ends included:
    /// whether a distance certified to `self` can lie in `other` (#305).
    #[must_use]
    pub fn meets(&self, other: &Self) -> bool {
        self.lower <= other.upper && other.lower <= self.upper
    }

    /// Whether it is a point interval.
    #[must_use]
    pub fn is_point(&self) -> bool {
        self.lower == self.upper
    }
}
