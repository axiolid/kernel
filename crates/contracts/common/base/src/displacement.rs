//! A bound on how far a rounding step moved a boundary.

/// A finite, non-negative distance bound in model units.
///
/// Never NaN, so `Eq` is sound.
#[derive(Debug, Clone, Copy, PartialEq, PartialOrd, Default)]
pub struct DisplacementBound(f64);

impl Eq for DisplacementBound {}

impl DisplacementBound {
    /// No displacement.
    pub const ZERO: Self = Self(0.0);

    /// The bound `value`; `None` unless it is finite and at least zero.
    /// `-0.0` becomes `+0.0`.
    #[must_use]
    pub fn new(value: f64) -> Option<Self> {
        (value.is_finite() && value >= 0.0).then_some(Self(value + 0.0))
    }

    /// The bound as a number.
    #[must_use]
    pub const fn get(self) -> f64 {
        self.0
    }

    /// The larger of two bounds.
    #[must_use]
    pub fn max(self, other: Self) -> Self {
        if other.0 > self.0 {
            other
        } else {
            self
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn only_finite_non_negative_values_are_bounds() {
        for bad in [f64::NAN, f64::INFINITY, f64::NEG_INFINITY, -1.0, -1e-300] {
            assert_eq!(DisplacementBound::new(bad), None, "{bad}");
        }
        assert_eq!(
            DisplacementBound::new(0.5).map(DisplacementBound::get),
            Some(0.5)
        );
    }

    #[test]
    fn negative_zero_becomes_positive_zero() {
        let zero = DisplacementBound::new(-0.0).unwrap();
        assert!(zero.get().is_sign_positive());
        assert_eq!(zero, DisplacementBound::ZERO);
        assert_eq!(DisplacementBound::default(), DisplacementBound::ZERO);
    }

    #[test]
    fn max_takes_the_larger() {
        let a = DisplacementBound::new(1.0).unwrap();
        let b = DisplacementBound::new(2.0).unwrap();
        assert_eq!(a.max(b), b);
        assert_eq!(b.max(a), b);
        assert_eq!(DisplacementBound::ZERO.max(a), a);
    }

    #[test]
    fn bounds_are_eq_and_ordered() {
        fn is_eq<T: Eq>() {}
        is_eq::<DisplacementBound>();
        assert!(DisplacementBound::ZERO < DisplacementBound::new(1e-300).unwrap());
    }
}
