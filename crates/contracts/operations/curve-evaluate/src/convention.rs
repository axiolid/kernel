//! Distance convention a provider evaluates against.

/// Which distance a [`CurveEvaluator`](crate::CurveEvaluator) measures.
///
/// Naming this is not bureaucracy. For a [`Curve3::Elevated`] the plan
/// distance and the true 3D arc length differ by `sqrt(1 + grade^2)`, so a
/// caller that assumes the wrong one misplaces a structure by 0.125 m per
/// 100 m at a 5% grade and 0.5 m per 100 m at 10%. A silent mismatch is
/// exactly the failure this contract exists to prevent, so the convention
/// is part of the contract rather than provider trivia.
///
/// [`Curve3::Elevated`]: axiolid_curve::Curve3::Elevated
#[non_exhaustive]
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub enum DistanceConvention {
    /// True arc length along the 3D curve.
    ///
    /// Distance `d` advances `d` metres along the curve as built.
    ArcLength3d,
    /// Distance along the horizontal projection.
    ///
    /// This is how an alignment is authored: both the plan and the
    /// vertical profile are functions of plan distance, and a station
    /// written on a drawing is a plan distance. It is NOT the distance a
    /// wheel travels, which is longer by the grade factor.
    PlanDistance,
    /// This provider cannot recover a distance for this curve.
    ///
    /// Reported rather than approximated: an ellipse needs elliptic
    /// integrals and a B-spline needs numeric inversion, and returning the
    /// native parameter as though it were a distance would be a lie a
    /// caller cannot detect.
    Unsupported,
}

impl DistanceConvention {
    /// Whether a distance can be evaluated at all.
    #[must_use]
    pub fn is_supported(self) -> bool {
        !matches!(self, Self::Unsupported)
    }
}
