//! Unbounded half-space representation and explicit mesh clipping policy.

use axiolid_core::{Plane3, Scalar};

/// One side of an infinite plane.
#[derive(Debug, Clone, Copy, PartialEq)]
#[cfg_attr(
    feature = "serde",
    derive(serde::Serialize, serde::Deserialize),
    serde(deny_unknown_fields)
)]
pub struct HalfSpace {
    /// Boundary plane.
    pub boundary: Plane3,
    /// `true` selects the normal side, `false` the opposite side.
    pub agreement: bool,
}

/// Relative amount by which a subject bound is expanded before an unbounded
/// half-space participates in a finite mesh operation.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct ClipMargin(Scalar);

impl ClipMargin {
    /// Construct a positive relative margin.
    pub fn new(factor: Scalar) -> Option<Self> {
        (factor.is_finite() && factor > 0.0).then_some(Self(factor))
    }

    /// Relative factor.
    pub const fn factor(self) -> Scalar {
        self.0
    }
}
