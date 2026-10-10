//! Small operation values shared by representation and algorithm contracts.

/// Regularized boolean set operation on solids.
///
/// The operand set matches `axiolid-overlay`'s planar contract so 2D and 3D
/// booleans describe the same algebra. It is deliberately *not* a mirror of any
/// one backend's operation enum: `SymmetricDifference` exists here because the
/// set algebra has it, not because a provider offered it.
///
/// Marked `#[non_exhaustive]` so a future operand cannot break downstream
/// `match` arms.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
#[non_exhaustive]
#[cfg_attr(
    feature = "serde",
    derive(serde::Serialize, serde::Deserialize),
    serde(deny_unknown_fields)
)]
pub enum BooleanOperator {
    /// Points in either operand.
    Union,
    /// Points in both operands.
    Intersection,
    /// Ordered: points in the subject and not in the tool.
    Difference,
    /// Points in exactly one operand, equal to `(A ∪ B) \ (A ∩ B)`.
    SymmetricDifference,
}

impl BooleanOperator {
    /// Every operand in a stable, declared order.
    ///
    /// Conformance suites iterate this so a new operand is automatically
    /// covered rather than silently untested.
    pub const ALL: [Self; 4] = [
        Self::Union,
        Self::Intersection,
        Self::Difference,
        Self::SymmetricDifference,
    ];

    /// Whether swapping the operands leaves the result unchanged.
    ///
    /// `Difference` is the only ordered operand; the identity is part of the
    /// public contract because callers rely on it to reorder work.
    pub const fn is_commutative(self) -> bool {
        !matches!(self, Self::Difference)
    }
}
