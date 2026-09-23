//! Explicit repair plans and reports.

/// One opt-in repair. There is deliberately no `All` variant.
#[non_exhaustive]
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub enum RepairAction {
    /// Merge vertices within tolerance.
    WeldVertices,
    /// Make connected-face orientation consistent.
    UnifyOrientation,
    /// Remove degenerate elements.
    DropDegenerateElements,
    /// Flip a closed shell that encloses negative volume.
    ///
    /// Separate from [`Self::UnifyOrientation`]: that one makes neighbours
    /// agree, this one decides which way round the agreed shell faces.
    OrientOutward,
}

/// Ordered caller-approved repair plan.
#[derive(Debug, Clone, Default, PartialEq, Eq)]
pub struct RepairPlan {
    /// Repairs to attempt in order.
    pub actions: Vec<RepairAction>,
}

use axiolid_mesh::AttributeFate;

/// Audit report returned with repaired geometry.
#[derive(Debug, Clone, Default, PartialEq, Eq)]
pub struct RepairReport {
    /// Repairs actually applied.
    pub applied: Vec<RepairAction>,
    /// Repairs requested but not applicable.
    pub skipped: Vec<RepairAction>,
    /// What happened to each named attribute channel, in input order.
    ///
    /// Every input channel appears exactly once, so a caller can tell a
    /// channel that survived from one that never existed (#114).
    pub attribute_fates: Vec<(String, AttributeFate)>,
}
