//! What a geometric audit found.

use axiolid_core::Scalar;

/// One geometric inconsistency, with enough detail to locate it.
#[derive(Debug, Clone, Copy, PartialEq)]
pub enum GeometricDefect {
    /// An edge endpoint does not lie on the edge's own 3D curve.
    VertexOffCurve {
        /// Index of the offending edge.
        edge: usize,
        /// Distance from the vertex to the curve point, in model units.
        error: Scalar,
    },
    /// A pcurve, lifted through its face surface, does not follow the edge.
    ///
    /// This is the defect a purely topological audit cannot see: the trim
    /// curve resolves and the loop closes, but the face boundary runs
    /// somewhere other than the edge it claims to trim.
    PcurveOffCurve {
        /// Index of the owning loop.
        loop_id: usize,
        /// Position of the edge use within that loop.
        use_index: usize,
        /// Worst sampled deviation, in model units.
        error: Scalar,
    },
    /// An edge's 3D curve could not be evaluated on its stated interval.
    UnevaluableCurve {
        /// Index of the offending edge.
        edge: usize,
    },
    /// A pcurve could not be evaluated on its stated interval.
    UnevaluablePcurve {
        /// Index of the owning loop.
        loop_id: usize,
        /// Position of the edge use within that loop.
        use_index: usize,
    },
    /// A face surface could not be evaluated at a pcurve sample.
    UnevaluableSurface {
        /// Index of the owning loop.
        loop_id: usize,
        /// Position of the edge use within that loop.
        use_index: usize,
    },
}

/// Result of a geometric audit.
///
/// Unlike the topological audit this is tolerance-dependent, so a clean
/// result is a statement about agreement AT A TOLERANCE, not an absolute one.
#[derive(Debug, Clone, Default, PartialEq)]
pub struct GeometricHealth {
    defects: Vec<GeometricDefect>,
}

impl GeometricHealth {
    /// Record one defect.
    pub(crate) fn push(&mut self, defect: GeometricDefect) {
        self.defects.push(defect);
    }

    /// Every defect found, in discovery order.
    #[must_use]
    pub fn defects(&self) -> &[GeometricDefect] {
        &self.defects
    }

    /// Whether the geometry is consistent at the audited tolerance.
    #[must_use]
    pub fn is_consistent(&self) -> bool {
        self.defects.is_empty()
    }

    /// The largest positional disagreement found, if any was measured.
    ///
    /// Useful for reporting how far off a rejected solid was, rather than
    /// only that it was rejected.
    #[must_use]
    pub fn worst_error(&self) -> Option<Scalar> {
        self.defects
            .iter()
            .filter_map(|defect| match defect {
                GeometricDefect::VertexOffCurve { error, .. }
                | GeometricDefect::PcurveOffCurve { error, .. } => Some(*error),
                _ => None,
            })
            .max_by(Scalar::total_cmp)
    }
}
