//! What a boolean actually did, alongside the mesh it produced.
//!
//! `axiolid-overlay` returns `OverlayResult`; `axiolid-field` returns
//! `FieldEvidence`. The 3D boolean previously returned a bare `TriMesh`, which
//! made the operation most in need of diagnostics the only one without any.
//! This module closes that gap with the same mental model: *what did the kernel
//! actually do to my geometry?*

use axiolid_mesh::{AttributeFate, DropReason, TriMesh};

/// Counters describing one boolean evaluation.
///
/// Every field is a fact about the computation, never a quality verdict. A
/// caller decides whether a given count is acceptable for its domain.
#[derive(Debug, Clone, PartialEq, Default)]
#[non_exhaustive]
pub struct BooleanEvidence {
    /// Triangles in the subject operand as supplied.
    pub subject_triangles: usize,
    /// Triangles across all tool operands as supplied.
    pub tool_triangles: usize,
    /// Triangles in the result.
    pub output_triangles: usize,
    /// Connected components in the result.
    ///
    /// A difference that splits a wall into two pieces reports `2`. Callers
    /// that expect a single solid can detect the split instead of discovering
    /// it downstream in a quantity takeoff.
    pub output_components: usize,
    /// Tool operands that did not intersect the subject at all.
    ///
    /// A no-op cut is usually a placement bug upstream, but it is not an error
    /// here, so it is reported rather than rejected.
    pub disjoint_tools: usize,
    /// Sub-operations executed, for composed operands.
    ///
    /// `SymmetricDifference` composed from union, intersection, and difference
    /// reports `3`; a native implementation reports `1`. This is how a caller
    /// tells a composed path from a primitive one.
    pub sub_operations: usize,
    /// What happened to each named attribute channel on the subject.
    ///
    /// A boolean creates vertices along the cut that have no preimage in
    /// either operand, so a channel cannot always survive. Reporting the
    /// fate per channel is what turns a silent loss into an answer: the
    /// caller learns the data is gone AND why, instead of comparing the
    /// mesh before and after to find out.
    pub attribute_fates: Vec<(String, AttributeFate)>,
    /// Whether the provider detected coincident faces between operands.
    ///
    /// Coincident faces are the dominant source of cross-kernel disagreement.
    /// Reporting the encounter lets a caller treat those results with more
    /// care without Axiolid choosing a policy on its behalf.
    pub coincident_faces_encountered: bool,
    /// Whether an analytic closed-form path produced this result instead of the
    /// general boolean solver.
    ///
    /// An analytic path is exact only for the operand shapes it recognises, and
    /// it produces a different (though equally valid) triangulation from the
    /// general solver. A caller comparing results across runs, or reproducing a
    /// result elsewhere, needs to know which machinery ran -- the same reason
    /// [`Self::sub_operations`] distinguishes a composed result from a
    /// primitive one.
    pub analytic_path: bool,
    /// Relative overlap between operands, when the provider measured it.
    ///
    /// The smallest operand-overlap extent divided by the operand size, so it
    /// is scale-free: two cubes overlapping by 1mm at metre scale and by 1um
    /// at millimetre scale report the same number. `None` means the provider
    /// did not measure conditioning, never that the input was well
    /// conditioned.
    ///
    /// Construction is f64 ([ADR 0045](../adr/0045-boolean-construction-arithmetic.md)),
    /// so accuracy degrades as this approaches zero: invisible above 1e-6,
    /// smooth between 1e-6 and 1e-12, severe below 1e-12. This is a fact
    /// about the computation, not a verdict -- a caller decides what its
    /// domain tolerates.
    pub relative_overlap: Option<f64>,
}

impl BooleanEvidence {
    /// Record one completed operation.
    ///
    /// A constructor rather than struct-literal syntax because the type is
    /// `#[non_exhaustive]`: out-of-tree providers must be able to build
    /// evidence without breaking when a counter is added.
    pub fn record(
        subject_triangles: usize,
        tool_triangles: usize,
        output_triangles: usize,
        output_components: usize,
    ) -> Self {
        Self {
            subject_triangles,
            tool_triangles,
            output_triangles,
            output_components,
            disjoint_tools: 0,
            sub_operations: 1,
            coincident_faces_encountered: false,
            analytic_path: false,
            relative_overlap: None,
            attribute_fates: Vec::new(),
        }
    }

    /// Record what happened to the subject's attribute channels.
    ///
    /// A provider that carries no attributes still calls this, reporting
    /// every channel as dropped: an empty list would be indistinguishable
    /// from a subject that had no channels to begin with.
    #[must_use]
    pub fn with_attribute_fates(mut self, fates: Vec<(String, AttributeFate)>) -> Self {
        self.attribute_fates = fates;
        self
    }

    /// Record the measured relative overlap between operands.
    ///
    /// Only a provider that actually measured conditioning may call this;
    /// leaving the field `None` is the honest default for one that did not.
    pub const fn with_relative_overlap(mut self, overlap: f64) -> Self {
        self.relative_overlap = Some(overlap);
        self
    }

    /// Record that an analytic closed-form path produced this result.
    pub const fn with_analytic_path(mut self, analytic: bool) -> Self {
        self.analytic_path = analytic;
        self
    }

    /// Set the count of tools that did not meet the subject.
    pub const fn with_disjoint_tools(mut self, count: usize) -> Self {
        self.disjoint_tools = count;
        self
    }

    /// Set how many sub-operations produced this result.
    pub const fn with_sub_operations(mut self, count: usize) -> Self {
        self.sub_operations = count;
        self
    }

    /// Record that coincident faces were encountered between operands.
    pub const fn with_coincident_faces(mut self, encountered: bool) -> Self {
        self.coincident_faces_encountered = encountered;
        self
    }

    /// Merge evidence from a sub-operation into a running total.
    ///
    /// Input counts come from the outermost call, so they are kept rather than
    /// summed; output counts and flags come from the final sub-operation.
    pub fn absorb(&mut self, other: Self) {
        self.output_triangles = other.output_triangles;
        self.output_components = other.output_components;
        self.disjoint_tools += other.disjoint_tools;
        self.sub_operations += other.sub_operations;
        self.coincident_faces_encountered |= other.coincident_faces_encountered;
        // Sticky: if any sub-operation took the analytic path, the composed
        // result is not purely a general-solver product and must not claim to
        // be.
        self.analytic_path |= other.analytic_path;
        self.attribute_fates = merge_fates(&self.attribute_fates, other.attribute_fates);
        // Worst case wins: a composed result is only as well conditioned as
        // its least well conditioned sub-operation. Taking the last or the
        // best would let a clean final step mask a degenerate earlier one.
        self.relative_overlap = match (self.relative_overlap, other.relative_overlap) {
            (Some(a), Some(b)) => Some(a.min(b)),
            (Some(a), None) => Some(a),
            (None, b) => b,
        };
    }
}

/// Compose per-channel fates across two sequential steps.
///
/// Keyed by name, in `earlier`'s order. A channel `earlier` tracked that
/// `later` does not mention was not on `later`'s input, so it was lost
/// there: dropped, keeping any earlier reason. When `earlier` is empty (the
/// first step of a fold over an evidence seeded without fates) `later` is
/// taken as is.
///
/// Public so providers composing their own batch paths merge identically.
pub fn merge_fates(
    earlier: &[(String, AttributeFate)],
    later: Vec<(String, AttributeFate)>,
) -> Vec<(String, AttributeFate)> {
    if earlier.is_empty() {
        return later;
    }
    earlier
        .iter()
        .map(|(name, fate)| {
            let next = later
                .iter()
                .find(|(n, _)| n == name)
                .map(|(_, f)| f.clone())
                .unwrap_or(AttributeFate::Dropped(DropReason::ProviderLimitation));
            (name.clone(), fate.clone().then(next))
        })
        .collect()
}

/// A boolean result: the mesh plus what was done to produce it.
#[derive(Debug, Clone, PartialEq)]
#[non_exhaustive]
pub struct BooleanOutcome {
    /// Resulting solid. An empty mesh is a valid answer, not a failure:
    /// `A ∩ B` for disjoint operands is legitimately empty.
    pub mesh: TriMesh,
    /// What the provider did.
    pub evidence: BooleanEvidence,
}

impl BooleanOutcome {
    /// Pair a mesh with its evidence.
    pub const fn new(mesh: TriMesh, evidence: BooleanEvidence) -> Self {
        Self { mesh, evidence }
    }

    /// Whether the operation produced no geometry.
    pub fn is_empty(&self) -> bool {
        self.mesh.indices.is_empty()
    }
}
