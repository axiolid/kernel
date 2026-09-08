//! The `MeshBoolean` implementation.

use crate::csg::{compute_boolean, OpType};
use axiolid_contracts::{
    Backend, BackendDescriptor, BackendId, CancellationGranularity, Determinism, ExecutionOptions,
    ExecutionTarget, GeomError, GeomResult, ScratchRequirement,
};
use axiolid_core::BooleanOperator;
use axiolid_mesh::{AttributeFate, TriMesh};
use axiolid_mesh_boolean_contract::{
    symmetric_difference_via_composition, BooleanEvidence, BooleanOutcome, MeshBoolean,
};

use crate::convert::{from_manifold, six_signed_volume, to_manifold};

/// Mesh boolean backed by `boolmesh` (pure Rust, `glam`-only, MPL-2.0).
///
/// Adopted rather than written: see `docs/adr/0014`. This type owns the
/// conversion and contract enforcement; the algorithm itself is upstream's.
#[derive(Debug, Clone, Copy, Default)]
pub struct BoolmeshBoolean;

impl BoolmeshBoolean {
    /// Stable identifier used in errors and explicit provider selection.
    pub const ID: BackendId = BackendId::new("boolmesh");

    /// Construct the provider.
    pub const fn new() -> Self {
        Self
    }

    /// Subtract axis-aligned box cutters from an axis-aligned box subject using
    /// an exact closed-form construction, bypassing the general solver.
    ///
    /// Returns `Ok(None)` when the operands are outside this path's competence,
    /// so the caller falls back to [`MeshBoolean::subtract_many`] rather than
    /// receiving an approximate answer. Declining cases:
    ///
    /// * subject or any tool is not an axis-aligned box (recognised
    ///   structurally -- triangle count, corner lattice, and face planes -- not
    ///   by bounding box, since every mesh has one of those)
    /// * no tools, or no tool overlapping the subject
    /// * the induced grid would exceed `max_cells`
    ///
    /// # Why opt-in rather than automatic
    ///
    /// Dispatching on shape automatically would make output topology depend on
    /// input in a way the caller cannot predict: the same wall would produce
    /// different triangle counts depending on whether its openings happened to
    /// be axis-aligned. Callers that want the speed ask for it and handle the
    /// `None`; callers that want one predictable code path never see this.
    ///
    /// # Guarantees
    ///
    /// The result is watertight by construction (adjacent cells address shared
    /// grid vertices by integer index) and deterministic across processes
    /// (vertex identity is an ordered map, not a randomly-seeded hash map).
    ///
    /// The returned evidence has [`BooleanEvidence::analytic_path`] set, so a
    /// caller can tell this result apart from a general-solver result.
    pub fn subtract_boxes_analytic(
        &self,
        subject: &TriMesh,
        tools: &[TriMesh],
        options: &ExecutionOptions,
        max_cells: usize,
    ) -> GeomResult<Option<BooleanOutcome>> {
        options.check_cancelled()?;

        // Tolerance is scale-relative: an absolute epsilon would reject a
        // building-scale box and accept a millimetre-scale non-box.
        let d = subject.bounds().diagonal();
        let eps = d.x.max(d.y).max(d.z) * 1e-9;

        let Some(host) = crate::box_detect::recognise(subject, eps) else {
            return Ok(None);
        };
        let mut cutters = Vec::with_capacity(tools.len());
        for tool in tools {
            let Some(b) = crate::box_detect::recognise(tool, eps) else {
                return Ok(None);
            };
            cutters.push(b);
        }

        let Some(result) = crate::cellular::subtract_boxes(&host, &cutters, max_cells) else {
            return Ok(None);
        };

        // The analytic path is exact, but it is still a provider result: hold it
        // to the same orientation contract as the general path rather than
        // trusting the construction.
        check_result(&result, BooleanOperator::Difference)?;

        let tool_refs: Vec<&TriMesh> = tools.iter().collect();
        let evidence = evidence_for(subject, &tool_refs, &result, 1).with_analytic_path(true);
        Ok(Some(BooleanOutcome::new(result, evidence)))
    }
}

impl Backend for BoolmeshBoolean {
    fn descriptor(&self) -> BackendDescriptor {
        BackendDescriptor::new(Self::ID, ExecutionTarget::PortableCpu)
    }
}

/// Verify the result against the contract before handing it back.
///
/// `boolmesh` reports its own manifoldness; we additionally check orientation,
/// because a returned inside-out solid would poison every subsequent operation
/// in a subtract chain. Blamed on the backend, not the caller: the inputs were
/// already validated on the way in.
fn check_result(result: &TriMesh, operation: BooleanOperator) -> GeomResult<()> {
    // An empty result is legitimate: subtracting a tool that fully contains the
    // subject leaves nothing. Orientation is undefined for it, so stop here.
    //
    // Removing this early return is behaviour-preserving (an empty mesh sums to
    // zero, not a negative), so no test can distinguish it. It is kept because
    // it states the intent at the point of the decision.
    if result.indices.is_empty() {
        return Ok(());
    }
    let six_volume = six_signed_volume(&result.positions, &result.indices);
    if !six_volume.is_finite() {
        return Err(GeomError::BackendContractViolation {
            backend: BoolmeshBoolean::ID,
            detail: format!("{operation:?} returned non-finite signed volume"),
        });
    }
    if six_volume < 0.0 {
        return Err(GeomError::BackendContractViolation {
            backend: BoolmeshBoolean::ID,
            detail: format!(
                "{operation:?} returned an inside-out solid (signed volume {:.6})",
                six_volume / 6.0
            ),
        });
    }
    Ok(())
}

impl MeshBoolean for BoolmeshBoolean {
    /// Honest about the general path, which is not byte-reproducible.
    ///
    /// Upstream `boolmesh` dedups vertices through a randomly-seeded
    /// `HashMap`, so its output ordering varies between processes even
    /// single-threaded. That is the same root cause documented in
    /// [`Self::subtract_boxes_analytic`], which avoids it with a `BTreeMap`. The general
    /// path therefore cannot claim [`Determinism::Bitwise`] at any thread
    /// count, and enabling `parallel` does not make it weaker than it
    /// already is.
    ///
    /// [`Determinism::Topological`] is the honest ceiling: the result's
    /// connectivity is stable, its vertex ordering is not. Callers needing
    /// byte reproducibility use [`Self::subtract_boxes_analytic`],
    /// whose ordered vertex identity makes it reproducible across processes.
    fn determinism(&self) -> Determinism {
        Determinism::Topological
    }

    /// Measured, not guessed.
    ///
    /// `boolmesh` builds a Morton collider and intersection tables sized by the
    /// combined input. It exposes no bound, so this was measured directly with
    /// a counting global allocator (`src/bin/scratch_probe.rs`) across all four
    /// operations at 24 to 1,536 input triangles:
    ///
    /// ```text
    ///  triangles   peak bytes   bytes/triangle   worst operation
    ///         24        63852             2660   SymmetricDifference
    ///         96       125580             1308   SymmetricDifference
    ///        384       438844             1142   SymmetricDifference
    ///       1536      1756060             1143   SymmetricDifference
    /// ```
    ///
    /// Consumption is linear in input size, converging to roughly 1.1 KiB per
    /// triangle; the higher ratio at small inputs is fixed setup cost being
    /// divided by few triangles. `SymmetricDifference` is worst because it is
    /// composed from three passes and holds intermediates alive.
    ///
    /// The declared bound is 4 KiB per triangle: above the worst observed
    /// ratio with roughly 1.5x headroom for allocator variance and future
    /// operand shapes. A declared bound that is occasionally too low is worse
    /// than `Unbounded`, because it makes a budget look enforced when it is not.
    fn scratch_requirement(&self) -> ScratchRequirement {
        ScratchRequirement::PerElement {
            bytes_per_element: 4096,
        }
    }

    /// Group mutually disjoint cutters and remove each group with one
    /// boolean, instead of one boolean per cutter.
    ///
    /// Rests on `(S \ A) \ B == S \ (A union B)` and on a concatenation of
    /// disjoint solids being their union. Bounding-box grouping over-separates
    /// but never wrongly fuses, so the result is identical to the sequential
    /// default -- gated by volume equality in `tests/batch.rs`.
    /// `boolmesh` takes no cancellation handle, so nothing can interrupt a
    /// single boolean once it starts. Declared honestly: the batch override
    /// polls between groups, which is the only real poll point available.
    fn cancellation_granularity(&self) -> CancellationGranularity {
        CancellationGranularity::BetweenOperations
    }

    fn subtract_many(
        &self,
        subject: &TriMesh,
        tools: &[TriMesh],
        options: &ExecutionOptions,
    ) -> GeomResult<BooleanOutcome> {
        self.subtract_grouped(subject, tools, options)
    }

    fn boolean(
        &self,
        subject: &TriMesh,
        tool: &TriMesh,
        operation: BooleanOperator,
        options: &ExecutionOptions,
    ) -> GeomResult<BooleanOutcome> {
        // `SymmetricDifference` has no `boolmesh` counterpart. Compose it from
        // the three primitives rather than pretending the backend's operation
        // set is the contract's; `sub_operations` in the evidence records that
        // the caller got a composed result.
        if operation == BooleanOperator::SymmetricDifference {
            return symmetric_difference_via_composition(self, subject, tool, options);
        }

        let op = match operation {
            BooleanOperator::Union => OpType::Add,
            BooleanOperator::Intersection => OpType::Intersect,
            BooleanOperator::Difference => OpType::Subtract,
            BooleanOperator::SymmetricDifference => unreachable!("composed above"),
            // A future operand added to the contract. Refusing is honest and
            // lets the registry fall back to a provider that supports it.
            _ => {
                return Err(GeomError::Unsupported {
                    backend: BoolmeshBoolean::ID,
                    operation: axiolid_contracts::Operation::MeshBoolean,
                })
            }
        };
        options.check_cancelled()?;
        let subject_manifold = to_manifold(subject, "subject")?;
        let tool_manifold = to_manifold(tool, "tool")?;

        let output = match compute_boolean(&subject_manifold, &tool_manifold, op) {
            Ok(output) => output,
            // An empty result is a legitimate value under the contract: the
            // intersection of disjoint solids is nothing. `boolmesh` signals
            // that by failing to build a mesh from an empty vertex matrix, so
            // translate it back into the empty solid the contract specifies
            // rather than surfacing a backend quirk as a geometry error.
            Err(reason) if is_empty_result(&reason) => {
                let empty = TriMesh::new(Vec::new(), Vec::new());
                let evidence = evidence_for(subject, &[tool], &empty, 1);
                return Ok(BooleanOutcome::new(empty, evidence));
            }
            Err(reason) => {
                return Err(GeomError::Degenerate(format!(
                    "boolmesh {operation:?} failed: {reason}"
                )))
            }
        };

        let result = from_manifold(&output);
        check_result(&result, operation)?;
        let evidence = evidence_for(subject, &[tool], &result, 1);
        Ok(BooleanOutcome::new(result, evidence))
    }
}

/// Whether a `boolmesh` failure actually means "the result is empty".
///
/// Matched on the message because the backend does not expose a typed empty
/// case. Narrow on purpose: any other failure stays a real error rather than
/// being silently converted into an empty solid, which would turn a genuine
/// fault into a plausible-looking zero-volume answer.
fn is_empty_result(reason: &impl std::fmt::Display) -> bool {
    let text = reason.to_string().to_lowercase();
    text.contains("empty pos matrix") || text.contains("empty mesh")
}

/// Build evidence for one completed operation.
///
/// `output_components` counts connected components by union-find over shared
/// vertex indices. A difference that severs a wall reports two components, so
/// the caller learns about the split here rather than in a later quantity.
fn evidence_for(
    subject: &TriMesh,
    tools: &[&TriMesh],
    result: &TriMesh,
    sub_operations: usize,
) -> BooleanEvidence {
    let disjoint = tools
        .iter()
        .filter(|tool| !subject.bounds().intersects(&tool.bounds()))
        .count();
    let evidence = BooleanEvidence::record(
        subject.triangle_count(),
        tools.iter().map(|t| t.triangle_count()).sum(),
        result.triangle_count(),
        axiolid_mesh::component_count(result),
    )
    .with_disjoint_tools(disjoint)
    .with_sub_operations(sub_operations)
    // `boolmesh` does not report coincident-face encounters, so claiming
    // detection would be a lie. Left false until a provider can answer.
    .with_coincident_faces(false)
    // A cut creates vertices with no preimage in either operand, and
    // `boolmesh` returns positions and indices only. Naming each channel
    // as dropped is what makes the loss an answer rather than a silence.
    .with_attribute_fates(
        subject
            .attributes
            .iter()
            .map(|channel| {
                let reason = match channel.blend {
                    axiolid_mesh::Blend::None => axiolid_mesh::DropReason::NotBlendable,
                    _ => axiolid_mesh::DropReason::ProviderLimitation,
                };
                (channel.name.clone(), AttributeFate::Dropped(reason))
            })
            .collect(),
    );

    match relative_overlap(subject, tools) {
        Some(ratio) => evidence.with_relative_overlap(ratio),
        None => evidence,
    }
}

/// Smallest relative overlap between the subject and any tool.
///
/// Measured from operand bounds, not from the result: the result cannot show
/// how thin the sliver that produced it was. For each intersecting tool, the
/// overlap box's shortest side is divided by the operand size, giving a
/// scale-free number (ADR 0045).
///
/// Disjoint tools are skipped — they contribute no intersection to condition.
/// `None` when nothing intersects or the operands are degenerate, which is
/// honest: no intersection was constructed, so none was ill conditioned.
fn relative_overlap(subject: &TriMesh, tools: &[&TriMesh]) -> Option<f64> {
    let subject_bounds = subject.bounds();
    let mut worst: Option<f64> = None;

    for tool in tools {
        let tool_bounds = tool.bounds();
        if !subject_bounds.intersects(&tool_bounds) {
            continue;
        }

        // The overlap box: componentwise max of mins, min of maxes.
        let lo = subject_bounds.min.max(tool_bounds.min);
        let hi = subject_bounds.max.min(tool_bounds.max);
        let overlap = hi - lo;

        // Scale: the larger operand's diagonal, so the ratio is relative to
        // the model rather than to whichever operand happens to be smaller.
        let scale = subject_bounds
            .diagonal()
            .length()
            .max(tool_bounds.diagonal().length());
        if !scale.is_finite() || scale <= 0.0 {
            continue;
        }

        // The shortest overlap side governs: a wide, thin sliver is exactly
        // as ill conditioned as its thin dimension makes it.
        let thinnest = overlap.x.min(overlap.y).min(overlap.z);
        if !thinnest.is_finite() {
            continue;
        }
        let ratio = (thinnest / scale).max(0.0);
        worst = Some(worst.map_or(ratio, |w: f64| w.min(ratio)));
    }

    worst
}

/// Batch difference: fuse mutually disjoint cutters, then subtract per group.
///
/// The sequential default runs N booleans, each against a subject that has
/// already been cut N-1 times. This override runs one boolean per GROUP of
/// mutually disjoint cutters, which for the dominant layout (a wall with
/// non-overlapping openings) collapses to a single boolean.
///
/// Correctness rests on `(S \ A) \ B == S \ (A union B)`, plus the fact that a
/// concatenation of disjoint solids IS their union. Grouping uses bounding
/// boxes, which over-separate but never wrongly fuse, so the result is
/// identical to the sequential path -- asserted by the volume gates.
impl BoolmeshBoolean {
    /// Subtract every tool, grouping disjoint ones into single operations.
    pub(crate) fn subtract_grouped(
        &self,
        subject: &TriMesh,
        tools: &[TriMesh],
        options: &ExecutionOptions,
    ) -> GeomResult<BooleanOutcome> {
        let bounds: Vec<_> = tools.iter().map(TriMesh::bounds).collect();
        let groups = crate::grouping::disjoint_groups(&bounds);

        let mut current = subject.clone();
        let mut sub_operations = 0;
        for group in &groups {
            // The only real poll point: between groups. Cancelling here returns
            // no mesh at all rather than a partially cut one.
            options.check_cancelled()?;
            sub_operations += 1;
            // A single-member group gains nothing from fusing, so skip the
            // copy and subtract the tool directly.
            if let [only] = group.as_slice() {
                current = self.difference(&current, &tools[*only], options)?.mesh;
                continue;
            }
            let members: Vec<&TriMesh> = group.iter().map(|&i| &tools[i]).collect();
            let fused = crate::grouping::fuse(&members);
            current = self.difference(&current, &fused, options)?.mesh;
        }
        let borrowed: Vec<&TriMesh> = tools.iter().collect();
        let evidence = evidence_for(subject, &borrowed, &current, sub_operations);
        Ok(BooleanOutcome::new(current, evidence))
    }

    /// One difference, routed through the same validation as `boolean`.
    fn difference(
        &self,
        subject: &TriMesh,
        tool: &TriMesh,
        options: &ExecutionOptions,
    ) -> GeomResult<BooleanOutcome> {
        self.boolean(subject, tool, BooleanOperator::Difference, options)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use axiolid_core::Point3;

    /// A tetrahedron with reversed winding: the shape a faulty backend would
    /// return if it inverted its output.
    fn inside_out_tetrahedron() -> TriMesh {
        let positions = vec![
            Point3::new(0.0, 0.0, 0.0),
            Point3::new(1.0, 0.0, 0.0),
            Point3::new(0.0, 1.0, 0.0),
            Point3::new(0.0, 0.0, 1.0),
        ];
        // Verified: signed volume -1/6, i.e. inward-facing normals.
        let indices = vec![0, 1, 2, 0, 3, 1, 0, 2, 3, 1, 3, 2];
        TriMesh::new(positions, indices)
    }

    fn positive_infinite_volume() -> TriMesh {
        let large = 1.0e308;
        let small = 1.0e-308;
        TriMesh::new(
            vec![
                Point3::new(small, small, 1.0),
                Point3::new(large, 1.0, small),
                Point3::new(small, large, 1.0),
            ],
            vec![0, 1, 2],
        )
    }

    /// `check_result` guards against an upstream regression that inverts its
    /// output. With validated inputs the current `boolmesh` release never does
    /// this -- verified by instrumenting the branch across the whole suite and
    /// observing zero hits -- so the guard is exercised directly here rather
    /// than left as untested defensive code.
    #[test]
    fn an_inside_out_result_is_blamed_on_the_backend() {
        let error = check_result(&inside_out_tetrahedron(), BooleanOperator::Difference)
            .expect_err("an inside-out result must be rejected");
        match error {
            GeomError::BackendContractViolation { backend, detail } => {
                assert_eq!(backend, BoolmeshBoolean::ID);
                assert!(detail.contains("inside-out"), "{detail}");
            }
            other => panic!("must blame the backend, not the caller: {other:?}"),
        }
    }

    #[test]
    fn a_non_finite_result_is_blamed_on_the_backend() {
        let error = check_result(&positive_infinite_volume(), BooleanOperator::Union)
            .expect_err("a non-finite result volume must be rejected");
        assert!(
            matches!(error, GeomError::BackendContractViolation { backend, ref detail }
                if backend == BoolmeshBoolean::ID && detail.contains("non-finite signed volume")),
            "must blame the backend, got {error:?}"
        );
    }

    /// An empty result is legitimate (tool fully contains subject) and must not
    /// be mistaken for an orientation fault.
    #[test]
    fn an_empty_result_is_accepted() {
        assert!(check_result(&TriMesh::default(), BooleanOperator::Difference).is_ok());
    }

    /// A correctly oriented result passes.
    #[test]
    fn an_outward_result_is_accepted() {
        let mut mesh = inside_out_tetrahedron();
        for corner in mesh.indices.chunks_exact_mut(3) {
            corner.swap(1, 2);
        }
        assert!(check_result(&mesh, BooleanOperator::Difference).is_ok());
    }
}
