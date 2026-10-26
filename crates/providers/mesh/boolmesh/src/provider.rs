//! The `MeshBoolean` implementation.
//!
//! # Whose fault an error is
//!
//! Input faults are the caller's: an operand that fails the orientation
//! gate (`InvalidInput`, `Degenerate`, naming which argument) or that the
//! kernel cannot build as a manifold (`NotManifold`). A result that fails
//! [`check_result`] is `BackendContractViolation`: the inputs were already
//! validated, so a bad result is this provider's defect, never blamed on
//! the caller. A failure inside the kernel's own solve (for example an odd
//! edge-point count in `pair_up`, #101) is `BackendContractViolation` for
//! the same reason: the operands passed every input gate, so the solve
//! failing on them is this provider's defect.

use std::borrow::Cow;

use crate::csg::{compute_boolean, OpType, SolveError};
use axiolid_contracts::{
    Backend, BackendDescriptor, BackendId, CancellationGranularity, Determinism, ExecutionOptions,
    ExecutionTarget, GeomError, GeomResult, ScratchRequirement,
};
use axiolid_core::BooleanOperator;
use axiolid_mesh::{AttributeChannel, AttributeFate, Blend, DropReason, TriMesh};
use axiolid_mesh_boolean_contract::{merge_fates, BooleanEvidence, BooleanOutcome, MeshBoolean};

use crate::attributes::{carry, sources_by_plane, FaceSource};
use crate::convert::{from_boolean_mesh, six_signed_volume, to_manifold};

/// Mesh boolean built on the algorithm absorbed from `boolmesh` (pure Rust,
/// `glam`-only, MPL-2.0).
///
/// Adopted in `docs/adr/0014` and absorbed into this crate's private `csg`
/// module by `docs/adr/0047`. This type owns conversion and contract
/// enforcement; no `csg` type is part of the public API, so the kernel can
/// be replaced without a consumer noticing.
#[derive(Debug, Clone, Copy, Default)]
pub struct BoolmeshBoolean;

/// Fixed scratch of one boolean; see [`MeshBoolean::scratch_requirement`].
pub(crate) const SCRATCH_BASE_BYTES: usize = 64 * 1024;
/// Scratch per input triangle, all operands together.
pub(crate) const SCRATCH_BYTES_PER_TRIANGLE: usize = 1536;
/// A rayon worker's bookkeeping and the work stolen onto it.
#[cfg(any(feature = "parallel", feature = "parallel-batch"))]
const WORKER_BYTES: usize = 16 * 1024;
#[cfg(not(any(feature = "parallel", feature = "parallel-batch")))]
const WORKER_BYTES: usize = 0;
/// Per worker thread. The booleans a batch runs concurrently are charged by
/// the batch, per boolean (`budget.rs`), not here: stolen work can put more
/// of them in flight than there are workers (#259).
pub(crate) const SCRATCH_BYTES_PER_WORKER: usize = WORKER_BYTES;

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

        let Some(mut result) = crate::cellular::subtract_boxes(&host, &cutters, max_cells) else {
            return Ok(None);
        };

        // The analytic path is exact, but it is still a provider result: hold it
        // to the same orientation contract as the general path rather than
        // trusting the construction.
        check_result(&result, BooleanOperator::Difference)?;

        let tool_refs: Vec<&TriMesh> = tools.iter().collect();
        let fates = carry_analytic(subject, tools, &mut result);
        let evidence = evidence_for(subject, &tool_refs, &result, 1)
            .with_analytic_path(true)
            .with_attribute_fates(fates);
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

    /// Measured, not guessed: [`ScratchRequirement::Affine`] over the input
    /// triangles of all operands and the worker threads the call runs on.
    ///
    /// `boolmesh` builds a Morton collider and intersection tables sized by the
    /// combined input. It exposes no bound, so this was measured directly with
    /// a counting global allocator (`src/bin/scratch_probe/`) across all four
    /// operations at 24 to 1,536 input triangles, after one discarded warmup
    /// call so the first operation measured is not charged for process
    /// startup (#110). With rayon each boolean runs in a pool started inside
    /// the measured window, with the join halves forced onto other workers
    /// (#226; the module docs of `scratch.rs` say how). Worst peak in bytes
    /// over the four operations and ten runs with `parallel`, six on an idle
    /// 20-core machine and four at load average 20 to 43:
    ///
    /// ```text
    ///  triangles   1 worker   4 workers   16 workers   64 workers
    ///         24      59752       65656       160136       471736
    ///         96     108232      112152       193504       460672
    ///        384     348792      357000       422472       680384
    ///       1536    1354692     1362900      1411548      1702180
    /// ```
    ///
    /// Three terms, each with its own evidence:
    ///
    /// - Per triangle, 1,536 bytes. The slope is linear and flat at about
    ///   880 bytes per triangle from 384 up to 24,576 triangles (one worker,
    ///   checked once at the two larger sizes), so this is 1.7x it.
    ///   `SymmetricDifference` is the worst operation because it is composed
    ///   from three passes and holds intermediates alive.
    /// - Base, 64 KiB: fixed setup that dominates small inputs. 24 triangles
    ///   peak at 59,752 bytes on one worker, 22 KiB above the per-triangle
    ///   term, so the base is 2.9x that excess.
    /// - Per worker, 16 KiB with `parallel`. A rayon worker's own
    ///   bookkeeping (deque, registry slot, thread-locals) is 4.8 KB at most
    ///   from 1 to 512 threads, and it is charged to whichever boolean runs
    ///   while the worker starts, which on a loaded machine can be any of
    ///   them. Work stolen onto it adds its scratch alongside the rest. The
    ///   measured excess over one worker is at most 6.7 KB per extra
    ///   worker, idle or loaded, so this is 2.4x it. With `parallel-batch`
    ///   alone the batch's pool workers pay the same bookkeeping. Without
    ///   either feature nothing runs off the calling thread and the term is
    ///   zero.
    ///
    /// # The result is charged in the solve (#308)
    ///
    /// This counts input triangles only, so it is the admission bound on
    /// the operands, not a bound on the call. Operands that cross many
    /// times have a far larger result: two perpendicular combs of 32 teeth,
    /// 1,048 triangles, cross 17,716 times, and their union peaked at 14.7
    /// MB, 8.8x this declaration (52.9 MB, 31x, for the symmetric
    /// difference). With a memory budget the solve therefore charges the
    /// result once its size is known: after the intersection pass,
    /// `boolean45` counts the result's half-edges exactly and refuses with
    /// `BudgetExceeded` before allocating the result when this declaration
    /// plus 448 bytes per half-edge beyond three per input triangle does not
    /// fit. Ordinary results stay within the three and are charged nothing
    /// more; the intersection pass stops as soon as its crossings alone
    /// could not fit. A budget that admits this declaration can still be
    /// refused, so, as for batches, the check before dispatch is a
    /// necessary condition, not a sufficient one. A composed
    /// `SymmetricDifference` charges each step the union and intersection
    /// it holds. `src/budget.rs` gives the measurements.
    ///
    /// # Batches
    ///
    /// The bound covers each boolean of [`MeshBoolean::subtract_many`] and
    /// [`MeshBoolean::union_many`] on that step's own operands, not the batch
    /// in its input triangles: intermediates can outgrow the inputs, and
    /// with `parallel-batch` stolen work can put every pair of a level in
    /// flight at once (#259). With a memory budget the batch charges each
    /// step when its operands exist (its own copies and intermediates, plus
    /// this bound for every boolean that can be live) and refuses with
    /// `BudgetExceeded` before the first step that does not fit, so the
    /// check before dispatch is a necessary condition for a batch, not a
    /// sufficient one. Each boolean then gets the budget the batch does not
    /// hold (split evenly between the pairs of a `parallel-batch` level)
    /// and charges its result against it (#308). `tests/scratch_bound.rs`
    /// holds the peak of both batches under forced schedules to the
    /// smallest budget they admit.
    ///
    /// The old `PerElement { bytes_per_element: 4096 }` was 2x the large-input
    /// slope but no bound at all for small inputs on a pool: 24 triangles
    /// on 16 workers peak at 125,688 to 160,136 bytes against its 98,304
    /// (#226).
    /// `tests/scratch_bound.rs` measures the same workloads on 1 to 64
    /// workers and fails if any peak exceeds this. A declared bound that is
    /// occasionally too low is worse than `Unbounded`, because it makes a
    /// budget look enforced when it is not.
    ///
    /// The worker count is only known where the call runs, so with a memory
    /// budget every boolean re-checks this bound against the pool it is in
    /// and refuses with `BudgetExceeded` when that pool is wider than the
    /// budget check before dispatch assumed. The re-check runs in every
    /// build, so a budget given to the provider directly, without dispatch,
    /// is enforced too.
    fn scratch_requirement(&self) -> ScratchRequirement {
        ScratchRequirement::Affine {
            base_bytes: SCRATCH_BASE_BYTES,
            bytes_per_element: SCRATCH_BYTES_PER_TRIANGLE,
            bytes_per_worker: SCRATCH_BYTES_PER_WORKER,
        }
    }

    /// `boolmesh` takes no cancellation handle, so nothing can interrupt a
    /// single boolean once it starts. Declared honestly: the batch override
    /// polls between groups, which is the only real poll point available.
    fn cancellation_granularity(&self) -> CancellationGranularity {
        CancellationGranularity::BetweenOperations
    }

    /// Group mutually disjoint cutters and remove each group with one
    /// boolean, instead of one boolean per cutter.
    ///
    /// Rests on `(S \ A) \ B == S \ (A union B)` and on a concatenation of
    /// disjoint solids being their union. Bounding-box grouping over-separates
    /// but never wrongly fuses, so the result is identical to the sequential
    /// default -- gated by volume equality in `tests/batch.rs`.
    fn subtract_many(
        &self,
        subject: &TriMesh,
        tools: &[TriMesh],
        options: &ExecutionOptions,
    ) -> GeomResult<BooleanOutcome> {
        self.subtract_grouped(subject, tools, options)
    }

    fn union_many(
        &self,
        solids: &[TriMesh],
        options: &ExecutionOptions,
    ) -> GeomResult<BooleanOutcome> {
        self.union_tree(solids, options)
    }

    fn boolean(
        &self,
        subject: &TriMesh,
        tool: &TriMesh,
        operation: BooleanOperator,
        options: &ExecutionOptions,
    ) -> GeomResult<BooleanOutcome> {
        self.boolean_within(
            subject,
            tool,
            operation,
            options,
            false,
            options.memory_budget_bytes(),
        )
    }
}

impl BoolmeshBoolean {
    /// Alternative to [`MeshBoolean::boolean`] using a cheaper winding-number
    /// classification for the same result.
    ///
    /// See `csg::boolean03::kernel03::winding03_fast`'s doc comment for the
    /// guarantee: not an approximation, but a bug in edge-break detection
    /// would mislabel a whole connected component instead of one vertex.
    /// Opt-in rather than the default until this has run against a wider
    /// correctness corpus than the differential tests already covering it.
    ///
    /// `SymmetricDifference` is refused rather than silently composed from
    /// three slow-path calls: a caller asking for the fast path should get
    /// it or a clear refusal, not a surprise fallback.
    pub fn boolean_fast(
        &self,
        subject: &TriMesh,
        tool: &TriMesh,
        operation: BooleanOperator,
        options: &ExecutionOptions,
    ) -> GeomResult<BooleanOutcome> {
        if operation == BooleanOperator::SymmetricDifference {
            return Err(GeomError::Unsupported {
                backend: BoolmeshBoolean::ID,
                operation: axiolid_contracts::Operation::MeshBoolean,
            });
        }
        self.boolean_within(
            subject,
            tool,
            operation,
            options,
            true,
            options.memory_budget_bytes(),
        )
    }

    /// One boolean with `budget` bytes for it alone (the caller's budget
    /// less whatever the caller holds), or none.
    ///
    /// Admitted on the declaration for the pool this call runs on (#226),
    /// then charged for its result inside the solve (#308, `budget.rs`).
    fn boolean_within(
        &self,
        subject: &TriMesh,
        tool: &TriMesh,
        operation: BooleanOperator,
        options: &ExecutionOptions,
        fast_winding: bool,
        budget: Option<usize>,
    ) -> GeomResult<BooleanOutcome> {
        let elements = subject.triangle_count() + tool.triangle_count();
        let result_budget = crate::budget::ResultBudget::admit(budget, elements)?;
        // `SymmetricDifference` has no `boolmesh` counterpart. Compose it from
        // the three primitives rather than pretending the backend's operation
        // set is the contract's; `sub_operations` in the evidence records that
        // the caller got a composed result.
        if operation == BooleanOperator::SymmetricDifference {
            return self.symmetric_difference(subject, tool, options, budget);
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

        let output = match compute_boolean(
            &subject_manifold,
            &tool_manifold,
            op,
            fast_winding,
            result_budget.as_ref(),
        ) {
            Ok(output) => output,
            Err(SolveError::OverBudget) => {
                return Err(GeomError::BudgetExceeded { resource: "memory" })
            }
            // An empty result is a legitimate value under the contract: the
            // intersection of disjoint solids is nothing. `boolmesh` signals
            // that by failing to build a mesh from an empty vertex matrix, so
            // translate it back into the empty solid the contract specifies
            // rather than surfacing a backend quirk as a geometry error.
            Err(SolveError::Failed(reason)) if is_empty_result(&reason) => {
                let mut empty = TriMesh::new(Vec::new(), Vec::new());
                // Zero faces carry every channel vacuously: nothing lost,
                // nothing derived. Keeps a union/fold over it consistent.
                let fates = carry(subject, tool, &mut empty, &[]);
                let evidence =
                    evidence_for(subject, &[tool], &empty, 1).with_attribute_fates(fates);
                return Ok(BooleanOutcome::new(empty, evidence));
            }
            Err(SolveError::Failed(reason)) => return Err(solve_refusal(operation, &reason)),
        };

        let mut result = from_boolean_mesh(&output);
        check_result(&result, operation)?;
        let sources: Vec<FaceSource> = output
            .src
            .iter()
            .map(|&(operand, triangle)| FaceSource { operand, triangle })
            .collect();
        let fates = carry(subject, tool, &mut result, &sources);
        let evidence = evidence_for(subject, &[tool], &result, 1).with_attribute_fates(fates);
        Ok(BooleanOutcome::new(result, evidence))
    }

    /// `A △ B = (A ∪ B) \ (A ∩ B)`, charging each step what the earlier
    /// steps leave alive (#308).
    ///
    /// The same composition as the contract's
    /// `symmetric_difference_via_composition`, with the same evidence, but
    /// that helper passes every step the whole budget while the union is
    /// held through the intersection and both through the difference. A
    /// union of crossing operands can be larger than either, so here each
    /// step gets the budget less the results it holds.
    fn symmetric_difference(
        &self,
        subject: &TriMesh,
        tool: &TriMesh,
        options: &ExecutionOptions,
        budget: Option<usize>,
    ) -> GeomResult<BooleanOutcome> {
        let held = |outcome: &BooleanOutcome| {
            crate::budget::heap_bytes(&outcome.mesh)
                + crate::budget::fates_bytes(&outcome.evidence.attribute_fates)
        };
        let less = |budget: Option<usize>, bytes: usize| budget.map(|b| b.saturating_sub(bytes));
        options.check_cancelled()?;
        let union = self.boolean_within(
            subject,
            tool,
            BooleanOperator::Union,
            options,
            false,
            budget,
        )?;
        options.check_cancelled()?;
        let budget = less(budget, held(&union));
        let intersection = self.boolean_within(
            subject,
            tool,
            BooleanOperator::Intersection,
            options,
            false,
            budget,
        )?;

        // A ∩ B empty means the operands are disjoint, so A △ B == A ∪ B,
        // and subtracting an empty solid is a degenerate operand.
        if intersection.mesh.indices.is_empty() {
            let mut evidence = union.evidence;
            evidence.sub_operations = 2;
            evidence.coincident_faces_encountered |=
                intersection.evidence.coincident_faces_encountered;
            return Ok(BooleanOutcome::new(union.mesh, evidence));
        }

        options.check_cancelled()?;
        let budget = less(budget, held(&intersection));
        let difference = self.boolean_within(
            &union.mesh,
            &intersection.mesh,
            BooleanOperator::Difference,
            options,
            false,
            budget,
        )?;

        let mut evidence = difference.evidence;
        // The result's channels came through union then difference: compose
        // those fates.
        evidence.attribute_fates =
            merge_fates(&union.evidence.attribute_fates, evidence.attribute_fates);
        evidence.subject_triangles = subject.triangle_count();
        evidence.tool_triangles = tool.triangle_count();
        evidence.sub_operations = 3;
        evidence.coincident_faces_encountered |= union.evidence.coincident_faces_encountered
            || intersection.evidence.coincident_faces_encountered;
        Ok(BooleanOutcome::new(difference.mesh, evidence))
    }
}

/// The error for a refusal inside the solve.
///
/// Every operand reaching the solve has passed the input gates, so a
/// refusal there (an unpairable edge-point count in `pair_up`, #101) is this
/// provider's defect: `BackendContractViolation`, never `Degenerate`, which
/// would blame operands that were admissible.
fn solve_refusal(operation: BooleanOperator, reason: &str) -> GeomError {
    GeomError::BackendContractViolation {
        backend: BoolmeshBoolean::ID,
        detail: format!("{operation:?} failed inside the solve: {reason}"),
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
    // Fallback fate, overwritten by every path that carries channels (all
    // of them since #116). What remains here is the answer when no
    // provenance exists: the `carry` refusal path, and nothing else.
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

/// Carry channels through the analytic box path.
///
/// The cellular result is rebuilt, not cut, so it has no `Tref`s. Every face
/// lies in a host face (as is) or a cutter face (reversed); plane lookup
/// recovers the source, and the pairwise sampler does the rest. The cutters
/// are fused so the sampler sees one tool, as in the general path.
fn carry_analytic(
    subject: &TriMesh,
    tools: &[TriMesh],
    result: &mut TriMesh,
) -> Vec<(String, AttributeFate)> {
    let members: Vec<&TriMesh> = tools.iter().collect();
    let tool = crate::grouping::fuse_with_channels(&members, &subject.attributes);
    match sources_by_plane(result, &[(subject, 1.0), (&tool, -1.0)]) {
        Some(sources) => carry(subject, &tool, result, &sources),
        // A face with no source means the construction is not what the
        // lookup assumes: attach nothing rather than a wrong value.
        None => dropped_fates(subject, DropReason::ProviderLimitation),
    }
}

/// Every subject channel as dropped, with `reason` (`NotBlendable` wins for
/// a `Blend::None` channel regardless: that is the data's own answer).
fn dropped_fates(subject: &TriMesh, reason: DropReason) -> Vec<(String, AttributeFate)> {
    subject
        .attributes
        .iter()
        .map(|c| {
            let r = match c.blend {
                Blend::None => DropReason::NotBlendable,
                _ => reason,
            };
            (c.name.clone(), AttributeFate::Dropped(r))
        })
        .collect()
}

/// Per-channel fates, in the subject's channel order.
type Fates = Vec<(String, AttributeFate)>;

/// Every subject channel `Preserved`: the starting point a composed path
/// merges each step's fates onto.
fn seed_fates(subject: &TriMesh) -> Fates {
    subject
        .attributes
        .iter()
        .map(|c| (c.name.clone(), AttributeFate::Preserved))
        .collect()
}

/// `mesh` carrying exactly `template`'s channels: its own where one matches
/// (name, width, blend), an all-`UNMAPPED` channel where none does.
///
/// Union is symmetric but channel reporting is keyed on the first operand,
/// so a tree node's left child must carry every channel the batch reports
/// on -- whichever solid it came from. Unmapped keeps the absence honest.
///
/// Built by fusing the one mesh, which copies positions and indices
/// unchanged, so no second copy is made only to take its channels.
fn conform(mesh: &TriMesh, template: &[AttributeChannel]) -> TriMesh {
    let mut out = crate::grouping::fuse_with_channels(&[mesh], template);
    out.normals.clone_from(&mesh.normals);
    out
}

/// Whether [`conform`] would return `mesh` unchanged: neither it nor the
/// template carries a channel.
fn conforms(mesh: &TriMesh, template: &[AttributeChannel]) -> bool {
    template.is_empty() && mesh.attributes.is_empty()
}

/// A `union_many` tree node: a solid, borrowed from the input until a union
/// produces it, and the fates of the path that built it.
type Node<'a> = (Cow<'a, TriMesh>, Fates);

/// `union_many`'s bookkeeping per input solid, by the most it can reach:
/// its node in the level and in the next, and its share of a level's list of
/// pairs.
const UNION_BOOKKEEPING_PER_SOLID: usize =
    2 * std::mem::size_of::<Node<'static>>() + std::mem::size_of::<&[Node<'static>]>();

/// `subtract_many`'s bookkeeping per tool, by the most it can reach: its
/// bounds, twice (the list and its group's), its index in a group and in the
/// step's member list, and the headers of a group of its own, doubled for
/// lists that grow by pushing.
const SUBTRACT_BOOKKEEPING_PER_TOOL: usize = 2
    * (2 * std::mem::size_of::<axiolid_core::Aabb>()
        + std::mem::size_of::<usize>()
        + std::mem::size_of::<&TriMesh>()
        + 2 * std::mem::size_of::<Vec<usize>>());

/// Heap bytes a tree node holds beyond the input.
fn node_bytes(node: &Node<'_>) -> usize {
    let mesh = match &node.0 {
        Cow::Borrowed(_) => 0,
        Cow::Owned(mesh) => crate::budget::heap_bytes(mesh),
    };
    mesh + crate::budget::fates_bytes(&node.1)
}

/// Remove every channel whose composed fate is `Dropped`.
fn strip_dropped(mesh: &mut TriMesh, fates: &[(String, AttributeFate)]) {
    mesh.attributes.retain(|c| {
        !fates
            .iter()
            .any(|(n, f)| n == &c.name && matches!(f, AttributeFate::Dropped(_)))
    });
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
        // The tool bounds and the grouping are held through the batch;
        // refused before they are built when even they do not fit.
        let bookkeeping = tools.len() * SUBTRACT_BOOKKEEPING_PER_TOOL;
        crate::budget::charge_step(options, bookkeeping, 0, 0)?;
        let bounds: Vec<_> = tools.iter().map(TriMesh::bounds).collect();
        let groups = crate::grouping::disjoint_groups(&bounds);

        // Borrowed until the first group cuts it: a copy of the subject
        // would be scratch held through the whole batch.
        let mut current = Cow::Borrowed(subject);
        let mut sub_operations = 0;
        // Seeded Preserved: the subject enters untouched, and each step's
        // fates compose onto this (`merge_fates`).
        let mut fates = seed_fates(subject);
        for group in &groups {
            // The only real poll point: between groups. Cancelling here returns
            // no mesh at all rather than a partially cut one.
            options.check_cancelled()?;
            sub_operations += 1;
            let members: Vec<&TriMesh> = group.iter().map(|&i| &tools[i]).collect();
            // Charged now, not up front: the subject grows with every group,
            // so this step's operands did not exist before the last one
            // (#259). A single-member group gains nothing from fusing and
            // holds no fused copy.
            let fused_bytes = match members.as_slice() {
                [_] => 0,
                _ => crate::budget::fused_bound(&members, &subject.attributes),
            };
            let held = match &current {
                Cow::Borrowed(_) => 0,
                Cow::Owned(mesh) => crate::budget::heap_bytes(mesh),
            } + crate::budget::fates_bytes(&fates)
                + fused_bytes
                + bookkeeping;
            let elements = current.triangle_count()
                + members.iter().map(|m| m.triangle_count()).sum::<usize>();
            crate::budget::charge_step(options, held, 1, elements)?;
            // The boolean gets what the batch does not hold, and charges its
            // own result against that (#308).
            let budget = crate::budget::remaining(options, held);
            let step = if let [only] = members.as_slice() {
                self.difference(&current, only, options, budget)?
            } else {
                let fused = crate::grouping::fuse_with_channels(&members, &subject.attributes);
                self.difference(&current, &fused, options, budget)?
            };
            fates = merge_fates(&fates, step.evidence.attribute_fates);
            current = Cow::Owned(step.mesh);
        }
        let current = current.into_owned();
        let borrowed: Vec<&TriMesh> = tools.iter().collect();
        let evidence =
            evidence_for(subject, &borrowed, &current, sub_operations).with_attribute_fates(fates);
        Ok(BooleanOutcome::new(current, evidence))
    }

    /// Union every solid by balanced pairwise reduction.
    ///
    /// The sequential fold the trait default performs makes step `i` union
    /// an accumulator already holding `i` solids against one more, so the
    /// subject is re-walked on every step and total work is quadratic in the
    /// operand count. Reducing in pairs -- union neighbours, then pairs of
    /// those, until one remains -- issues the SAME number of booleans, but
    /// the operands stay small until the final levels.
    ///
    /// Measured on a k^3 grid of icospheres (see the `benchmarks` sibling
    /// repo, `sphere_grid`): 1.4x at 8 solids, 7.8x at 125, 28.9x at 512
    /// (44.4s to 1.5s). The gap widens with n because the difference is
    /// complexity, not a constant factor.
    ///
    /// # Known cost: intermediates are rebuilt at every level
    ///
    /// Each level hands the next one a [`TriMesh`], so every intermediate is
    /// converted back from the internal half-edge form and rebuilt by
    /// `to_manifold` on the level above. Instrumenting `to_manifold` against
    /// `compute_boolean` on a disjoint sphere grid:
    ///
    /// ```text
    ///   n   wall_ms  build_ms  build%  triangles rebuilt vs input
    ///   8       8.9       3.9   43.8%   3.00x
    ///  27      51.2      22.2   43.4%   4.85x
    ///  64     145.7      65.2   44.7%   6.00x
    /// 125     344.1     147.5   42.9%   6.98x
    /// ```
    ///
    /// Construction is ~43% of wall time and the redundancy grows with `n`:
    /// at 125 solids the same triangles are rebuilt seven times. Threading
    /// the built structure through the levels instead would have a ceiling
    /// of about **1.6x** for this path (perfect reuse, zero conversion
    /// cost), so the real figure would be lower.
    ///
    /// It has NOT been done, deliberately. It requires a new internal type
    /// boundary carrying a collider and planar grid between levels, which
    /// makes intermediates heavier in memory, and it touches the code path
    /// with the recorded ordering nondeterminism -- while this path is
    /// currently STABLE at 1000 solids. That is a real property to risk for
    /// a bounded win.
    ///
    /// Note the ceiling applies to BATCH paths only. A single
    /// [`MeshBoolean::boolean`] call has no intermediates -- both operands
    /// are leaves -- so construction there is irreducible and reuse saves
    /// exactly nothing. It is not a lead on the single-boolean gap against
    /// upstream Manifold, which profiling showed to be flat rather than
    /// hotspot-shaped.
    ///
    /// # Memory budget (#259)
    ///
    /// Each level is charged before its pairs start, on what then exists:
    /// the nodes it holds, and the declared bound of every pair that can be
    /// live. With `parallel-batch` that is all of them at once, since a
    /// worker waiting in one boolean's join can steal another pair; without
    /// it, one pair at a time on top of the results so far. Each union gets
    /// the budget the level does not hold, shared evenly between pairs that
    /// run at once, and charges its result against it (#308). See
    /// `budget.rs`.
    ///
    /// Union is associative and commutative, so any reduction order yields
    /// the same solid; unlike `subtract_grouped` this needs no disjointness
    /// precondition and no fusing, and therefore has no correctness cliff.
    /// Floating-point differences remain -- a differently ordered triangle
    /// list sums to a marginally different volume -- so the gates compare
    /// volumes to a relative tolerance, not bitwise.
    pub(crate) fn union_tree(
        &self,
        solids: &[TriMesh],
        options: &ExecutionOptions,
    ) -> GeomResult<BooleanOutcome> {
        // The union of nothing is nothing. Returning an empty solid keeps
        // this total rather than making every caller special-case it.
        if solids.is_empty() {
            let borrowed: Vec<&TriMesh> = Vec::new();
            let empty = TriMesh::default();
            let evidence = evidence_for(&empty, &borrowed, &empty, 0);
            return Ok(BooleanOutcome::new(empty, evidence));
        }

        // Every solid carries solids[0]'s channel set, so a channel the
        // result keeps is defined (or explicitly unmapped) on every piece.
        // Each node carries the fates of the path that built it. A solid that
        // already has that set (none, the common case) is borrowed: a copy
        // would be scratch held through the first level.
        let template = &solids[0].attributes;
        let seed = seed_fates(&solids[0]);
        let copies: usize = solids
            .iter()
            .filter(|s| !conforms(s, template))
            .map(|s| crate::budget::fused_bound(&[s], template))
            .sum();
        // The node lists of a level and the next, held through the batch.
        let bookkeeping = solids.len() * UNION_BOOKKEEPING_PER_SOLID;
        crate::budget::charge_step(
            options,
            copies + solids.len() * crate::budget::fates_bytes(&seed) + bookkeeping,
            0,
            0,
        )?;
        let mut level: Vec<Node<'_>> = solids
            .iter()
            .map(|s| {
                let mesh = if conforms(s, template) {
                    Cow::Borrowed(s)
                } else {
                    Cow::Owned(conform(s, template))
                };
                (mesh, seed.clone())
            })
            .collect();
        let mut sub_operations = 0;
        while level.len() > 1 {
            // An odd trailing solid rides to the next level uncombined
            // rather than being unioned into an already-merged neighbour,
            // which would unbalance the tree this exists to keep balanced.
            // Moved out, not copied, so it is held once.
            let carried = if level.len() % 2 == 1 {
                level.pop()
            } else {
                None
            };
            // What the level holds while its pairs run: its own copies and
            // the previous level's results, with every node's history.
            let held: usize = level
                .iter()
                .chain(carried.iter())
                .map(node_bytes)
                .sum::<usize>()
                + bookkeeping;
            // Every pair at one level is independent: no pair reads another
            // pair's output, so the level is a pure map. That is the whole
            // reason this axis can be threaded while intra-solve threading
            // cannot pay for itself -- there is no coordination inside the
            // level, only a join at the end of it.
            let pairs = level.as_chunks::<2>().0;
            let step = |pair: &[Node<'_>], budget: Option<usize>| -> GeomResult<(TriMesh, Fates)> {
                let out = self.union(&pair[0].0, &pair[1].0, options, budget)?;
                // Both halves' histories, then this union.
                let history = merge_fates(&pair[0].1, pair[1].1.clone());
                Ok((
                    out.mesh,
                    merge_fates(&history, out.evidence.attribute_fates),
                ))
            };
            let triangles =
                |pair: &[Node<'_>]| pair[0].0.triangle_count() + pair[1].0.triangle_count();

            // Cancellation is polled ONCE per level rather than per pair.
            // Under rayon the pairs are not ordered, so a per-pair poll
            // would abort at a schedule-dependent point; per level, the
            // cut point is deterministic and the partial work discarded is
            // the same either way.
            options.check_cancelled()?;

            #[cfg(feature = "parallel-batch")]
            let mut next: Vec<Node<'_>> = {
                use rayon::prelude::*;
                // Charged for every pair at once (#259): the pairs run
                // concurrently, and a worker waiting inside one boolean's
                // join can steal another pair and run it on top, so up to
                // every pair of the level is in flight however few workers
                // there are. A finished pair's result is covered by its own
                // term.
                let chunks: Vec<&[Node<'_>]> =
                    pairs.iter().map(<[Node<'_>; 2]>::as_slice).collect();
                let elements = chunks.iter().map(|pair| triangles(pair)).sum();
                crate::budget::charge_step(options, held, chunks.len(), elements)?;
                // Each pair charges its result against its own share of the
                // budget (#308), so pairs running at once cannot together
                // overrun it.
                let share = |pair: &[Node<'_>]| {
                    crate::budget::concurrent_share(
                        options,
                        held,
                        chunks.len(),
                        elements,
                        triangles(pair),
                    )
                };
                // `par_iter` over an indexed slice, NOT `par_bridge`: the
                // bridge does not preserve order, which would permute the
                // level and silently change the tree's shape from run to
                // run. An indexed parallel iterator collects positionally,
                // so the output is identical to the sequential path.
                let merged: Result<Vec<Node<'_>>, GeomError> = chunks
                    .par_iter()
                    .map(|pair| {
                        step(pair, share(pair)).map(|(mesh, fates)| (Cow::Owned(mesh), fates))
                    })
                    .collect();
                merged?
            };

            #[cfg(not(feature = "parallel-batch"))]
            let mut next: Vec<Node<'_>> = {
                let mut acc: Vec<Node<'_>> = Vec::with_capacity(level.len() / 2 + 1);
                for pair in pairs.iter().map(<[Node<'_>; 2]>::as_slice) {
                    // One boolean at a time, on top of the level and the
                    // results so far.
                    let done: usize = acc.iter().map(node_bytes).sum();
                    crate::budget::charge_step(options, held + done, 1, triangles(pair))?;
                    // The union gets what the batch does not hold, and
                    // charges its own result against that (#308).
                    let budget = crate::budget::remaining(options, held + done);
                    let (mesh, fates) = step(pair, budget)?;
                    acc.push((Cow::Owned(mesh), fates));
                }
                acc
            };

            sub_operations += next.len();
            next.extend(carried);
            level = next;
        }

        let (result, fates) = level.into_iter().next().expect("non-empty input");
        let mut result = result.into_owned();
        // A channel some step dropped is not on the result, or only partly:
        // strip it so the mesh never carries a channel its fate calls lost.
        strip_dropped(&mut result, &fates);
        let borrowed: Vec<&TriMesh> = solids.iter().collect();
        // Evidence names the first operand as the subject and the rest as
        // tools, matching how a caller reads a batch: one solid grown by the
        // others. The reduction order is an implementation detail.
        let (subject, tools) = borrowed.split_first().expect("non-empty input");
        let evidence =
            evidence_for(subject, tools, &result, sub_operations).with_attribute_fates(fates);
        Ok(BooleanOutcome::new(result, evidence))
    }

    /// One union, routed through the same validation as `boolean`, with
    /// `budget` bytes for it alone.
    fn union(
        &self,
        subject: &TriMesh,
        tool: &TriMesh,
        options: &ExecutionOptions,
        budget: Option<usize>,
    ) -> GeomResult<BooleanOutcome> {
        self.boolean_within(
            subject,
            tool,
            BooleanOperator::Union,
            options,
            false,
            budget,
        )
    }

    /// One difference, routed through the same validation as `boolean`,
    /// with `budget` bytes for it alone.
    fn difference(
        &self,
        subject: &TriMesh,
        tool: &TriMesh,
        options: &ExecutionOptions,
        budget: Option<usize>,
    ) -> GeomResult<BooleanOutcome> {
        self.boolean_within(
            subject,
            tool,
            BooleanOperator::Difference,
            options,
            false,
            budget,
        )
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

    /// A refusal inside the solve blames this provider, not the operands.
    ///
    /// Replaces the grid-union reproducer that reached the refusal until
    /// #203 was fixed; no admissible input is known to reach it now, so the
    /// mapping is pinned at the seam where `compute_boolean`'s error lands.
    #[test]
    fn a_refusal_inside_the_solve_is_a_backend_contract_violation() {
        let error = solve_refusal(
            BooleanOperator::Union,
            "boolean produced 3 edge points, which cannot be paired into halfedges",
        );
        assert!(
            matches!(
                &error,
                GeomError::BackendContractViolation { backend, detail }
                    if *backend == BoolmeshBoolean::ID
                        && detail.contains("Union failed inside the solve")
                        && detail.contains("cannot be paired")
            ),
            "{error:?}"
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
        for corner in mesh.indices.as_chunks_mut::<3>().0 {
            corner.swap(1, 2);
        }
        assert!(check_result(&mesh, BooleanOperator::Difference).is_ok());
    }
}
