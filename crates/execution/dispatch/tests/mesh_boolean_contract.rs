//! Contract tests for the solid boolean surface (ADR 0017 sections 1-4).
//!
//! These replace four retired `csg_deferral` guards. They assert the landed
//! behaviour rather than the absence of it, and they are written against the
//! contract types only -- no provider is named here, so they stay valid when
//! the backend changes.

#![cfg(feature = "mesh-boolean")]

use axiolid_contracts::{
    Backend, BackendDescriptor, BackendId, CancellationGranularity, CancellationToken,
    ExecutionOptions, ExecutionTarget, GeomError, GeomResult, ScratchRequirement,
};
use axiolid_core::{BooleanOperator, Tolerance};
use axiolid_dispatch::MeshBooleanRegistry;
use axiolid_mesh::TriMesh;
use axiolid_mesh_boolean_contract::{BooleanEvidence, BooleanOutcome, MeshBoolean};
use axiolid_mesh_contracts::SolidRequirements;

/// Unit cube at the origin, outward-oriented. The minimal admissible operand.
fn cube() -> TriMesh {
    let positions = vec![
        [0.0, 0.0, 0.0].into(),
        [1.0, 0.0, 0.0].into(),
        [1.0, 1.0, 0.0].into(),
        [0.0, 1.0, 0.0].into(),
        [0.0, 0.0, 1.0].into(),
        [1.0, 0.0, 1.0].into(),
        [1.0, 1.0, 1.0].into(),
        [0.0, 1.0, 1.0].into(),
    ];
    let indices = vec![
        0, 2, 1, 0, 3, 2, // bottom, inward normal is -z
        4, 5, 6, 4, 6, 7, // top
        0, 1, 5, 0, 5, 4, // front
        1, 2, 6, 1, 6, 5, // right
        2, 3, 7, 2, 7, 6, // back
        3, 0, 4, 3, 4, 7, // left
    ];
    TriMesh::new(positions, indices)
}

fn options() -> ExecutionOptions {
    ExecutionOptions::new(Tolerance::METRE)
}

/// A provider that records what it was asked to do.
#[derive(Debug, Default)]
struct Recorder {
    granularity: Option<CancellationGranularity>,
}

impl Backend for Recorder {
    fn descriptor(&self) -> BackendDescriptor {
        BackendDescriptor::new(BackendId::new("recorder"), ExecutionTarget::PortableCpu)
    }
}

impl MeshBoolean for Recorder {
    fn scratch_requirement(&self) -> ScratchRequirement {
        ScratchRequirement::None
    }

    fn cancellation_granularity(&self) -> CancellationGranularity {
        self.granularity.unwrap_or(CancellationGranularity::None)
    }

    fn boolean(
        &self,
        subject: &TriMesh,
        _tool: &TriMesh,
        _operation: BooleanOperator,
        options: &ExecutionOptions,
    ) -> GeomResult<BooleanOutcome> {
        options.check_cancelled()?;
        let evidence = BooleanEvidence::record(subject.triangle_count(), 0, 0, 0);
        Ok(BooleanOutcome::new(subject.clone(), evidence))
    }
}

fn registry_with(provider: Recorder) -> MeshBooleanRegistry {
    let mut registry = MeshBooleanRegistry::new();
    registry.register(0, provider);
    registry
}

// --- Section 1: operation set -----------------------------------------

#[test]
fn the_operand_set_matches_the_planar_contract() {
    // Four regularized operations, the same algebra `axiolid-overlay` exposes.
    assert_eq!(BooleanOperator::ALL.len(), 4);
    assert!(BooleanOperator::ALL.contains(&BooleanOperator::SymmetricDifference));
}

#[test]
fn difference_is_the_only_ordered_operand() {
    assert!(!BooleanOperator::Difference.is_commutative());
    for operation in BooleanOperator::ALL {
        if operation != BooleanOperator::Difference {
            assert!(
                operation.is_commutative(),
                "{operation:?} should be commutative"
            );
        }
    }
}

// --- Section 2: precondition ownership --------------------------------

#[test]
fn the_contract_rejects_inadmissible_operands_before_dispatch() {
    let registry = registry_with(Recorder::default());
    let good = cube();

    // An inside-out operand would silently invert the operation.
    let mut flipped = cube();
    for triangle in flipped.indices.as_chunks_mut::<3>().0 {
        triangle.swap(0, 1);
    }
    let error = registry
        .boolean(&good, &flipped, BooleanOperator::Difference, &options())
        .unwrap_err();
    assert!(
        matches!(error, GeomError::InvalidInput(ref detail) if detail.contains("tool[0]")),
        "rejection must name the offending operand, got {error:?}"
    );

    // An empty operand is malformed data, not geometry.
    let empty = TriMesh::new(Vec::new(), Vec::new());
    assert!(matches!(
        registry.boolean(&good, &empty, BooleanOperator::Union, &options()),
        Err(GeomError::InvalidInput(_))
    ));
}

#[test]
fn admissibility_levels_are_cumulative() {
    let good = cube();
    // A flat sheet is structurally fine but encloses nothing.
    let sheet = TriMesh::new(
        vec![
            [0.0, 0.0, 0.0].into(),
            [1.0, 0.0, 0.0].into(),
            [1.0, 1.0, 0.0].into(),
        ],
        vec![0, 1, 2],
    );

    assert!(SolidRequirements::Structural.validate(&sheet, "s").is_ok());
    assert!(matches!(
        SolidRequirements::Enclosing.validate(&sheet, "s"),
        Err(GeomError::Degenerate(_))
    ));
    assert!(SolidRequirements::Oriented.validate(&good, "s").is_ok());
}

#[test]
fn solid_admission_rejects_positive_infinite_signed_volume() {
    let magnitude = 1.0e308;
    let epsilon = 1.0e-308;
    let overflow = TriMesh::new(
        vec![
            [epsilon, epsilon, 1.0].into(),
            [magnitude, 1.0, epsilon].into(),
            [epsilon, magnitude, 1.0].into(),
        ],
        vec![0, 1, 2],
    );

    let six_volume = overflow.positions[0].dot(overflow.positions[1].cross(overflow.positions[2]));
    assert!(six_volume == f64::INFINITY, "fixture produced {six_volume}");

    let error = SolidRequirements::Oriented
        .validate(&overflow, "subject")
        .unwrap_err();
    assert!(
        matches!(error, GeomError::InvalidInput(ref detail) if detail.contains("non-finite signed volume")),
        "positive-infinite signed volume must fail closed, got {error:?}"
    );
}

// --- Section 3: evidence ----------------------------------------------

#[test]
fn every_boolean_reports_what_it_did() {
    let registry = registry_with(Recorder::default());
    let cube = cube();
    let outcome = registry
        .boolean(&cube, &cube, BooleanOperator::Union, &options())
        .expect("admissible operands");

    assert_eq!(outcome.evidence.subject_triangles, cube.triangle_count());
    assert_eq!(outcome.evidence.sub_operations, 1);
    // The recorder echoes its subject, so the outcome carries that geometry.
    // The point is that mesh and evidence arrive together.
    assert_eq!(outcome.mesh.triangle_count(), cube.triangle_count());
    assert!(!outcome.is_empty());
}

#[test]
fn evidence_absorbs_sub_operations_without_losing_input_counts() {
    let mut total = BooleanEvidence::record(12, 6, 12, 1);
    total.absorb(BooleanEvidence::record(0, 0, 8, 2).with_sub_operations(3));

    // Input counts belong to the outermost call and survive the merge.
    assert_eq!(total.subject_triangles, 12);
    assert_eq!(total.tool_triangles, 6);
    // Output counts come from the final sub-operation; work accumulates.
    assert_eq!(total.output_triangles, 8);
    assert_eq!(total.output_components, 2);
    assert_eq!(total.sub_operations, 4);
}

// --- Section 4: cancellation ------------------------------------------

#[test]
fn a_cancelled_token_stops_the_operation() {
    let registry = registry_with(Recorder::default());
    let cube = cube();
    let token = CancellationToken::new();
    token.cancel();

    let error = registry
        .boolean(
            &cube,
            &cube,
            BooleanOperator::Union,
            &options().with_cancellation(token),
        )
        .unwrap_err();
    assert_eq!(error, GeomError::Cancelled);
}

#[test]
fn cancellation_is_opt_in_and_absent_by_default() {
    assert!(options().cancellation().is_none());
    assert!(options().check_cancelled().is_ok());

    let token = CancellationToken::new();
    let with = options().with_cancellation(token.clone());
    assert!(with.check_cancelled().is_ok(), "not cancelled yet");
    token.cancel();
    assert!(with.check_cancelled().is_err(), "clones share one flag");
}

#[test]
fn token_equality_is_identity_not_current_state() {
    // Two independently created tokens are different cancellation sources even
    // though both are uncancelled right now. Comparing state instead of
    // identity would call them equal and then silently diverge.
    let first = CancellationToken::new();
    let second = CancellationToken::new();
    assert!(!first.is_cancelled() && !second.is_cancelled());
    assert_ne!(first, second, "distinct sources must not compare equal");

    // A clone IS the same source, regardless of state.
    let shared = first.clone();
    assert_eq!(first, shared);
    first.cancel();
    assert_eq!(first, shared, "cancelling must not change identity");
    assert!(shared.is_cancelled(), "the flag is genuinely shared");
}

#[test]
fn a_provider_declares_its_real_polling_granularity() {
    // Silence is the honest default: a provider that has not said it polls is
    // assumed not to.
    assert_eq!(
        Recorder::default().cancellation_granularity(),
        CancellationGranularity::None
    );
    let batching = Recorder {
        granularity: Some(CancellationGranularity::BetweenOperations),
    };
    assert_eq!(
        batching.cancellation_granularity(),
        CancellationGranularity::BetweenOperations
    );
}
