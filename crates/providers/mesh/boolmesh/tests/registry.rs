//! The provider must work through `MeshBooleanRegistry`, not just directly.
//!
//! This is what makes the seam real: `ifc-geometry` sees the registry, never
//! this crate.

mod support;

use axiolid_contracts::ExecutionOptions;
use axiolid_core::{BooleanOperator, Tolerance};
use axiolid_dispatch::MeshBooleanRegistry;
use axiolid_mesh_boolean_boolmesh::BoolmeshBoolean;
use support::{boxx, volume};

fn registry() -> MeshBooleanRegistry {
    let mut registry = MeshBooleanRegistry::new();
    registry.register(0, BoolmeshBoolean::new());
    registry
}

/// Dispatching through the registry produces the same geometry as calling the
/// provider directly.
#[test]
fn the_registry_dispatches_to_the_provider() {
    let wall = boxx(2.0, 0.1, 0.0, 4.0, 0.2, 3.0, 0.0);
    let opening = boxx(1.5, 0.1, 0.3, 1.0, 0.4, 1.2, 0.0);

    let result = registry()
        .boolean(
            &wall,
            &opening,
            BooleanOperator::Difference,
            &ExecutionOptions::new(Tolerance::METRE),
        )
        .expect("registry dispatch")
        .mesh;

    assert!(volume(&result) < volume(&wall));
}

/// A caller with a hard memory budget must be refused rather than allowed to
/// allocate past it: this provider declares `Unbounded` scratch.
#[test]
fn a_bounded_memory_budget_refuses_this_provider() {
    let options = ExecutionOptions::new(Tolerance::METRE).with_memory_budget(1024);
    let error = registry()
        .boolean(
            &boxx(2.0, 0.1, 0.0, 4.0, 0.2, 3.0, 0.0),
            &boxx(1.5, 0.1, 0.3, 1.0, 0.4, 1.2, 0.0),
            BooleanOperator::Difference,
            &options,
        )
        .expect_err("an unbounded provider cannot fit a declared budget");

    assert!(
        matches!(error, axiolid_contracts::GeomError::BudgetExceeded { .. }),
        "expected BudgetExceeded, got {error:?}"
    );
}

/// The per-worker scratch term is charged for the pool the boolean actually
/// runs on, not the width the budget check before dispatch assumed (#226).
///
/// The budget fits the declared bound for one worker exactly. Asked for
/// `Serial`, the registry admits the call; run inside a two-worker pool the
/// provider must refuse, because there the bound is one worker's term larger.
#[cfg(feature = "parallel")]
#[test]
fn a_pool_wider_than_the_budget_allows_is_refused() {
    use axiolid_contracts::{Parallelism, ScratchRequirement};
    use axiolid_mesh_boolean_contract::MeshBoolean;

    let subject = boxx(2.0, 0.1, 0.0, 4.0, 0.2, 3.0, 0.0);
    let tool = boxx(1.5, 0.1, 0.3, 1.0, 0.4, 1.2, 0.0);
    let elements = subject.triangle_count() + tool.triangle_count();
    let requirement = BoolmeshBoolean::new().scratch_requirement();
    assert!(matches!(
        requirement,
        ScratchRequirement::Affine { bytes_per_worker, .. } if bytes_per_worker > 0
    ));
    let one_worker = requirement
        .upper_bound_bytes_on(elements, 1)
        .expect("finite bound");
    let options = ExecutionOptions::new(Tolerance::METRE)
        .with_memory_budget(one_worker)
        .with_parallelism(Parallelism::Serial)
        .expect("serial");
    let run_on = |workers| {
        rayon::ThreadPoolBuilder::new()
            .num_threads(workers)
            .build()
            .expect("pool")
            .install(|| registry().boolean(&subject, &tool, BooleanOperator::Difference, &options))
    };

    run_on(1).expect("one worker fits the budget exactly");
    let error = run_on(2).expect_err("two workers do not fit");
    assert!(
        matches!(error, axiolid_contracts::GeomError::BudgetExceeded { .. }),
        "expected BudgetExceeded, got {error:?}"
    );
}
