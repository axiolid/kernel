//! Contract gates for the execution policy surface.

use axiolid_contracts::{
    BackendId, DataResidency, Determinism, ExecutionOptions, GeomError, Parallelism, Residency,
    ScratchRequirement,
};
use axiolid_core::Tolerance;

fn options() -> ExecutionOptions {
    ExecutionOptions::new(Tolerance::METRE)
}

// --- R9: determinism is three contracts, not one flag -------------------

#[test]
fn determinism_levels_form_a_strength_lattice() {
    use Determinism::{BestEffort, Bitwise, NumericallyBounded, Topological};

    // Every level satisfies itself and everything weaker.
    for (stronger, weaker) in [
        (Topological, BestEffort),
        (NumericallyBounded, Topological),
        (Bitwise, NumericallyBounded),
        (Bitwise, BestEffort),
    ] {
        assert!(
            stronger.satisfies(weaker),
            "{stronger:?} must satisfy {weaker:?}"
        );
        assert!(
            !weaker.satisfies(stronger),
            "{weaker:?} must NOT satisfy {stronger:?}"
        );
    }
    for level in [BestEffort, Topological, NumericallyBounded, Bitwise] {
        assert!(level.satisfies(level), "{level:?} must satisfy itself");
    }
}

/// The whole point of splitting the enum: stable ordering and bitwise identity
/// are different promises. A backend offering only the former must be rejected
/// for the latter, which a single `Required` flag could not express.
#[test]
fn topological_determinism_does_not_promise_bitwise_identity() {
    assert!(!Determinism::Topological.satisfies(Determinism::Bitwise));
    assert!(!Determinism::NumericallyBounded.satisfies(Determinism::Bitwise));
}

/// Each of the three named contracts must be independently selectable and
/// survive a round trip through the builder.
#[test]
fn every_determinism_contract_is_distinctly_selectable() {
    let levels = [
        Determinism::BestEffort,
        Determinism::Topological,
        Determinism::NumericallyBounded,
        Determinism::Bitwise,
    ];
    for level in levels {
        assert_eq!(
            options().with_determinism(level).determinism(),
            level,
            "{level:?} must round-trip unchanged"
        );
    }
    // Distinct variants, not aliases of one another.
    for (i, a) in levels.iter().enumerate() {
        for (j, b) in levels.iter().enumerate() {
            assert_eq!(i == j, a == b, "{a:?} vs {b:?} must be distinct");
        }
    }
}

#[test]
fn the_default_policy_is_numerically_bounded_not_best_effort() {
    // A geometry kernel defaulting to BestEffort would silently permit
    // unstable topology; defaulting to Bitwise would ban every parallel
    // reduction. NumericallyBounded is the honest middle.
    assert_eq!(options().determinism(), Determinism::NumericallyBounded);
}

// --- R8: residency is part of the execution plan ------------------------

#[test]
fn residency_defaults_to_host_and_is_explicitly_overridable() {
    assert_eq!(options().residency(), DataResidency::HOST);

    let gpu = BackendId::new("test-gpu");
    let plan = DataResidency::new(Residency::Device(gpu), Residency::Host);
    let configured = options().with_residency(plan);

    assert_eq!(configured.residency().input(), Residency::Device(gpu));
    assert_eq!(configured.residency().output(), Residency::Host);
}

#[test]
fn residency_distinguishes_device_local_from_host_readable() {
    let gpu = BackendId::new("gpu-a");
    let other = BackendId::new("gpu-b");

    assert!(Residency::Device(gpu).is_local_to(gpu));
    assert!(!Residency::Device(gpu).is_local_to(other));
    assert!(!Residency::Device(gpu).is_host_readable());

    // Unified memory is readable from both sides without a copy.
    assert!(Residency::Unified(gpu).is_local_to(gpu));
    assert!(Residency::Unified(gpu).is_host_readable());

    // Host data is host-readable but never device-local: reaching it from a
    // device costs a transfer, which is exactly what routing must see.
    assert!(Residency::Host.is_host_readable());
    assert!(!Residency::Host.is_local_to(gpu));
}

#[test]
fn a_transfer_free_plan_requires_both_ends_on_the_device() {
    let gpu = BackendId::new("gpu-a");
    let resident = DataResidency::new(Residency::Device(gpu), Residency::Device(gpu));
    assert!(resident.is_transfer_free_on(gpu));

    // Device-resident input but host-wanted output still costs a readback.
    let readback = DataResidency::new(Residency::Device(gpu), Residency::Host);
    assert!(!readback.is_transfer_free_on(gpu));
    assert!(!DataResidency::HOST.is_transfer_free_on(gpu));
}

// --- R4: scratch is declared and the budget is enforced -----------------

#[test]
fn scratch_requirements_report_bounds() {
    assert_eq!(ScratchRequirement::None.upper_bound_bytes(1_000), Some(0));
    assert_eq!(
        ScratchRequirement::Fixed { bytes: 4_096 }.upper_bound_bytes(1_000),
        Some(4_096)
    );
    assert_eq!(
        ScratchRequirement::PerElement {
            bytes_per_element: 32
        }
        .upper_bound_bytes(100),
        Some(3_200)
    );
    assert_eq!(
        ScratchRequirement::Unbounded.upper_bound_bytes(1),
        None,
        "unbounded scratch must not report a bound"
    );
}

/// A per-element bound must not silently wrap into a small number on overflow;
/// that would turn an enormous requirement into an apparently affordable one.
#[test]
fn per_element_bounds_saturate_to_unbounded_on_overflow() {
    let huge = ScratchRequirement::PerElement {
        bytes_per_element: usize::MAX,
    };
    assert_eq!(huge.upper_bound_bytes(2), None);
    assert!(!huge.fits_budget(&options().with_memory_budget(1_024), 2));
}

#[test]
fn an_unbounded_requirement_never_fits_a_declared_budget() {
    let budgeted = options().with_memory_budget(1_000_000);
    assert!(!ScratchRequirement::Unbounded.fits_budget(&budgeted, 1));

    // With no budget declared, anything is permitted.
    assert!(ScratchRequirement::Unbounded.fits_budget(&options(), 1));
}

#[test]
fn budgets_admit_what_fits_and_reject_what_does_not() {
    let budgeted = options().with_memory_budget(1_000);
    let per_element = ScratchRequirement::PerElement {
        bytes_per_element: 10,
    };

    assert!(per_element.fits_budget(&budgeted, 100), "1000 <= 1000");
    assert!(!per_element.fits_budget(&budgeted, 101), "1010 > 1000");
}

/// Every term of an affine bound counts, the worker term per worker (#226).
#[test]
fn affine_bounds_charge_base_elements_and_workers() {
    let affine = ScratchRequirement::Affine {
        base_bytes: 1_000,
        bytes_per_element: 10,
        bytes_per_worker: 100,
    };
    assert_eq!(affine.upper_bound_bytes_on(0, 0), Some(1_000), "base alone");
    assert_eq!(affine.upper_bound_bytes_on(5, 0), Some(1_050));
    assert_eq!(affine.upper_bound_bytes_on(5, 3), Some(1_350));
    // A worker term cannot be bounded without a worker count.
    assert_eq!(affine.upper_bound_bytes(5), None);
    let unthreaded = ScratchRequirement::Affine {
        base_bytes: 1_000,
        bytes_per_element: 10,
        bytes_per_worker: 0,
    };
    assert_eq!(unthreaded.upper_bound_bytes(5), Some(1_050));
    assert_eq!(unthreaded.upper_bound_bytes_on(5, 64), Some(1_050));
}

#[test]
fn affine_bounds_refuse_on_overflow_in_any_term() {
    for (base_bytes, bytes_per_element, bytes_per_worker) in [
        (usize::MAX, 0, 1),
        (0, usize::MAX, 1),
        (0, 1, usize::MAX),
        (usize::MAX, 1, 0),
    ] {
        let affine = ScratchRequirement::Affine {
            base_bytes,
            bytes_per_element,
            bytes_per_worker,
        };
        assert_eq!(affine.upper_bound_bytes_on(2, 2), None, "{affine:?}");
        assert!(!affine.fits_budget_on(&options().with_memory_budget(usize::MAX), 2, 2));
    }
}

/// The budget check charges the worker term for the requested parallelism,
/// and `fits_budget_on` for the width the caller says the call runs on.
#[test]
fn affine_budgets_follow_the_worker_count() {
    let affine = ScratchRequirement::Affine {
        base_bytes: 1_000,
        bytes_per_element: 10,
        bytes_per_worker: 100,
    };
    let budget = |parallelism| {
        options()
            .with_memory_budget(1_250)
            .with_parallelism(parallelism)
            .expect("valid parallelism")
    };
    // 1000 + 10*5 + 100*2 = 1250.
    assert!(affine.fits_budget(&budget(Parallelism::Threads(2)), 5));
    assert!(!affine.fits_budget(&budget(Parallelism::Threads(3)), 5));
    assert!(affine.fits_budget(&budget(Parallelism::Serial), 5));
    assert!(
        !affine.fits_budget(&budget(Parallelism::Serial), 16),
        "1000+160+100"
    );
    assert!(affine.fits_budget_on(&budget(Parallelism::Serial), 5, 2));
    assert!(!affine.fits_budget_on(&budget(Parallelism::Serial), 5, 3));
    let auto = budget(Parallelism::Auto);
    let width = Parallelism::Auto.worker_bound();
    assert_eq!(
        affine.fits_budget(&auto, 5),
        affine.fits_budget_on(&auto, 5, width)
    );
    // No budget admits anything, whatever the width.
    assert!(affine.fits_budget_on(&options(), usize::MAX, usize::MAX));
}

#[test]
fn worker_bounds_follow_the_parallelism_preference() {
    assert_eq!(Parallelism::Serial.worker_bound(), 1);
    assert_eq!(Parallelism::Threads(7).worker_bound(), 7);
    assert!(Parallelism::Auto.worker_bound() >= 1);
}

#[test]
fn charging_scratch_reports_a_structured_budget_error() {
    let budgeted = options().with_memory_budget(512);
    assert!(budgeted.charge_scratch(512).is_ok());
    assert!(matches!(
        budgeted.charge_scratch(513),
        Err(GeomError::BudgetExceeded { resource: "memory" })
    ));
    // No budget means no ceiling.
    assert!(options().charge_scratch(usize::MAX).is_ok());
}
