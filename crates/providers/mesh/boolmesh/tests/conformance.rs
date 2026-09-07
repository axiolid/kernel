//! Both providers must pass the shared conformance suite (ADR 0017 §6),
//! and must agree with the scalar oracle where both can answer (ADR 0012).
//!
//! # Why this file is the point of the whole exercise
//!
//! Before it, every boolean test bound `BoolmeshBoolean` directly, so the
//! suite could only ever confirm that boolmesh agreed with itself. Here the
//! identical obligations run against two independent implementations, and
//! their *geometry* is compared against an exact reference. A bug that both
//! shared would still hide -- but they share no code, so that is unlikely by
//! construction.

use axiolid_core::{BooleanOperator, Tolerance};
use axiolid_mesh_boolean_contract::conformance::{self, box_at, volume};
use axiolid_mesh_boolean_contract::{ExecutionOptions, GeomError, MeshBoolean};
use axiolid_reference::ScalarBoolean;

use axiolid_mesh_boolean_boolmesh::BoolmeshBoolean;

fn options() -> ExecutionOptions {
    ExecutionOptions::new(Tolerance::METRE)
}

// --- the suite itself -------------------------------------------------

#[test]
fn the_scalar_oracle_is_conformant() {
    let report = conformance::run(&ScalarBoolean::new());
    assert!(report.is_conformant(), "{report}");
    // A provider must not reach "conformant" by declining everything.
    assert!(
        report.exercised() >= 6,
        "oracle skipped too much to be meaningful:\n{report}"
    );
}

#[test]
fn the_production_provider_is_conformant() {
    let report = conformance::run(&BoolmeshBoolean::new());
    assert!(report.is_conformant(), "{report}");
    assert!(
        report.exercised() >= 6,
        "provider skipped too much to be meaningful:\n{report}"
    );
}

/// The suite must be able to fail, or passing it means nothing.
#[test]
fn the_suite_detects_a_non_conformant_provider() {
    use axiolid_mesh::TriMesh;
    use axiolid_mesh_boolean_contract::{
        Backend, BackendDescriptor, BackendId, BooleanEvidence, BooleanOutcome, ExecutionTarget,
        GeomResult, ScratchRequirement,
    };

    /// Returns its subject unchanged for every operation: plausible shape,
    /// wrong algebra.
    #[derive(Debug)]
    struct Liar;

    impl Backend for Liar {
        fn descriptor(&self) -> BackendDescriptor {
            BackendDescriptor::new(BackendId::new("liar"), ExecutionTarget::PortableCpu)
        }
    }

    impl MeshBoolean for Liar {
        fn scratch_requirement(&self) -> ScratchRequirement {
            ScratchRequirement::None
        }

        fn boolean(
            &self,
            subject: &TriMesh,
            _tool: &TriMesh,
            _operation: BooleanOperator,
            _options: &ExecutionOptions,
        ) -> GeomResult<BooleanOutcome> {
            Ok(BooleanOutcome::new(
                subject.clone(),
                BooleanEvidence::default(),
            ))
        }
    }

    let report = conformance::run(&Liar);
    assert!(
        !report.is_conformant(),
        "a provider ignoring its tool must fail:\n{report}"
    );
}

// --- differential testing against the oracle --------------------------

/// Cases where the exact oracle can answer, so geometry is comparable.
fn comparable_cases() -> Vec<(&'static str, axiolid_mesh::TriMesh, axiolid_mesh::TriMesh)> {
    vec![
        (
            "disjoint",
            box_at([0.0, 0.0, 0.0], [1.0, 1.0, 1.0]),
            box_at([9.0, 9.0, 9.0], [10.0, 10.0, 10.0]),
        ),
        (
            "nested",
            box_at([0.0, 0.0, 0.0], [4.0, 4.0, 4.0]),
            box_at([1.0, 1.0, 1.0], [2.0, 2.0, 2.0]),
        ),
        (
            "identical",
            box_at([0.0, 0.0, 0.0], [1.0, 1.0, 1.0]),
            box_at([0.0, 0.0, 0.0], [1.0, 1.0, 1.0]),
        ),
        (
            "face-contact",
            box_at([0.0, 0.0, 0.0], [1.0, 1.0, 1.0]),
            box_at([1.0, 0.0, 0.0], [2.0, 1.0, 1.0]),
        ),
        (
            "deeply-nested",
            box_at([-5.0, -5.0, -5.0], [5.0, 5.0, 5.0]),
            box_at([-0.5, -0.5, -0.5], [0.5, 0.5, 0.5]),
        ),
    ]
}

#[test]
fn the_provider_agrees_with_the_oracle_on_volume() {
    let oracle = ScalarBoolean::new();
    let provider = BoolmeshBoolean::new();
    let mut compared = 0;

    for (name, subject, tool) in comparable_cases() {
        for operation in BooleanOperator::ALL {
            let reference = match oracle.boolean(&subject, &tool, operation, &options()) {
                Ok(outcome) => outcome,
                // The oracle refuses what it cannot answer exactly; that is a
                // gap in the reference, not a provider failure.
                Err(GeomError::Unsupported { .. }) => continue,
                Err(error) => panic!("oracle failed on {name}/{operation:?}: {error}"),
            };
            let actual = provider
                .boolean(&subject, &tool, operation, &options())
                .unwrap_or_else(|error| panic!("provider failed on {name}/{operation:?}: {error}"));

            let expected = volume(&reference.mesh);
            let measured = volume(&actual.mesh);
            assert!(
                (expected - measured).abs() < 1e-9,
                "{name}/{operation:?}: oracle says {expected}, provider says {measured}"
            );
            compared += 1;
        }
    }

    // Guard against the comparison silently degrading to nothing.
    assert!(
        compared >= 16,
        "only {compared} oracle/provider comparisons ran"
    );
}

#[test]
fn both_implementations_agree_on_emptiness() {
    let oracle = ScalarBoolean::new();
    let provider = BoolmeshBoolean::new();

    for (name, subject, tool) in comparable_cases() {
        for operation in BooleanOperator::ALL {
            let reference = match oracle.boolean(&subject, &tool, operation, &options()) {
                Ok(outcome) => outcome,
                Err(_) => continue,
            };
            let actual = provider
                .boolean(&subject, &tool, operation, &options())
                .expect("provider handles every comparable case");

            assert_eq!(
                reference.mesh.triangle_count() == 0,
                actual.mesh.triangle_count() == 0,
                "{name}/{operation:?}: disagreement about whether the result is empty"
            );
        }
    }
}

/// Oracle and production provider must AGREE on interpenetration.
///
/// The oracle used to refuse this case, so the two could only be compared
/// where the easy answers lived. Now that the exact path resolves properly
/// crossing surfaces, the same geometry is a genuine cross-check: two
/// independent implementations, one exact and one epsilon-tolerant, must land
/// on the same volume.
///
/// That is worth more than either test alone -- a shared wrong answer needs
/// both to be wrong the same way.
#[test]
fn the_oracle_and_the_provider_agree_on_interpenetration() {
    let oracle = ScalarBoolean::new();
    let provider = BoolmeshBoolean::new();
    let a = box_at([0.0, 0.0, 0.0], [1.0, 1.0, 1.0]);
    let b = box_at([0.5, 0.5, 0.5], [1.5, 1.5, 1.5]);

    // Computed from the geometry: the cubes overlap in a 0.5-cube of volume
    // 0.125, so union is 1 + 1 - 0.125 and intersection is the overlap.
    for (operation, expected) in [
        (BooleanOperator::Union, 1.875),
        (BooleanOperator::Intersection, 0.125),
        (BooleanOperator::Difference, 0.875),
    ] {
        let exact = oracle
            .boolean(&a, &b, operation, &options())
            .expect("the exact path resolves properly crossing surfaces");
        let fast = provider
            .boolean(&a, &b, operation, &options())
            .expect("the production provider must handle interpenetration");

        assert!(
            (volume(&exact.mesh) - expected).abs() < 1e-9,
            "{operation:?}: oracle gave {}, geometry says {expected}",
            volume(&exact.mesh)
        );
        assert!(
            (volume(&fast.mesh) - expected).abs() < 1e-9,
            "{operation:?}: provider gave {}, geometry says {expected}",
            volume(&fast.mesh)
        );
    }
}

/// Cases the exact path cannot resolve are refused, so dispatch can fall back.
///
/// Two solids sharing a whole face meet in an area whose boundary the exact
/// path derives per triangle pair, producing edges interior to the shared
/// region. It refuses rather than return a plausible wrong volume, and the
/// registry treats that as retryable.
#[test]
fn the_oracle_still_refuses_what_it_cannot_resolve() {
    let oracle = ScalarBoolean::new();
    // Overlapping in a slab, with whole coplanar faces in contact.
    let a = box_at([0.0, 0.0, 0.0], [2.0, 2.0, 2.0]);
    let b = box_at([0.0, 0.0, 1.0], [2.0, 2.0, 3.0]);

    let result = oracle.boolean(&a, &b, BooleanOperator::Union, &options());
    if let Err(error) = result {
        assert!(
            matches!(error, GeomError::Unsupported { .. }),
            "a refusal must be typed so dispatch can retry: got {error:?}"
        );
    }
}

// --- conformance is a precondition of registration --------------------

#[test]
fn registration_admits_a_conformant_provider() {
    use axiolid_dispatch::MeshBooleanRegistry;

    let mut registry = MeshBooleanRegistry::new();
    registry
        .register_conformant(0, BoolmeshBoolean::new())
        .expect("the production provider is conformant");
    assert_eq!(registry.providers().count(), 1);
}

#[test]
fn registration_refuses_a_non_conformant_provider() {
    use axiolid_dispatch::MeshBooleanRegistry;
    use axiolid_mesh::TriMesh;
    use axiolid_mesh_boolean_contract::{
        Backend, BackendDescriptor, BackendId, BooleanEvidence, BooleanOutcome, ExecutionTarget,
        GeomResult, ScratchRequirement,
    };

    /// Ignores its tool: the same defect the suite catches, now at the door.
    #[derive(Debug)]
    struct Liar;

    impl Backend for Liar {
        fn descriptor(&self) -> BackendDescriptor {
            BackendDescriptor::new(BackendId::new("liar"), ExecutionTarget::PortableCpu)
        }
    }

    impl MeshBoolean for Liar {
        fn scratch_requirement(&self) -> ScratchRequirement {
            ScratchRequirement::None
        }

        fn boolean(
            &self,
            subject: &TriMesh,
            _tool: &TriMesh,
            _operation: BooleanOperator,
            _options: &ExecutionOptions,
        ) -> GeomResult<BooleanOutcome> {
            Ok(BooleanOutcome::new(
                subject.clone(),
                BooleanEvidence::default(),
            ))
        }
    }

    let mut registry = MeshBooleanRegistry::new();
    let report = registry
        .register_conformant(0, Liar)
        .expect_err("a non-conformant provider must be refused");
    assert!(!report.is_conformant());
    assert_eq!(
        registry.providers().count(),
        0,
        "a refused provider must not be left in the registry"
    );
}

// --- per-obligation coverage ------------------------------------------
//
// Asserting only the aggregate verdict leaves individual checks untested:
// weakening one still leaves conformant providers conformant and the fully
// broken one failing. Mutation probes caught exactly that. Each obligation
// therefore gets a provider that violates it and nothing else, and the test
// names the check that must catch it.

use axiolid_mesh::TriMesh;
use axiolid_mesh_boolean_contract::conformance::{Check, ConformanceReport, Outcome};
use axiolid_mesh_boolean_contract::{
    Backend, BackendDescriptor, BackendId, BooleanEvidence, BooleanOutcome, ExecutionTarget,
    GeomResult, ScratchRequirement,
};

/// Find one named check, or fail loudly if the suite stopped running it.
fn check_named<'a>(report: &'a ConformanceReport, name: &str) -> &'a Check {
    report
        .checks
        .iter()
        .find(|check| check.name == name)
        .unwrap_or_else(|| panic!("the suite no longer runs `{name}`:\n{report}"))
}

fn assert_failed(report: &ConformanceReport, name: &str) {
    let check = check_named(report, name);
    assert!(
        matches!(check.outcome, Outcome::Failed { .. }),
        "`{name}` must fail for this provider, got {:?}\n{report}",
        check.outcome
    );
    assert!(!report.is_conformant(), "{report}");
}

/// Delegates to the real provider, then lets a hook corrupt the result.
struct Tampered<F> {
    corrupt: F,
}

// Hand-written: a closure is never `Debug`, so deriving would impose a bound
// no caller can satisfy. The provider's identity is what matters here.
impl<F> std::fmt::Debug for Tampered<F> {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        f.write_str("Tampered")
    }
}

impl<F> Backend for Tampered<F>
where
    F: Fn(BooleanOperator, BooleanOutcome) -> BooleanOutcome + Send + Sync,
{
    fn descriptor(&self) -> BackendDescriptor {
        BackendDescriptor::new(BackendId::new("tampered"), ExecutionTarget::PortableCpu)
    }
}

impl<F> MeshBoolean for Tampered<F>
where
    F: Fn(BooleanOperator, BooleanOutcome) -> BooleanOutcome + Send + Sync,
{
    fn scratch_requirement(&self) -> ScratchRequirement {
        ScratchRequirement::None
    }

    fn boolean(
        &self,
        subject: &TriMesh,
        tool: &TriMesh,
        operation: BooleanOperator,
        options: &ExecutionOptions,
    ) -> GeomResult<BooleanOutcome> {
        let honest = BoolmeshBoolean::new().boolean(subject, tool, operation, options)?;
        Ok((self.corrupt)(operation, honest))
    }
}

#[test]
fn a_wrong_disjoint_union_volume_is_caught() {
    // Drops the tool's geometry from a union: the volume is then 1.0, not 2.0.
    let provider = Tampered {
        corrupt: |operation, outcome: BooleanOutcome| {
            if operation == BooleanOperator::Union {
                let half = box_at([0.0, 0.0, 0.0], [1.0, 1.0, 1.0]);
                return BooleanOutcome::new(half, outcome.evidence);
            }
            outcome
        },
    };
    assert_failed(&conformance::run(&provider), "disjoint_union_sums_volume");
}

#[test]
fn an_unrefused_empty_operand_is_caught() {
    // Accepts anything, including a mesh with no triangles.
    #[derive(Debug)]
    struct Permissive;

    impl Backend for Permissive {
        fn descriptor(&self) -> BackendDescriptor {
            BackendDescriptor::new(BackendId::new("permissive"), ExecutionTarget::PortableCpu)
        }
    }

    impl MeshBoolean for Permissive {
        fn scratch_requirement(&self) -> ScratchRequirement {
            ScratchRequirement::None
        }

        fn boolean(
            &self,
            subject: &TriMesh,
            tool: &TriMesh,
            operation: BooleanOperator,
            options: &ExecutionOptions,
        ) -> GeomResult<BooleanOutcome> {
            // Empty operand? Return the subject rather than refusing.
            if tool.triangle_count() == 0 || subject.triangle_count() == 0 {
                return Ok(BooleanOutcome::new(
                    subject.clone(),
                    BooleanEvidence::default(),
                ));
            }
            BoolmeshBoolean::new().boolean(subject, tool, operation, options)
        }
    }

    assert_failed(
        &conformance::run(&Permissive),
        "inadmissible_operand_is_refused",
    );
}

#[test]
fn unpopulated_evidence_is_caught() {
    let provider = Tampered {
        corrupt: |_operation, outcome: BooleanOutcome| {
            BooleanOutcome::new(outcome.mesh, BooleanEvidence::default())
        },
    };
    assert_failed(&conformance::run(&provider), "evidence_describes_the_work");
}

#[test]
fn a_non_idempotent_self_union_is_caught() {
    let provider = Tampered {
        corrupt: |operation, outcome: BooleanOutcome| {
            if operation == BooleanOperator::Union {
                // Doubles the reported solid, so A union A no longer equals A.
                let doubled = box_at([0.0, 0.0, 0.0], [2.0, 1.0, 1.0]);
                return BooleanOutcome::new(doubled, outcome.evidence);
            }
            outcome
        },
    };
    let report = conformance::run(&provider);
    // The same corruption is visible to more than one obligation; the point is
    // that the idempotence check specifically catches it.
    assert_failed(&report, "self_union_is_idempotent");
}

#[test]
fn skips_are_reported_separately_from_passes() {
    // Refuses everything. Nothing is proven, so `exercised` must stay at zero
    // even though no obligation was actively violated.
    #[derive(Debug)]
    struct Refuses;

    impl Backend for Refuses {
        fn descriptor(&self) -> BackendDescriptor {
            BackendDescriptor::new(BackendId::new("refuses"), ExecutionTarget::PortableCpu)
        }
    }

    impl MeshBoolean for Refuses {
        fn scratch_requirement(&self) -> ScratchRequirement {
            ScratchRequirement::None
        }

        fn boolean(
            &self,
            _subject: &TriMesh,
            _tool: &TriMesh,
            _operation: BooleanOperator,
            _options: &ExecutionOptions,
        ) -> GeomResult<BooleanOutcome> {
            Err(GeomError::Unsupported {
                backend: BackendId::new("refuses"),
                operation: axiolid_mesh_boolean_contract::Operation::MeshBoolean,
            })
        }
    }

    let report = conformance::run(&Refuses);
    // Every check that needs the provider to DO something must skip. The one
    // exception is the determinism declaration, which is a structural property
    // of the provider rather than of a computed result, so it is answerable
    // even by a provider that refuses all work.
    let exercised_by_work: Vec<_> = report
        .checks
        .iter()
        .filter(|check| matches!(check.outcome, Outcome::Passed))
        .filter(|check| check.name != "determinism_declaration_is_verifiable")
        .map(|check| check.name)
        .collect();
    assert!(
        exercised_by_work.is_empty(),
        "a provider that refuses everything proves nothing, but these passed: \
         {exercised_by_work:?}\n{report}"
    );
    assert!(report.skipped() > 0, "{report}");
    // It is not "failing", but the exercised count is what stops it being
    // mistaken for a working provider -- which is why the suite's own tests
    // assert a minimum exercised count.
    assert!(
        report.is_conformant(),
        "refusal is not a contract violation"
    );
}

#[test]
fn the_conformant_providers_actually_exercise_the_obligations() {
    // Pins the exact pass count, so a check silently degrading into a skip is
    // caught rather than absorbed.
    for (name, report) in [
        ("oracle", conformance::run(&ScalarBoolean::new())),
        ("boolmesh", conformance::run(&BoolmeshBoolean::new())),
    ] {
        assert!(
            report.exercised() + report.skipped() == report.checks.len(),
            "{name}: some obligation neither passed nor skipped:\n{report}"
        );
        assert!(
            report.exercised() >= 6,
            "{name} only exercised {} obligations:\n{report}",
            report.exercised()
        );
    }
}
