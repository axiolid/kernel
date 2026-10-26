//! Conformance suite every `MeshBoolean` provider must pass (ADR 0017 §6).
//!
//! # Why this is library code, not a test file
//!
//! Before this existed, all provider tests bound the concrete `BoolmeshBoolean`
//! type, so they tested *that provider*, not *the contract*. A second provider
//! inherited zero obligations. This suite is generic over `impl MeshBoolean`
//! and exported, so an out-of-tree provider can run the identical checks.
//!
//! # Usage
//!
//! ```no_run
//! # use axiolid_mesh_boolean_contract::conformance::{self, ConformanceReport};
//! # fn check(provider: &impl axiolid_mesh_boolean_contract::MeshBoolean) {
//! let report = conformance::run(provider);
//! assert!(report.is_conformant(), "{report}");
//! # }
//! ```
//!
//! # What it does and does not prove
//!
//! It checks *contract* obligations: operand algebra, admissibility, evidence,
//! empty-result handling, and determinism. It does not check that geometry is
//! numerically correct -- that is what differential testing against
//! `axiolid-reference`'s oracle is for. A provider passing this suite is
//! well-behaved, not necessarily accurate.
//!
//! # Skips are not passes
//!
//! A provider may legitimately refuse work (`Unsupported`). Such a case is
//! recorded as [`Outcome::Skipped`] with its reason and reported separately,
//! so a provider cannot reach "conformant" by refusing everything: callers can
//! see exactly what was actually exercised.

use core::fmt;

use axiolid_core::{BooleanOperator, Tolerance};

use crate::{BooleanOutcome, MeshBoolean};
use axiolid_contracts::{Determinism, ExecutionOptions, GeomError};
use axiolid_mesh::TriMesh;

/// Result of one conformance check.
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum Outcome {
    /// The provider satisfied the obligation.
    Passed,
    /// The provider declined the input; not a failure, but not proof either.
    Skipped {
        /// Why the provider declined.
        reason: String,
    },
    /// The provider violated the contract.
    Failed {
        /// What went wrong, in terms of the obligation.
        detail: String,
    },
}

/// One named obligation and how the provider answered it.
#[derive(Debug, Clone)]
pub struct Check {
    /// Stable identifier for the obligation.
    pub name: &'static str,
    /// What the provider did.
    pub outcome: Outcome,
}

/// Full conformance result for one provider.
#[derive(Debug, Clone, Default)]
pub struct ConformanceReport {
    /// Every obligation, in execution order.
    pub checks: Vec<Check>,
}

impl ConformanceReport {
    /// Whether the provider violated no obligation.
    ///
    /// Skips do not fail conformance -- a provider is allowed to decline work.
    /// Read [`Self::exercised`] to see how much was actually proven.
    #[must_use]
    pub fn is_conformant(&self) -> bool {
        !self
            .checks
            .iter()
            .any(|check| matches!(check.outcome, Outcome::Failed { .. }))
    }

    /// How many obligations the provider actually satisfied.
    #[must_use]
    pub fn exercised(&self) -> usize {
        self.checks
            .iter()
            .filter(|check| check.outcome == Outcome::Passed)
            .count()
    }

    /// How many obligations the provider declined.
    #[must_use]
    pub fn skipped(&self) -> usize {
        self.checks
            .iter()
            .filter(|check| matches!(check.outcome, Outcome::Skipped { .. }))
            .count()
    }

    fn record(&mut self, name: &'static str, outcome: Outcome) {
        self.checks.push(Check { name, outcome });
    }
}

impl fmt::Display for ConformanceReport {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        writeln!(
            f,
            "conformance: {} passed, {} skipped, {} failed",
            self.exercised(),
            self.skipped(),
            self.checks.len() - self.exercised() - self.skipped()
        )?;
        for check in &self.checks {
            match &check.outcome {
                Outcome::Passed => writeln!(f, "  pass  {}", check.name)?,
                Outcome::Skipped { reason } => {
                    writeln!(f, "  skip  {} ({reason})", check.name)?;
                }
                Outcome::Failed { detail } => {
                    writeln!(f, "  FAIL  {} -- {detail}", check.name)?;
                }
            }
        }
        Ok(())
    }
}

/// Outward-oriented axis-aligned box.
#[must_use]
pub fn box_at(min: [f64; 3], max: [f64; 3]) -> TriMesh {
    let [x0, y0, z0] = min;
    let [x1, y1, z1] = max;
    let positions = vec![
        [x0, y0, z0].into(),
        [x1, y0, z0].into(),
        [x1, y1, z0].into(),
        [x0, y1, z0].into(),
        [x0, y0, z1].into(),
        [x1, y0, z1].into(),
        [x1, y1, z1].into(),
        [x0, y1, z1].into(),
    ];
    let indices = vec![
        0, 2, 1, 0, 3, 2, 4, 5, 6, 4, 6, 7, 0, 1, 5, 0, 5, 4, 1, 2, 6, 1, 6, 5, 2, 3, 7, 2, 7, 6,
        3, 0, 4, 3, 4, 7,
    ];
    TriMesh::new(positions, indices)
}

/// Enclosed volume by the divergence theorem.
#[must_use]
pub fn volume(mesh: &TriMesh) -> f64 {
    mesh.indices
        .as_chunks::<3>()
        .0
        .iter()
        .map(|t| {
            let a = mesh.positions[t[0] as usize];
            let b = mesh.positions[t[1] as usize];
            let c = mesh.positions[t[2] as usize];
            a.dot(b.cross(c))
        })
        .sum::<f64>()
        / 6.0
}

fn options() -> ExecutionOptions {
    ExecutionOptions::new(Tolerance::METRE)
}

/// Run an operation, mapping a refusal to a skip rather than a failure.
fn attempt(
    provider: &impl MeshBoolean,
    subject: &TriMesh,
    tool: &TriMesh,
    operation: BooleanOperator,
) -> Result<BooleanOutcome, Outcome> {
    match provider.boolean(subject, tool, operation, &options()) {
        Ok(outcome) => Ok(outcome),
        Err(GeomError::Unsupported { .. }) => Err(Outcome::Skipped {
            reason: format!("provider does not support {operation:?}"),
        }),
        Err(error) => Err(Outcome::Failed {
            detail: format!("{operation:?} on admissible operands failed: {error}"),
        }),
    }
}

/// Run every conformance obligation against `provider`.
///
/// Registration should be gated on the resulting report being conformant.
#[must_use]
pub fn run(provider: &impl MeshBoolean) -> ConformanceReport {
    let mut report = ConformanceReport::default();

    check_disjoint_algebra(provider, &mut report);
    check_identical_operands(provider, &mut report);
    check_difference_is_ordered(provider, &mut report);
    check_empty_result_is_a_value(provider, &mut report);
    check_inadmissible_operands_rejected(provider, &mut report);
    check_evidence_is_populated(provider, &mut report);
    check_determinism(provider, &mut report);
    check_determinism_is_not_overstated(provider, &mut report);
    check_declared_granularity_is_honoured(provider, &mut report);

    report
}

/// Disjoint operands: union adds, intersection empties, difference preserves.
fn check_disjoint_algebra(provider: &impl MeshBoolean, report: &mut ConformanceReport) {
    let a = box_at([0.0, 0.0, 0.0], [1.0, 1.0, 1.0]);
    let b = box_at([9.0, 9.0, 9.0], [10.0, 10.0, 10.0]);

    let outcome = match attempt(provider, &a, &b, BooleanOperator::Union) {
        Ok(outcome) => outcome,
        Err(skip_or_fail) => return report.record("disjoint_union_sums_volume", skip_or_fail),
    };
    let measured = volume(&outcome.mesh);
    report.record(
        "disjoint_union_sums_volume",
        if (measured - 2.0).abs() < 1e-9 {
            Outcome::Passed
        } else {
            Outcome::Failed {
                detail: format!("expected volume 2.0 for two disjoint unit cubes, got {measured}"),
            }
        },
    );

    match attempt(provider, &a, &b, BooleanOperator::Difference) {
        Ok(outcome) => {
            let measured = volume(&outcome.mesh);
            report.record(
                "disjoint_difference_preserves_subject",
                if (measured - 1.0).abs() < 1e-9 {
                    Outcome::Passed
                } else {
                    Outcome::Failed {
                        detail: format!("A \\ B must equal A when disjoint, got volume {measured}"),
                    }
                },
            );
        }
        Err(skip_or_fail) => {
            report.record("disjoint_difference_preserves_subject", skip_or_fail);
        }
    }
}

/// `A ∪ A = A` and `A \ A = ∅`.
fn check_identical_operands(provider: &impl MeshBoolean, report: &mut ConformanceReport) {
    let a = box_at([0.0, 0.0, 0.0], [1.0, 1.0, 1.0]);

    match attempt(provider, &a, &a, BooleanOperator::Union) {
        Ok(outcome) => {
            let measured = volume(&outcome.mesh);
            report.record(
                "self_union_is_idempotent",
                if (measured - 1.0).abs() < 1e-9 {
                    Outcome::Passed
                } else {
                    Outcome::Failed {
                        detail: format!("A union A must equal A, got volume {measured}"),
                    }
                },
            );
        }
        Err(skip_or_fail) => report.record("self_union_is_idempotent", skip_or_fail),
    }

    match attempt(provider, &a, &a, BooleanOperator::Difference) {
        Ok(outcome) => {
            let measured = volume(&outcome.mesh).abs();
            report.record(
                "self_difference_annihilates",
                if measured < 1e-9 {
                    Outcome::Passed
                } else {
                    Outcome::Failed {
                        detail: format!("A minus A must be empty, got volume {measured}"),
                    }
                },
            );
        }
        Err(skip_or_fail) => report.record("self_difference_annihilates", skip_or_fail),
    }
}

/// Difference is the only non-commutative operand, and must behave that way.
fn check_difference_is_ordered(provider: &impl MeshBoolean, report: &mut ConformanceReport) {
    let outer = box_at([0.0, 0.0, 0.0], [4.0, 4.0, 4.0]);
    let inner = box_at([1.0, 1.0, 1.0], [2.0, 2.0, 2.0]);

    let forward = attempt(provider, &outer, &inner, BooleanOperator::Difference);
    let reverse = attempt(provider, &inner, &outer, BooleanOperator::Difference);
    match (forward, reverse) {
        (Ok(forward), Ok(reverse)) => {
            let (big, small) = (volume(&forward.mesh), volume(&reverse.mesh).abs());
            report.record(
                "difference_respects_operand_order",
                if (big - 63.0).abs() < 1e-9 && small < 1e-9 {
                    Outcome::Passed
                } else {
                    Outcome::Failed {
                        detail: format!(
                            "outer minus inner should be 63.0 and inner minus outer empty, \
                             got {big} and {small}"
                        ),
                    }
                },
            );
        }
        (Err(skip_or_fail), _) | (_, Err(skip_or_fail)) => {
            report.record("difference_respects_operand_order", skip_or_fail);
        }
    }
}

/// An empty result is a value, never an error.
fn check_empty_result_is_a_value(provider: &impl MeshBoolean, report: &mut ConformanceReport) {
    let a = box_at([0.0, 0.0, 0.0], [1.0, 1.0, 1.0]);
    let b = box_at([9.0, 9.0, 9.0], [10.0, 10.0, 10.0]);

    match attempt(provider, &a, &b, BooleanOperator::Intersection) {
        Ok(outcome) => report.record(
            "empty_intersection_is_a_value_not_an_error",
            if outcome.mesh.triangle_count() == 0 {
                Outcome::Passed
            } else {
                Outcome::Failed {
                    detail: format!(
                        "disjoint intersection must be empty, got {} triangles",
                        outcome.mesh.triangle_count()
                    ),
                }
            },
        ),
        Err(skip_or_fail) => {
            report.record("empty_intersection_is_a_value_not_an_error", skip_or_fail);
        }
    }
}

/// Inadmissible operands must be refused, not silently processed.
fn check_inadmissible_operands_rejected(
    provider: &impl MeshBoolean,
    report: &mut ConformanceReport,
) {
    let good = box_at([0.0, 0.0, 0.0], [1.0, 1.0, 1.0]);
    let empty = TriMesh::new(Vec::new(), Vec::new());

    let result = provider.boolean(&good, &empty, BooleanOperator::Union, &options());
    report.record(
        "inadmissible_operand_is_refused",
        match result {
            Err(GeomError::InvalidInput(_) | GeomError::Degenerate(_)) => Outcome::Passed,
            Err(GeomError::Unsupported { .. }) => Outcome::Skipped {
                reason: "provider declined before validating".into(),
            },
            Err(other) => Outcome::Failed {
                detail: format!("expected InvalidInput or Degenerate, got {other:?}"),
            },
            Ok(_) => Outcome::Failed {
                detail: "an empty mesh is not a solid and must be refused".into(),
            },
        },
    );
}

/// Evidence must describe the actual work, not be left at defaults.
fn check_evidence_is_populated(provider: &impl MeshBoolean, report: &mut ConformanceReport) {
    let a = box_at([0.0, 0.0, 0.0], [1.0, 1.0, 1.0]);
    let b = box_at([9.0, 9.0, 9.0], [10.0, 10.0, 10.0]);

    match attempt(provider, &a, &b, BooleanOperator::Union) {
        Ok(outcome) => {
            let evidence = &outcome.evidence;
            let mut problems = Vec::new();
            if evidence.subject_triangles != a.triangle_count() {
                problems.push(format!(
                    "subject_triangles {} != {}",
                    evidence.subject_triangles,
                    a.triangle_count()
                ));
            }
            if evidence.output_triangles != outcome.mesh.triangle_count() {
                problems.push(format!(
                    "output_triangles {} != actual {}",
                    evidence.output_triangles,
                    outcome.mesh.triangle_count()
                ));
            }
            if evidence.sub_operations == 0 {
                problems.push("sub_operations must be at least 1".into());
            }
            report.record(
                "evidence_describes_the_work",
                if problems.is_empty() {
                    Outcome::Passed
                } else {
                    Outcome::Failed {
                        detail: problems.join("; "),
                    }
                },
            );
        }
        Err(skip_or_fail) => report.record("evidence_describes_the_work", skip_or_fail),
    }
}

/// The same inputs must produce the same output, every time.
fn check_determinism(provider: &impl MeshBoolean, report: &mut ConformanceReport) {
    let a = box_at([0.0, 0.0, 0.0], [4.0, 4.0, 4.0]);
    let b = box_at([1.0, 1.0, 1.0], [2.0, 2.0, 2.0]);

    let first = attempt(provider, &a, &b, BooleanOperator::Difference);
    let second = attempt(provider, &a, &b, BooleanOperator::Difference);
    match (first, second) {
        (Ok(first), Ok(second)) => report.record(
            "repeated_calls_are_deterministic",
            if first.mesh.positions == second.mesh.positions
                && first.mesh.indices == second.mesh.indices
            {
                Outcome::Passed
            } else {
                Outcome::Failed {
                    detail: "identical inputs produced different geometry".into(),
                }
            },
        ),
        (Err(skip_or_fail), _) | (_, Err(skip_or_fail)) => {
            report.record("repeated_calls_are_deterministic", skip_or_fail);
        }
    }
}

/// A provider must not declare a reproducibility level it cannot honour.
///
/// `check_determinism` above compares two calls in ONE process, so it cannot
/// see the failure mode that matters here: `std`'s `RandomState` seeds each
/// `HashMap` per process, so a provider deduping through one is stable within
/// a run and unstable across runs. Two identical calls agree; a rebuild
/// disagrees.
///
/// So the declaration is checked structurally instead. `Determinism::Bitwise`
/// is the level that promises cross-process byte equality, and it is exactly
/// the level a same-process test cannot verify. A provider claiming it must
/// have earned it by construction — ordered vertex identity, no hash
/// iteration, no thread-order-dependent accumulation — and the conformance
/// suite has no way to confirm that from the outside.
///
/// Refusing the claim here is the fail-closed choice: an overstated level
/// makes `Plan::admit` accept a step it should have refused, which is worse
/// than an understated one that merely refuses work the provider could have
/// done.
fn check_determinism_is_not_overstated(
    provider: &impl MeshBoolean,
    report: &mut ConformanceReport,
) {
    let declared = provider.determinism();
    report.record(
        "determinism_declaration_is_verifiable",
        if declared == Determinism::Bitwise {
            Outcome::Failed {
                detail: "provider declares Determinism::Bitwise, which this suite cannot \
                         verify: cross-process byte equality is not observable from a \
                         single-process check. Declare the strongest level the provider \
                         can prove by construction instead"
                    .into(),
            }
        } else {
            Outcome::Passed
        },
    );
}

/// A pre-cancelled token must stop the operation, if the provider polls at all.
fn check_declared_granularity_is_honoured(
    provider: &impl MeshBoolean,
    report: &mut ConformanceReport,
) {
    use axiolid_contracts::{CancellationGranularity, CancellationToken};

    let granularity = provider.cancellation_granularity();
    if granularity == CancellationGranularity::None {
        report.record(
            "declared_cancellation_is_honoured",
            Outcome::Skipped {
                reason: "provider declares it does not poll".into(),
            },
        );
        return;
    }

    let a = box_at([0.0, 0.0, 0.0], [1.0, 1.0, 1.0]);
    let b = box_at([9.0, 9.0, 9.0], [10.0, 10.0, 10.0]);
    let token = CancellationToken::new();
    token.cancel();
    let cancelled = options().with_cancellation(token);

    let result = provider.boolean(&a, &b, BooleanOperator::Union, &cancelled);
    report.record(
        "declared_cancellation_is_honoured",
        match result {
            Err(GeomError::Cancelled) => Outcome::Passed,
            Err(GeomError::Unsupported { .. }) => Outcome::Skipped {
                reason: "provider declined the operation".into(),
            },
            Ok(_) => Outcome::Failed {
                detail: format!(
                    "provider declares {granularity:?} polling but ran to completion \
                     with a pre-cancelled token"
                ),
            },
            Err(other) => Outcome::Failed {
                detail: format!("expected Cancelled, got {other:?}"),
            },
        },
    );
}
