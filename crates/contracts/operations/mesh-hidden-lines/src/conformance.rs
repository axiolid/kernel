//! Shared conformance checks for portable hidden-line providers.
//!
//! The suite calls the provider directly, so a provider must refuse bad
//! input itself and not rely on a registry to do it. Each drawing case
//! states its expected counts and lengths per edge kind and visibility, and
//! the exact pieces of chosen edges; the scenes and the arithmetic behind
//! them are in the (private) fixtures, checked by their own unit tests.

use core::fmt;

use axiolid_contracts::ExecutionOptions;
use axiolid_core::{Point2, Tolerance};

use crate::fixtures::{cases, refusals, reversed, Case, Fixture, Piece, LENGTH_TOLERANCE};
use crate::{
    EdgeKind, HiddenLineError, HiddenLineMesh, HiddenLineOutcome, HiddenLineRequest,
    MeshHiddenLines, MeshId, SourceEdge, Visibility,
};

/// One piece of an edge, as expected or as returned.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct EdgePiece {
    pub kind: EdgeKind,
    pub visibility: Visibility,
    pub start: Point2,
    pub end: Point2,
}

#[derive(Debug, Clone, PartialEq)]
pub enum ConformanceFailure {
    /// The provider refused a request it must draw, or refused a request
    /// that it must refuse with a different error.
    ReturnedError {
        case: &'static str,
        error: HiddenLineError,
    },
    /// The provider drew a request it must refuse.
    MissingRefusal { case: &'static str },
    /// The refusal names the wrong mesh, or none.
    WrongRefusalMesh {
        case: &'static str,
        expected: Option<MeshId>,
        actual: Option<MeshId>,
    },
    /// The total length of the segments of one kind and visibility is wrong.
    WrongLength {
        case: &'static str,
        kind: EdgeKind,
        visibility: Visibility,
        expected: f64,
        actual: f64,
    },
    /// The number of segments of one kind and visibility is wrong.
    WrongSegmentCount {
        case: &'static str,
        kind: EdgeKind,
        visibility: Visibility,
        expected: usize,
        actual: usize,
    },
    /// The pieces of one source edge differ from the expected ones.
    WrongEdge {
        case: &'static str,
        edge: SourceEdge,
        expected: Vec<EdgePiece>,
        actual: Vec<EdgePiece>,
    },
    /// Evidence disagrees with the request, the outcome or the scene.
    IncorrectEvidence { case: &'static str },
    /// Two identical requests returned different outcomes.
    NonDeterministic { case: &'static str },
    /// Reversing every triangle changed the outcome.
    OrientationDependent { case: &'static str },
}

#[derive(Debug, Clone, Default, PartialEq)]
pub struct ConformanceReport {
    pub failures: Vec<ConformanceFailure>,
}

impl ConformanceReport {
    pub fn is_success(&self) -> bool {
        self.failures.is_empty()
    }
}

impl fmt::Display for ConformanceReport {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        if self.is_success() {
            write!(f, "conformant")
        } else {
            write!(f, "{:?}", self.failures)
        }
    }
}

pub struct ConformanceSuite;

impl ConformanceSuite {
    pub fn run(provider: &dyn MeshHiddenLines) -> ConformanceReport {
        let mut report = ConformanceReport::default();
        let options = ExecutionOptions::new(Tolerance::new(1e-9, 1e-9).expect("valid tolerance"));

        for case in cases() {
            let first = draw(provider, &case, &case.meshes, &options);
            let second = draw(provider, &case, &case.meshes, &options);
            match (first, second) {
                (Ok(a), Ok(b)) => {
                    check_outcome(&mut report, &case, &a);
                    if a != b {
                        report
                            .failures
                            .push(ConformanceFailure::NonDeterministic { case: case.name });
                    }
                    if case.orientation_independent {
                        check_orientation(provider, &mut report, &case, &a, &options);
                    }
                }
                (Err(error), _) | (_, Err(error)) => {
                    report.failures.push(ConformanceFailure::ReturnedError {
                        case: case.name,
                        error,
                    })
                }
            }
        }

        for refusal in refusals() {
            let meshes: Vec<HiddenLineMesh<'_>> = refusal
                .meshes
                .iter()
                .map(|f| HiddenLineMesh::new(f.id, f.role, &f.mesh))
                .collect();
            let request =
                HiddenLineRequest::new(&meshes, refusal.view, refusal.edges, refusal.limits);
            match provider.hidden_lines(&request, &options) {
                Ok(_) => report
                    .failures
                    .push(ConformanceFailure::MissingRefusal { case: refusal.name }),
                Err(error) if !(refusal.accepts)(&error.error) => {
                    report.failures.push(ConformanceFailure::ReturnedError {
                        case: refusal.name,
                        error,
                    })
                }
                Err(error) if error.mesh != refusal.mesh => {
                    report.failures.push(ConformanceFailure::WrongRefusalMesh {
                        case: refusal.name,
                        expected: refusal.mesh,
                        actual: error.mesh,
                    })
                }
                Err(_) => {}
            }
        }

        report
    }
}

fn draw(
    provider: &dyn MeshHiddenLines,
    case: &Case,
    fixtures: &[Fixture],
    options: &ExecutionOptions,
) -> Result<HiddenLineOutcome, HiddenLineError> {
    let meshes: Vec<HiddenLineMesh<'_>> = fixtures
        .iter()
        .map(|f| HiddenLineMesh::new(f.id, f.role, &f.mesh))
        .collect();
    let mut request = HiddenLineRequest::new(&meshes, case.view, case.edges, case.limits);
    if let Some(region) = case.region {
        request = request.with_region(region);
    }
    provider.hidden_lines(&request, options)
}

fn check_orientation(
    provider: &dyn MeshHiddenLines,
    report: &mut ConformanceReport,
    case: &Case,
    expected: &HiddenLineOutcome,
    options: &ExecutionOptions,
) {
    let flipped: Vec<Fixture> = case
        .meshes
        .iter()
        .map(|f| Fixture {
            id: f.id,
            role: f.role,
            mesh: reversed(&f.mesh),
        })
        .collect();
    match draw(provider, case, &flipped, options) {
        Ok(outcome) if same_drawing(&outcome, expected) => {}
        Ok(_) => report
            .failures
            .push(ConformanceFailure::OrientationDependent { case: case.name }),
        Err(error) => report.failures.push(ConformanceFailure::ReturnedError {
            case: case.name,
            error,
        }),
    }
}

/// Equal segments and evidence, except exact_escalations, which depends on
/// the order in which the vertices are met.
fn same_drawing(a: &HiddenLineOutcome, b: &HiddenLineOutcome) -> bool {
    let mut evidence = a.evidence;
    evidence.exact_escalations = b.evidence.exact_escalations;
    a.view == b.view && a.segments == b.segments && evidence == b.evidence
}

fn length(start: Point2, end: Point2) -> f64 {
    (end - start).length()
}

fn check_outcome(report: &mut ConformanceReport, case: &Case, outcome: &HiddenLineOutcome) {
    let name = case.name;
    let count = |visibility: Visibility| {
        outcome
            .segments
            .iter()
            .filter(|s| s.visibility == visibility)
            .count()
    };
    let evidence = &outcome.evidence;
    let expect = &case.expect;
    let triangles: usize = case.meshes.iter().map(|f| f.mesh.triangle_count()).sum();
    let evidence_ok = outcome.view == case.view
        && evidence.is_derived_from_input_mesh()
        && evidence.source_meshes == case.meshes.len()
        && evidence.source_triangles == triangles
        && evidence.visible_segments == count(Visibility::Visible)
        && evidence.hidden_segments == count(Visibility::Hidden)
        && (case.edges.hidden || evidence.hidden_segments == 0)
        && expect
            .selected_edges
            .is_none_or(|n| evidence.selected_edges == n)
        && expect.point_edges.is_none_or(|n| evidence.point_edges == n)
        && expect
            .clipped_edges
            .is_none_or(|n| evidence.clipped_edges == n);
    if !evidence_ok {
        report
            .failures
            .push(ConformanceFailure::IncorrectEvidence { case: name });
    }

    if let Some(rows) = &expect.rows {
        for kind in [EdgeKind::Boundary, EdgeKind::Silhouette, EdgeKind::Crease] {
            for visibility in [Visibility::Visible, Visibility::Hidden] {
                let (want_count, want_length) = rows
                    .iter()
                    .find(|r| r.kind == kind && r.visibility == visibility)
                    .map_or((0, 0.0), |r| (r.count, r.length));
                let picked = outcome
                    .segments
                    .iter()
                    .filter(|s| s.kind == kind && s.visibility == visibility);
                let (got_count, got_length) = picked.fold((0, 0.0), |(n, sum), s| {
                    (n + 1, sum + length(s.start, s.end))
                });
                if got_count != want_count {
                    report.failures.push(ConformanceFailure::WrongSegmentCount {
                        case: name,
                        kind,
                        visibility,
                        expected: want_count,
                        actual: got_count,
                    });
                }
                if (got_length - want_length).abs() > LENGTH_TOLERANCE {
                    report.failures.push(ConformanceFailure::WrongLength {
                        case: name,
                        kind,
                        visibility,
                        expected: want_length,
                        actual: got_length,
                    });
                }
            }
        }
    }

    for wanted in &expect.edges {
        let expected: Vec<EdgePiece> = wanted
            .pieces
            .iter()
            .map(|p| to_piece(wanted.kind, p))
            .collect();
        let actual: Vec<EdgePiece> = outcome
            .segments
            .iter()
            .filter(|s| s.source == wanted.edge)
            .map(|s| EdgePiece {
                kind: s.kind,
                visibility: s.visibility,
                start: s.start,
                end: s.end,
            })
            .collect();
        let same = expected.len() == actual.len()
            && expected.iter().zip(&actual).all(|(e, a)| {
                e.kind == a.kind
                    && e.visibility == a.visibility
                    && length(e.start, a.start) <= LENGTH_TOLERANCE
                    && length(e.end, a.end) <= LENGTH_TOLERANCE
            });
        if !same {
            report.failures.push(ConformanceFailure::WrongEdge {
                case: name,
                edge: wanted.edge,
                expected,
                actual,
            });
        }
    }
}

fn to_piece(kind: EdgeKind, piece: &Piece) -> EdgePiece {
    EdgePiece {
        kind,
        visibility: piece.visibility,
        start: Point2::new(piece.start[0], piece.start[1]),
        end: Point2::new(piece.end[0], piece.end[1]),
    }
}

#[cfg(test)]
mod tests {
    use axiolid_contracts::{Backend, BackendDescriptor, BackendId, ExecutionTarget, GeomError};

    use super::*;

    /// Refuses everything: the offending mesh (the second of two) for the
    /// per-mesh refusals, or the whole request when `blame_offender` is off.
    #[derive(Debug)]
    struct Refuser {
        blame_offender: bool,
    }

    impl Backend for Refuser {
        fn descriptor(&self) -> BackendDescriptor {
            BackendDescriptor::new(
                BackendId::new("conformance-refuser"),
                ExecutionTarget::PortableCpu,
            )
        }
    }

    impl MeshHiddenLines for Refuser {
        fn hidden_lines(
            &self,
            request: &HiddenLineRequest<'_>,
            _options: &ExecutionOptions,
        ) -> Result<HiddenLineOutcome, HiddenLineError> {
            if request.meshes.len() == 2 {
                let error = GeomError::NotManifold("refused".into());
                return Err(if self.blame_offender {
                    HiddenLineError::mesh(request.meshes[1].id, error)
                } else {
                    HiddenLineError::request(error)
                });
            }
            if request.limits.max_output_segments == 1 {
                return Err(GeomError::BudgetExceeded {
                    resource: "hidden-line output segments",
                }
                .into());
            }
            if request.limits.max_candidate_pairs == 0 {
                return Err(GeomError::BudgetExceeded {
                    resource: "hidden-line candidate pairs",
                }
                .into());
            }
            Err(GeomError::Cancelled.into())
        }
    }

    const REFUSAL_CASES: [&str; 5] = [
        "three triangles on one edge",
        "inconsistent winding",
        "sheet declared solid",
        "output limit",
        "candidate limit",
    ];

    fn refusal_failures(report: &ConformanceReport) -> Vec<&ConformanceFailure> {
        report
            .failures
            .iter()
            .filter(|failure| match failure {
                ConformanceFailure::ReturnedError { case, .. }
                | ConformanceFailure::MissingRefusal { case }
                | ConformanceFailure::WrongRefusalMesh { case, .. } => REFUSAL_CASES.contains(case),
                _ => false,
            })
            .collect()
    }

    #[test]
    fn a_provider_that_names_the_offending_mesh_passes_the_refusal_cases() {
        let report = ConformanceSuite::run(&Refuser {
            blame_offender: true,
        });
        assert!(refusal_failures(&report).is_empty(), "{report}");
        // It draws nothing, so every drawing case fails.
        assert!(!report.is_success());
    }

    #[test]
    fn a_refusal_that_loses_the_mesh_id_fails_exactly_the_mesh_cases() {
        let report = ConformanceSuite::run(&Refuser {
            blame_offender: false,
        });
        let wrong = refusal_failures(&report);
        assert_eq!(wrong.len(), 3, "{wrong:?}");
        for failure in wrong {
            assert!(matches!(
                failure,
                ConformanceFailure::WrongRefusalMesh {
                    expected: Some(MeshId(2)),
                    actual: None,
                    ..
                }
            ));
        }
    }

    /// Draws nothing and refuses nothing.
    #[derive(Debug)]
    struct Silent;

    impl Backend for Silent {
        fn descriptor(&self) -> BackendDescriptor {
            BackendDescriptor::new(
                BackendId::new("conformance-silent"),
                ExecutionTarget::PortableCpu,
            )
        }
    }

    impl MeshHiddenLines for Silent {
        fn hidden_lines(
            &self,
            request: &HiddenLineRequest<'_>,
            _options: &ExecutionOptions,
        ) -> Result<HiddenLineOutcome, HiddenLineError> {
            let triangles = request.meshes.iter().map(|m| m.mesh.triangle_count()).sum();
            Ok(HiddenLineOutcome::new(
                request.view,
                Vec::new(),
                crate::HiddenLineEvidence::input_mesh(request.meshes.len(), triangles),
            ))
        }
    }

    #[test]
    fn an_empty_drawing_and_a_missing_refusal_are_reported() {
        let report = ConformanceSuite::run(&Silent);
        let missing = report
            .failures
            .iter()
            .filter(|f| matches!(f, ConformanceFailure::MissingRefusal { .. }))
            .count();
        assert_eq!(missing, REFUSAL_CASES.len());
        assert!(report
            .failures
            .iter()
            .any(|f| matches!(f, ConformanceFailure::WrongSegmentCount { .. })));
        assert!(report
            .failures
            .iter()
            .any(|f| matches!(f, ConformanceFailure::WrongEdge { .. })));
    }
}
