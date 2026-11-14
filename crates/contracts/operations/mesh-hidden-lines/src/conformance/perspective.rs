//! Conformance checks for providers that draw perspective views (ADR 0107).
//!
//! The base suites cannot grow a field ([`ConformanceReport`]
//! is exhaustive), so perspective has a suite of its own, as the region
//! suite does for section requests (ADR 0091). The claim is inferred from
//! `supports_projection`, never from probing: a provider that answers `false`
//! for [`ViewProjection::Perspective`] does not claim it, and must then refuse a
//! valid perspective request with `Unsupported` or `UnsupportedInput`, so
//! that a registry can fall back past it. A provider that answers `true`
//! claims it, and the suite then draws the perspective fixtures (a box, a box
//! cut by the near plane, an eye inside a box, an exact boundary, an edge
//! along the eye ray, an eye inside an edge-on triangle) and checks the drawing and the refusal of a bad near
//! depth. Hidden-line fixtures exist; the visible-face claim has none yet and
//! fails closed with [`PerspectiveConformanceFailure::ClaimUnverifiable`]
//! until they land.

use core::fmt;

use axiolid_contracts::{ExecutionOptions, GeomError};
use axiolid_core::Tolerance;

use super::faces::FaceConformanceFailure;
use super::{check_outcome, draw, mesh_of, same_drawing, ConformanceFailure, ConformanceReport};
use crate::fixtures::{
    box_mesh, limits, perspective_cases, perspective_refusals, plan_frame, reversed, Fixture,
};
use crate::{
    EdgeSelection, HiddenLineError, HiddenLineMesh, HiddenLineRequest, MeshHiddenLines, MeshId,
    MeshRole, MeshVisibleFaces, View, ViewProjection, VisibleFaceRequest,
};

/// One failure of the perspective suite.
#[derive(Debug, Clone, PartialEq)]
#[non_exhaustive]
pub enum PerspectiveConformanceFailure {
    /// The provider fails the hidden-line suite. Reported by a registry that
    /// runs both suites, never by [`PerspectiveConformanceSuite`] itself.
    HiddenLinesBase(ConformanceFailure),
    /// The provider fails the visible-face suite; see
    /// [`HiddenLinesBase`](Self::HiddenLinesBase).
    VisibleFacesBase(FaceConformanceFailure),
    /// The provider does not claim perspective but drew a valid perspective
    /// request, which it would have to ignore the projection to do.
    UnclaimedButDrawn {
        /// The operation, `"hidden lines"` or `"visible faces"`.
        operation: &'static str,
    },
    /// The provider does not claim perspective and refused a valid
    /// perspective request with something other than `Unsupported` or
    /// `UnsupportedInput`, which a registry would not fall back past.
    UnclaimedWrongRefusal {
        /// The operation, `"hidden lines"` or `"visible faces"`.
        operation: &'static str,
        /// The refusal.
        error: HiddenLineError,
    },
    /// A hidden-line provider that claims perspective failed a perspective
    /// fixture; the failure names the fixture.
    HiddenLines(ConformanceFailure),
    /// The provider claims perspective, which this suite cannot verify.
    ClaimUnverifiable {
        /// The operation, `"hidden lines"` or `"visible faces"`.
        operation: &'static str,
    },
}

/// Result of [`PerspectiveConformanceSuite`].
#[derive(Debug, Clone, PartialEq)]
#[non_exhaustive]
pub struct PerspectiveConformanceReport {
    /// Every failure found.
    pub failures: Vec<PerspectiveConformanceFailure>,
    /// Whether the provider claims perspective views, as
    /// `supports_projection(ViewProjection::Perspective)` says.
    pub claimed: bool,
}

impl PerspectiveConformanceReport {
    /// Whether no check failed.
    pub fn is_success(&self) -> bool {
        self.failures.is_empty()
    }
}

impl fmt::Display for PerspectiveConformanceReport {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        if self.is_success() {
            write!(
                f,
                "{}",
                if self.claimed {
                    "conformant"
                } else {
                    "not claimed"
                }
            )
        } else {
            write!(f, "{:?}", self.failures)
        }
    }
}

/// The perspective conformance suite for `MeshHiddenLines` and
/// `MeshVisibleFaces`.
pub struct PerspectiveConformanceSuite;

impl PerspectiveConformanceSuite {
    /// Check a hidden-line provider's perspective claim.
    pub fn run(provider: &dyn MeshHiddenLines) -> PerspectiveConformanceReport {
        let claimed = provider.supports_projection(ViewProjection::Perspective);
        let mut report = PerspectiveConformanceReport {
            failures: Vec::new(),
            claimed,
        };
        const OPERATION: &str = "hidden lines";
        if claimed {
            report.failures = claimed_failures(provider)
                .into_iter()
                .map(PerspectiveConformanceFailure::HiddenLines)
                .collect();
            return report;
        }
        let mesh = box_mesh([-1.0, -1.0, -4.0], [1.0, 1.0, -2.0]);
        let meshes = [HiddenLineMesh::new(MeshId(1), MeshRole::Solid, &mesh)];
        let request = HiddenLineRequest::new(
            &meshes,
            perspective_view(),
            EdgeSelection::new(Some(0.5)),
            limits(),
        );
        let result = provider.hidden_lines(&request, &options());
        record_unclaimed(&mut report, OPERATION, result.map(|_| ()));
        report
    }

    /// Check a visible-face provider's perspective claim.
    pub fn run_faces(provider: &dyn MeshVisibleFaces) -> PerspectiveConformanceReport {
        let claimed = provider.supports_projection(ViewProjection::Perspective);
        let mut report = PerspectiveConformanceReport {
            failures: Vec::new(),
            claimed,
        };
        const OPERATION: &str = "visible faces";
        if claimed {
            report
                .failures
                .push(PerspectiveConformanceFailure::ClaimUnverifiable {
                    operation: OPERATION,
                });
            return report;
        }
        let mesh = box_mesh([-1.0, -1.0, -4.0], [1.0, 1.0, -2.0]);
        let meshes = [HiddenLineMesh::new(MeshId(1), MeshRole::Solid, &mesh)];
        let request =
            VisibleFaceRequest::new(&meshes, perspective_view(), crate::face_fixtures::limits());
        let result = provider.visible_faces(&request, &options());
        record_unclaimed(&mut report, OPERATION, result.map(|_| ()));
        report
    }
}

/// The failures of a provider that claims perspective against the
/// perspective fixtures. The provider is called directly, as in the base
/// suite, so it must refuse a bad near depth itself.
fn claimed_failures(provider: &dyn MeshHiddenLines) -> Vec<ConformanceFailure> {
    let options = options();
    let mut report = ConformanceReport::default();
    for case in perspective_cases() {
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
                    let flipped: Vec<Fixture> = case
                        .meshes
                        .iter()
                        .map(|f| Fixture {
                            id: f.id,
                            role: f.role,
                            mesh: reversed(&f.mesh),
                            groups: f.groups.clone(),
                            planar: f.planar,
                        })
                        .collect();
                    match draw(provider, &case, &flipped, &options) {
                        Ok(outcome) if same_drawing(&outcome, &a) => {}
                        Ok(_) => report
                            .failures
                            .push(ConformanceFailure::OrientationDependent { case: case.name }),
                        Err(error) => report.failures.push(ConformanceFailure::ReturnedError {
                            case: case.name,
                            error,
                        }),
                    }
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
    for refusal in perspective_refusals() {
        let meshes: Vec<HiddenLineMesh<'_>> = refusal.meshes.iter().map(mesh_of).collect();
        let request = HiddenLineRequest::new(&meshes, refusal.view, refusal.edges, refusal.limits);
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
            Err(_) => {}
        }
    }
    report.failures
}

fn perspective_view() -> View {
    View::perspective(plan_frame(), -1.0).expect("a negative finite near")
}

fn options() -> ExecutionOptions {
    ExecutionOptions::new(Tolerance::new(1e-9, 1e-9).expect("valid tolerance"))
}

fn record_unclaimed(
    report: &mut PerspectiveConformanceReport,
    operation: &'static str,
    result: Result<(), HiddenLineError>,
) {
    match result {
        Ok(()) => report
            .failures
            .push(PerspectiveConformanceFailure::UnclaimedButDrawn { operation }),
        Err(HiddenLineError {
            error: GeomError::Unsupported { .. } | GeomError::UnsupportedInput { .. },
            ..
        }) => {}
        Err(error) => report
            .failures
            .push(PerspectiveConformanceFailure::UnclaimedWrongRefusal { operation, error }),
    }
}

#[cfg(test)]
mod tests {
    use axiolid_contracts::{Backend, BackendDescriptor, BackendId, ExecutionTarget, Operation};

    use super::*;
    use crate::{HiddenLineOutcome, VisibleFaceOutcome};

    /// A provider whose behaviour on a perspective request is chosen.
    #[derive(Debug)]
    struct Probe {
        claims: bool,
        answer: Answer,
    }

    #[derive(Debug, Clone, Copy)]
    enum Answer {
        Unsupported,
        UnsupportedInput,
        Invalid,
        Draw,
    }

    impl Backend for Probe {
        fn descriptor(&self) -> BackendDescriptor {
            BackendDescriptor::new(
                BackendId::new("perspective-probe"),
                ExecutionTarget::PortableCpu,
            )
        }
    }

    impl Probe {
        fn refusal(&self, operation: Operation) -> Option<HiddenLineError> {
            let backend = BackendId::new("perspective-probe");
            match self.answer {
                Answer::Unsupported => Some(GeomError::Unsupported { backend, operation }.into()),
                Answer::UnsupportedInput => Some(
                    GeomError::UnsupportedInput {
                        backend,
                        operation,
                        input: "perspective view",
                    }
                    .into(),
                ),
                Answer::Invalid => Some(GeomError::InvalidInput("no".into()).into()),
                Answer::Draw => None,
            }
        }
    }

    impl MeshHiddenLines for Probe {
        fn supports_projection(&self, projection: ViewProjection) -> bool {
            projection == ViewProjection::Orthographic || self.claims
        }

        fn hidden_lines(
            &self,
            request: &HiddenLineRequest<'_>,
            _options: &ExecutionOptions,
        ) -> Result<HiddenLineOutcome, HiddenLineError> {
            match self.refusal(Operation::MeshHiddenLines) {
                Some(error) => Err(error),
                None => Ok(HiddenLineOutcome::new(
                    request.view,
                    Vec::new(),
                    crate::HiddenLineEvidence::input_mesh(0, 0),
                )),
            }
        }
    }

    impl MeshVisibleFaces for Probe {
        fn supports_projection(&self, projection: ViewProjection) -> bool {
            projection == ViewProjection::Orthographic || self.claims
        }

        fn visible_faces(
            &self,
            request: &VisibleFaceRequest<'_>,
            _options: &ExecutionOptions,
        ) -> Result<VisibleFaceOutcome, HiddenLineError> {
            match self.refusal(Operation::MeshVisibleFaces) {
                Some(error) => Err(error),
                None => Ok(VisibleFaceOutcome::new(
                    request.view,
                    Vec::new(),
                    crate::VisibleFaceEvidence::input_mesh(0, 0),
                )),
            }
        }
    }

    fn both(probe: &Probe) -> [PerspectiveConformanceReport; 2] {
        [
            PerspectiveConformanceSuite::run(probe),
            PerspectiveConformanceSuite::run_faces(probe),
        ]
    }

    #[test]
    fn a_provider_that_refuses_by_name_passes_unclaimed() {
        for answer in [Answer::Unsupported, Answer::UnsupportedInput] {
            for report in both(&Probe {
                claims: false,
                answer,
            }) {
                assert!(!report.claimed);
                assert!(report.is_success(), "{report}");
                assert_eq!(report.to_string(), "not claimed");
            }
        }
    }

    #[test]
    fn an_unclaimed_provider_that_draws_or_refuses_otherwise_fails() {
        let [lines, faces] = both(&Probe {
            claims: false,
            answer: Answer::Draw,
        });
        assert_eq!(
            lines.failures,
            [PerspectiveConformanceFailure::UnclaimedButDrawn {
                operation: "hidden lines"
            }]
        );
        assert_eq!(
            faces.failures,
            [PerspectiveConformanceFailure::UnclaimedButDrawn {
                operation: "visible faces"
            }]
        );
        for report in both(&Probe {
            claims: false,
            answer: Answer::Invalid,
        }) {
            assert!(matches!(
                report.failures[..],
                [PerspectiveConformanceFailure::UnclaimedWrongRefusal { .. }]
            ));
        }
    }

    #[test]
    fn a_claim_is_inferred_from_supports_projection_and_verified_by_the_fixtures() {
        let probe = Probe {
            claims: true,
            answer: Answer::Draw,
        };
        let lines = PerspectiveConformanceSuite::run(&probe);
        assert!(lines.claimed);
        // The probe draws nothing, so every fixture with segments fails and
        // every bad near depth is drawn instead of refused.
        assert!(!lines.is_success());
        assert!(lines
            .failures
            .iter()
            .all(|f| matches!(f, PerspectiveConformanceFailure::HiddenLines(_))));
        for case in [
            "perspective box",
            "perspective box cut by the near plane",
            "perspective eye inside a box",
        ] {
            assert!(
                lines.failures.iter().any(|f| matches!(
                    f,
                    PerspectiveConformanceFailure::HiddenLines(
                        ConformanceFailure::WrongSegmentCount { case: c, .. }
                    ) if *c == case
                )),
                "{case}: {lines}"
            );
        }
        assert!(lines.failures.iter().any(|f| matches!(
            f,
            PerspectiveConformanceFailure::HiddenLines(ConformanceFailure::MissingRefusal {
                case: "perspective view without a near depth"
            })
        )));
    }

    #[test]
    fn a_visible_face_claim_still_fails_closed() {
        let faces = PerspectiveConformanceSuite::run_faces(&Probe {
            claims: true,
            answer: Answer::Draw,
        });
        assert!(faces.claimed);
        assert!(matches!(
            faces.failures[..],
            [PerspectiveConformanceFailure::ClaimUnverifiable {
                operation: "visible faces"
            }]
        ));
    }
}
