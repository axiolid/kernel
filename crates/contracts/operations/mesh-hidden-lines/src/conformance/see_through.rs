//! Conformance checks for providers that honour see-through triangles (ADR
//! 0111).
//!
//! The base reports are exhaustive, so see-through triangles have a suite of
//! their own, as perspective views do (ADR 0107). The claim is inferred from
//! `supports_occlusion`, never from probing: a provider that answers `false`
//! for [`Occlusion::SeeThrough`] does not claim it, and must then refuse a
//! see-through request with `Unsupported` or `UnsupportedInput`, so that a
//! registry can fall back past it. A provider that answers `true` claims it,
//! and the suite then draws the see-through fixtures (panes over a floor,
//! panes over panes, a pane under an opaque sheet, a see-through box, edges
//! under a pane's border and diagonal, the occluding cap of a cut
//! see-through solid and, when the provider also claims perspective, a
//! see-through slab seen from an eye), checks each piece's
//! [`seen_through`](crate::HiddenLineSegment::seen_through) count and the
//! evidence, requires a short occlusion slice to be refused naming its mesh,
//! and requires base scenes flagged all opaque to be drawn exactly as without
//! flags. A visible-face claim must reproduce the regions and `behind` lists
//! of eight scenes (a pane over a box, two panes over a box, a depth cycle of
//! three panes, a see-through box, a pane tied with a floor, a stack stopped
//! by an opaque sheet, the cap of a cut see-through solid and planar groups
//! split by occlusion), each also with every triangle reversed, refuse a mixed
//! face group, a short occlusion slice and a stack over
//! `max_behind_surfaces`, and answer base scenes flagged all opaque as
//! without flags.

use core::fmt;

use axiolid_contracts::{ExecutionOptions, GeomError};
use axiolid_core::Tolerance;

use super::faces::{
    answer, check_expected_evidence, check_orientation, check_shape, same_result,
    FaceConformanceFailure, FaceConformanceReport, Scene,
};
use super::{check_outcome, draw, mesh_of, same_drawing, ConformanceFailure, ConformanceReport};
use crate::face_fixtures::{
    cases as face_cases, see_through_face_cases, see_through_face_refusals, FaceFixture,
    NEUTRAL_FACE_CASES,
};
use crate::fixtures::{
    cases, see_through_cases, see_through_refusals, Case, Fixture, LENGTH_TOLERANCE, NEUTRAL_CASES,
};
use crate::{
    EdgeKind, HiddenLineError, HiddenLineMesh, HiddenLineOutcome, HiddenLineRequest,
    MeshHiddenLines, MeshVisibleFaces, Occlusion, SourceEdge, ViewProjection, Visibility,
    VisibleFaceRequest,
};

/// One failure of the see-through suite.
#[derive(Debug, Clone, PartialEq)]
#[non_exhaustive]
pub enum SeeThroughConformanceFailure {
    /// The provider does not claim see-through triangles but drew a
    /// see-through request, which it would have to let them hide to do.
    UnclaimedButDrawn {
        /// The operation, `"hidden lines"` or `"visible faces"`.
        operation: &'static str,
    },
    /// The provider does not claim see-through triangles and refused a
    /// see-through request with something other than `Unsupported` or
    /// `UnsupportedInput`, which a registry would not fall back past.
    UnclaimedWrongRefusal {
        /// The operation, `"hidden lines"` or `"visible faces"`.
        operation: &'static str,
        /// The refusal.
        error: HiddenLineError,
    },
    /// A hidden-line provider that claims see-through triangles failed a
    /// see-through fixture or refusal in its geometry, visibility, base
    /// evidence, determinism or orientation; the failure names the case.
    HiddenLines(ConformanceFailure),
    /// A visible-face provider that claims see-through triangles failed a
    /// see-through fixture or refusal: its regions (`behind` included),
    /// shape, base evidence, determinism or orientation; the failure names
    /// the case.
    VisibleFaces(FaceConformanceFailure),
    /// The segments of one kind, visibility and `seen_through` count differ
    /// in number or total length from the expected ones.
    WrongSeenThroughRow {
        /// The fixture.
        case: &'static str,
        /// The edge kind.
        kind: EdgeKind,
        /// The visibility.
        visibility: Visibility,
        /// The `seen_through` count.
        seen_through: u32,
        /// Expected number of segments.
        expected_count: usize,
        /// Returned number of segments.
        actual_count: usize,
        /// Expected total image length.
        expected_length: f64,
        /// Returned total image length.
        actual_length: f64,
    },
    /// The `seen_through` counts of the pieces of one source edge, in order,
    /// differ from the expected ones.
    WrongSeenThrough {
        /// The fixture.
        case: &'static str,
        /// The edge.
        edge: SourceEdge,
        /// Expected counts.
        expected: Vec<u32>,
        /// Returned counts.
        actual: Vec<u32>,
    },
    /// The see-through evidence disagrees with the request or the outcome, or
    /// a hidden piece has a nonzero count.
    IncorrectEvidence {
        /// The fixture.
        case: &'static str,
    },
    /// A base fixture flagged all opaque was not drawn as without flags.
    NotNeutral {
        /// The base fixture.
        case: &'static str,
    },
}

/// Result of [`SeeThroughConformanceSuite`].
#[derive(Debug, Clone, PartialEq)]
#[non_exhaustive]
pub struct SeeThroughConformanceReport {
    /// Every failure found.
    pub failures: Vec<SeeThroughConformanceFailure>,
    /// Whether the provider claims see-through triangles, as
    /// `supports_occlusion(Occlusion::SeeThrough)` says.
    pub claimed: bool,
}

impl SeeThroughConformanceReport {
    /// Whether no check failed.
    pub fn is_success(&self) -> bool {
        self.failures.is_empty()
    }
}

impl fmt::Display for SeeThroughConformanceReport {
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

/// The see-through conformance suite for `MeshHiddenLines` and
/// `MeshVisibleFaces`.
pub struct SeeThroughConformanceSuite;

impl SeeThroughConformanceSuite {
    /// Check a hidden-line provider's see-through claim.
    pub fn run(provider: &dyn MeshHiddenLines) -> SeeThroughConformanceReport {
        let claimed = provider.supports_occlusion(Occlusion::SeeThrough);
        let mut report = SeeThroughConformanceReport {
            failures: Vec::new(),
            claimed,
        };
        if claimed {
            report.failures = claimed_failures(provider);
            return report;
        }
        let case = first_case();
        let result = draw(provider, &case, &case.meshes, &options());
        record_unclaimed(&mut report, "hidden lines", result.map(|_| ()));
        report
    }

    /// Check a visible-face provider's see-through claim.
    pub fn run_faces(provider: &dyn MeshVisibleFaces) -> SeeThroughConformanceReport {
        let claimed = provider.supports_occlusion(Occlusion::SeeThrough);
        let mut report = SeeThroughConformanceReport {
            failures: Vec::new(),
            claimed,
        };
        const OPERATION: &str = "visible faces";
        if claimed {
            report.failures = claimed_face_failures(provider);
            return report;
        }
        let case = first_case();
        let meshes: Vec<HiddenLineMesh<'_>> = case.meshes.iter().map(mesh_of).collect();
        let request = VisibleFaceRequest::new(
            &meshes,
            case.view,
            crate::face_fixtures::limits().with_max_behind_surfaces(16),
        );
        let result = provider.visible_faces(&request, &options());
        record_unclaimed(&mut report, OPERATION, result.map(|_| ()));
        report
    }
}

/// The pane over a floor: what an unclaiming provider must refuse.
fn first_case() -> Case {
    see_through_cases()
        .into_iter()
        .next()
        .expect("a see-through case")
}

/// The failures of a provider that claims see-through triangles. The
/// provider is called directly, as in the base suite, so it must refuse a
/// short occlusion slice itself.
fn claimed_failures(provider: &dyn MeshHiddenLines) -> Vec<SeeThroughConformanceFailure> {
    let options = options();
    let mut base = ConformanceReport::default();
    let mut failures = Vec::new();
    let perspective = provider.supports_projection(ViewProjection::Perspective);
    for case in see_through_cases() {
        if case.view.projection == ViewProjection::Perspective && !perspective {
            continue;
        }
        let first = draw(provider, &case, &case.meshes, &options);
        let second = draw(provider, &case, &case.meshes, &options);
        match (first, second) {
            (Ok(a), Ok(b)) => {
                check_outcome(&mut base, &case, &a);
                check_seen_through(&mut failures, &case, &a);
                if a != b {
                    base.failures
                        .push(ConformanceFailure::NonDeterministic { case: case.name });
                }
                let flipped: Vec<Fixture> = case.meshes.iter().map(Fixture::flipped).collect();
                match draw(provider, &case, &flipped, &options) {
                    Ok(outcome) if same_drawing(&outcome, &a) => {}
                    Ok(_) => base
                        .failures
                        .push(ConformanceFailure::OrientationDependent { case: case.name }),
                    Err(error) => base.failures.push(ConformanceFailure::ReturnedError {
                        case: case.name,
                        error,
                    }),
                }
            }
            (Err(error), _) | (_, Err(error)) => {
                base.failures.push(ConformanceFailure::ReturnedError {
                    case: case.name,
                    error,
                })
            }
        }
    }
    for refusal in see_through_refusals() {
        let meshes: Vec<HiddenLineMesh<'_>> = refusal.meshes.iter().map(mesh_of).collect();
        let request = HiddenLineRequest::new(&meshes, refusal.view, refusal.edges, refusal.limits);
        match provider.hidden_lines(&request, &options) {
            Ok(_) => base
                .failures
                .push(ConformanceFailure::MissingRefusal { case: refusal.name }),
            Err(error) if !(refusal.accepts)(&error.error) => {
                base.failures.push(ConformanceFailure::ReturnedError {
                    case: refusal.name,
                    error,
                })
            }
            Err(error) if error.mesh != refusal.mesh => {
                base.failures.push(ConformanceFailure::WrongRefusalMesh {
                    case: refusal.name,
                    expected: refusal.mesh,
                    actual: error.mesh,
                })
            }
            Err(_) => {}
        }
    }
    for case in cases()
        .into_iter()
        .filter(|case| NEUTRAL_CASES.contains(&case.name))
    {
        let opaque: Vec<Fixture> = case
            .meshes
            .iter()
            .map(|f| Fixture {
                id: f.id,
                role: f.role,
                mesh: f.mesh.clone(),
                groups: f.groups.clone(),
                planar: f.planar,
                occlusion: vec![Occlusion::Opaque; f.mesh.triangle_count()],
            })
            .collect();
        match (
            draw(provider, &case, &case.meshes, &options),
            draw(provider, &case, &opaque, &options),
        ) {
            (Ok(plain), Ok(flagged)) if same_drawing(&plain, &flagged) => {}
            (Ok(_), Ok(_)) => {
                failures.push(SeeThroughConformanceFailure::NotNeutral { case: case.name })
            }
            (Err(error), _) | (_, Err(error)) => {
                base.failures.push(ConformanceFailure::ReturnedError {
                    case: case.name,
                    error,
                })
            }
        }
    }
    base.failures
        .into_iter()
        .map(SeeThroughConformanceFailure::HiddenLines)
        .chain(failures)
        .collect()
}

/// The failures of a visible-face provider that claims see-through
/// triangles. The provider is called directly, so it must refuse a mixed
/// group, a short occlusion slice and a stack over the limit itself.
fn claimed_face_failures(provider: &dyn MeshVisibleFaces) -> Vec<SeeThroughConformanceFailure> {
    let options = options();
    let mut base = FaceConformanceReport::default();
    let mut failures = Vec::new();
    for case in see_through_face_cases() {
        let scene = Scene {
            name: case.name,
            meshes: &case.meshes,
            view: case.view,
        };
        let ask = || {
            answer(
                provider,
                &case.meshes,
                case.view,
                case.region,
                case.limits,
                &options,
            )
        };
        match (ask(), ask()) {
            (Ok(a), Ok(b)) => {
                check_shape(&mut base, &scene, &a);
                if a.regions != case.regions {
                    base.failures.push(FaceConformanceFailure::WrongRegions {
                        case: case.name,
                        expected: case.regions.clone(),
                        actual: a.regions.clone(),
                    });
                }
                check_expected_evidence(&mut base, case.name, &case.evidence, &a);
                let flags: usize = case.meshes.iter().map(see_through_flags).sum();
                let behind: usize = a.regions.iter().map(|r| r.behind.len()).sum();
                if a.evidence.see_through_triangles != flags || a.evidence.behind_surfaces != behind
                {
                    failures
                        .push(SeeThroughConformanceFailure::IncorrectEvidence { case: case.name });
                }
                if a != b {
                    base.failures
                        .push(FaceConformanceFailure::NonDeterministic { case: case.name });
                }
                check_orientation(
                    provider,
                    &mut base,
                    &scene,
                    case.region,
                    case.limits,
                    &a,
                    &options,
                );
            }
            (Err(error), _) | (_, Err(error)) => {
                base.failures.push(FaceConformanceFailure::ReturnedError {
                    case: case.name,
                    error,
                })
            }
        }
    }
    for refusal in see_through_face_refusals() {
        match answer(
            provider,
            &refusal.meshes,
            refusal.view,
            None,
            refusal.limits,
            &options,
        ) {
            Ok(_) => base
                .failures
                .push(FaceConformanceFailure::MissingRefusal { case: refusal.name }),
            Err(error) if !(refusal.accepts)(&error.error) => {
                base.failures.push(FaceConformanceFailure::ReturnedError {
                    case: refusal.name,
                    error,
                })
            }
            Err(error) if error.mesh != refusal.mesh => {
                base.failures
                    .push(FaceConformanceFailure::WrongRefusalMesh {
                        case: refusal.name,
                        expected: refusal.mesh,
                        actual: error.mesh,
                    })
            }
            Err(_) => {}
        }
    }
    for case in face_cases()
        .into_iter()
        .filter(|case| NEUTRAL_FACE_CASES.contains(&case.name))
    {
        let opaque: Vec<FaceFixture> = case
            .meshes
            .iter()
            .map(|f| FaceFixture {
                id: f.id,
                role: f.role,
                mesh: f.mesh.clone(),
                groups: f.groups.clone(),
                planar: f.planar,
                occlusion: vec![Occlusion::Opaque; f.mesh.triangle_count()],
            })
            .collect();
        let ask = |meshes: &[FaceFixture]| {
            answer(
                provider,
                meshes,
                case.view,
                case.region,
                case.limits,
                &options,
            )
        };
        match (ask(&case.meshes), ask(&opaque)) {
            (Ok(plain), Ok(flagged)) if same_result(&plain, &flagged) => {}
            (Ok(_), Ok(_)) => {
                failures.push(SeeThroughConformanceFailure::NotNeutral { case: case.name })
            }
            (Err(error), _) | (_, Err(error)) => {
                base.failures.push(FaceConformanceFailure::ReturnedError {
                    case: case.name,
                    error,
                })
            }
        }
    }
    base.failures
        .into_iter()
        .map(SeeThroughConformanceFailure::VisibleFaces)
        .chain(failures)
        .collect()
}

fn see_through_flags(fixture: &FaceFixture) -> usize {
    fixture
        .occlusion
        .iter()
        .filter(|&&o| o == Occlusion::SeeThrough)
        .count()
}

/// The `seen_through` counts and the see-through evidence of one drawing.
/// Geometry and visibility are the base checks' (`check_outcome`).
fn check_seen_through(
    failures: &mut Vec<SeeThroughConformanceFailure>,
    case: &Case,
    outcome: &HiddenLineOutcome,
) {
    let name = case.name;
    let flags: usize = case
        .meshes
        .iter()
        .map(|f| {
            f.occlusion
                .iter()
                .filter(|&&o| o == Occlusion::SeeThrough)
                .count()
        })
        .sum();
    let counted = outcome
        .segments
        .iter()
        .filter(|s| s.seen_through > 0)
        .count();
    let hidden_counted = outcome
        .segments
        .iter()
        .any(|s| s.visibility == Visibility::Hidden && s.seen_through > 0);
    if outcome.evidence.see_through_triangles != flags
        || outcome.evidence.seen_through_segments != counted
        || hidden_counted
    {
        failures.push(SeeThroughConformanceFailure::IncorrectEvidence { case: name });
    }

    if let Some(rows) = &case.expect.rows {
        let mut keys: Vec<(EdgeKind, Visibility, u32)> = Vec::new();
        let returned = outcome
            .segments
            .iter()
            .map(|s| (s.kind, s.visibility, s.seen_through));
        let expected = rows.iter().map(|r| (r.kind, r.visibility, r.seen_through));
        for key in expected.chain(returned) {
            if !keys.contains(&key) {
                keys.push(key);
            }
        }
        for (kind, visibility, seen_through) in keys {
            let (expected_count, expected_length) = rows
                .iter()
                .filter(|r| {
                    (r.kind, r.visibility, r.seen_through) == (kind, visibility, seen_through)
                })
                .fold((0, 0.0), |(n, sum), r| (n + r.count, sum + r.length));
            let (actual_count, actual_length) = outcome
                .segments
                .iter()
                .filter(|s| {
                    (s.kind, s.visibility, s.seen_through) == (kind, visibility, seen_through)
                })
                .fold((0, 0.0), |(n, sum), s| {
                    (n + 1, sum + (s.end - s.start).length())
                });
            if expected_count != actual_count
                || (expected_length - actual_length).abs() > LENGTH_TOLERANCE
            {
                failures.push(SeeThroughConformanceFailure::WrongSeenThroughRow {
                    case: name,
                    kind,
                    visibility,
                    seen_through,
                    expected_count,
                    actual_count,
                    expected_length,
                    actual_length,
                });
            }
        }
    }

    for wanted in &case.expect.edges {
        let expected: Vec<u32> = wanted.pieces.iter().map(|p| p.seen_through).collect();
        let actual: Vec<u32> = outcome
            .segments
            .iter()
            .filter(|s| s.source == wanted.edge)
            .map(|s| s.seen_through)
            .collect();
        if expected != actual {
            failures.push(SeeThroughConformanceFailure::WrongSeenThrough {
                case: name,
                edge: wanted.edge,
                expected,
                actual,
            });
        }
    }
}

fn options() -> ExecutionOptions {
    ExecutionOptions::new(Tolerance::new(1e-9, 1e-9).expect("valid tolerance"))
}

fn record_unclaimed(
    report: &mut SeeThroughConformanceReport,
    operation: &'static str,
    result: Result<(), HiddenLineError>,
) {
    match result {
        Ok(()) => report
            .failures
            .push(SeeThroughConformanceFailure::UnclaimedButDrawn { operation }),
        Err(HiddenLineError {
            error: GeomError::Unsupported { .. } | GeomError::UnsupportedInput { .. },
            ..
        }) => {}
        Err(error) => report
            .failures
            .push(SeeThroughConformanceFailure::UnclaimedWrongRefusal { operation, error }),
    }
}

#[cfg(test)]
mod tests {
    use axiolid_contracts::{Backend, BackendDescriptor, BackendId, ExecutionTarget, Operation};

    use super::*;
    use crate::{HiddenLineEvidence, VisibleFaceEvidence, VisibleFaceOutcome};

    /// A provider whose claim and answer to a see-through request are chosen.
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
                BackendId::new("see-through-probe"),
                ExecutionTarget::PortableCpu,
            )
        }
    }

    impl Probe {
        fn refusal(&self, operation: Operation) -> Option<HiddenLineError> {
            let backend = BackendId::new("see-through-probe");
            match self.answer {
                Answer::Unsupported => Some(GeomError::Unsupported { backend, operation }.into()),
                Answer::UnsupportedInput => Some(
                    GeomError::UnsupportedInput {
                        backend,
                        operation,
                        input: "see-through triangles",
                    }
                    .into(),
                ),
                Answer::Invalid => Some(GeomError::InvalidInput("no".into()).into()),
                Answer::Draw => None,
            }
        }
    }

    impl MeshHiddenLines for Probe {
        fn supports_occlusion(&self, occlusion: Occlusion) -> bool {
            occlusion == Occlusion::Opaque || self.claims
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
                    HiddenLineEvidence::input_mesh(0, 0),
                )),
            }
        }
    }

    impl MeshVisibleFaces for Probe {
        fn supports_occlusion(&self, occlusion: Occlusion) -> bool {
            occlusion == Occlusion::Opaque || self.claims
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
                    VisibleFaceEvidence::input_mesh(0, 0),
                )),
            }
        }
    }

    /// A provider that keeps every default.
    #[derive(Debug)]
    struct Plain;

    impl Backend for Plain {
        fn descriptor(&self) -> BackendDescriptor {
            BackendDescriptor::new(BackendId::new("plain"), ExecutionTarget::PortableCpu)
        }
    }

    impl MeshHiddenLines for Plain {
        fn hidden_lines(
            &self,
            _request: &HiddenLineRequest<'_>,
            _options: &ExecutionOptions,
        ) -> Result<HiddenLineOutcome, HiddenLineError> {
            Err(GeomError::InvalidInput("not called".into()).into())
        }
    }

    impl MeshVisibleFaces for Plain {
        fn visible_faces(
            &self,
            _request: &VisibleFaceRequest<'_>,
            _options: &ExecutionOptions,
        ) -> Result<VisibleFaceOutcome, HiddenLineError> {
            Err(GeomError::InvalidInput("not called".into()).into())
        }
    }

    #[test]
    fn the_default_claims_opaque_triangles_only() {
        let lines: &dyn MeshHiddenLines = &Plain;
        let faces: &dyn MeshVisibleFaces = &Plain;
        assert!(lines.supports_occlusion(Occlusion::Opaque));
        assert!(!lines.supports_occlusion(Occlusion::SeeThrough));
        assert!(faces.supports_occlusion(Occlusion::Opaque));
        assert!(!faces.supports_occlusion(Occlusion::SeeThrough));
    }

    fn both(probe: &Probe) -> [SeeThroughConformanceReport; 2] {
        [
            SeeThroughConformanceSuite::run(probe),
            SeeThroughConformanceSuite::run_faces(probe),
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
            [SeeThroughConformanceFailure::UnclaimedButDrawn {
                operation: "hidden lines"
            }]
        );
        assert_eq!(
            faces.failures,
            [SeeThroughConformanceFailure::UnclaimedButDrawn {
                operation: "visible faces"
            }]
        );
        for report in both(&Probe {
            claims: false,
            answer: Answer::Invalid,
        }) {
            assert!(matches!(
                report.failures[..],
                [SeeThroughConformanceFailure::UnclaimedWrongRefusal { .. }]
            ));
        }
    }

    #[test]
    fn a_claim_that_ignores_the_flags_fails_the_fixtures() {
        let lines = SeeThroughConformanceSuite::run(&Probe {
            claims: true,
            answer: Answer::Draw,
        });
        assert!(lines.claimed);
        // The probe draws nothing: every fixture fails, the short slice is
        // drawn instead of refused, and its evidence counts no flag.
        assert!(lines.failures.iter().any(|f| matches!(
            f,
            SeeThroughConformanceFailure::HiddenLines(ConformanceFailure::WrongEdge {
                case: "pane over a floor",
                ..
            })
        )));
        assert!(lines.failures.iter().any(|f| matches!(
            f,
            SeeThroughConformanceFailure::HiddenLines(ConformanceFailure::MissingRefusal {
                case: "occlusion flags shorter than the triangles"
            })
        )));
        assert!(lines.failures.iter().any(|f| matches!(
            f,
            SeeThroughConformanceFailure::WrongSeenThrough {
                case: "edge under a pane's diagonal",
                ..
            }
        )));
        assert!(lines
            .failures
            .contains(&SeeThroughConformanceFailure::IncorrectEvidence {
                case: "pane over a floor"
            }));
        // The perspective case is drawn only for a perspective claim.
        assert!(!lines.to_string().contains("perspective see-through slab"));
    }

    #[test]
    fn a_visible_face_claim_that_ignores_the_flags_fails_the_fixtures() {
        let faces = SeeThroughConformanceSuite::run_faces(&Probe {
            claims: true,
            answer: Answer::Draw,
        });
        assert!(faces.claimed);
        // The probe answers nothing: every case has the wrong regions, and
        // every refusal is answered instead.
        for case in ["pane over a box", "a stack stops at an opaque sheet"] {
            assert!(
                faces.failures.iter().any(|f| matches!(
                    f,
                    SeeThroughConformanceFailure::VisibleFaces(
                        FaceConformanceFailure::WrongRegions { case: c, .. }
                    ) if *c == case
                )),
                "{case}: {faces}"
            );
        }
        for case in [
            "a face group mixing see-through and opaque triangles",
            "behind surfaces over the limit",
        ] {
            assert!(faces.failures.iter().any(|f| matches!(
                f,
                SeeThroughConformanceFailure::VisibleFaces(
                    FaceConformanceFailure::MissingRefusal { case: c }
                ) if *c == case
            )));
        }
        assert!(faces
            .failures
            .contains(&SeeThroughConformanceFailure::IncorrectEvidence {
                case: "pane over a box"
            }));
    }
}
