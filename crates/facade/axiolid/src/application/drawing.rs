//! Hidden-line drawings and visible faces on [`Application`] (ADR 0105).
//!
//! The exact scalar providers are stateless and their conformance suites
//! deterministic, so each suite runs at most once per process, on the first
//! call that reaches the provider ([`verified`]). `Application::portable()`
//! therefore stays cheap (the face suite alone takes about half a second,
//! and the C ABI builds an application per context). A provider that fails
//! its suite refuses every call with
//! [`ApplicationErrorSource::ProviderConformance`]; there is no output.
//!
//! The hidden-line provider also claims perspective views (ADR 0107). That
//! claim has its own suite, run once and only by a call whose view is
//! perspective, so an orthographic call never pays for it; a failure is
//! folded into the same refusal as the case `perspective claim`, as the
//! dispatch registry does. Both providers claim see-through triangles (ADR
//! 0111) in the same way: that suite runs once and only by a call with a
//! see-through triangle, and a failure is the case `see-through claim`.

use std::sync::OnceLock;

use axiolid_contracts::{
    capability_ids, BackendId, CapabilityDescriptor, Exactness, ExecutionOptions, GeomError,
    IntegrationDescriptor, Operation, Representation,
};
use axiolid_hidden_lines::{ScalarHiddenLines, ScalarVisibleFaces};
use axiolid_mesh_hidden_lines_contract::conformance::{
    faces::{FaceConformanceFailure, FaceConformanceReport, FaceConformanceSuite},
    perspective::{PerspectiveConformanceReport, PerspectiveConformanceSuite},
    see_through::{SeeThroughConformanceReport, SeeThroughConformanceSuite},
    ConformanceFailure, ConformanceReport, ConformanceSuite,
};
use axiolid_mesh_hidden_lines_contract::{
    HiddenLineError, HiddenLineMesh, HiddenLineOutcome, HiddenLineRequest, RefusedMesh,
    ViewProjection, VisibleFaceOutcome, VisibleFaceRequest,
};

use super::{
    Application, ApplicationBuilder, ApplicationError, ApplicationErrorSource, CallContext,
    ProviderConformanceReport, APPLICATION_ID, MESH_INPUT,
};

type Verdict = OnceLock<Result<(), ProviderConformanceReport>>;

/// The hidden-line provider's verdict, tied to its suite.
fn hidden_lines_verdict() -> Result<(), &'static ProviderConformanceReport> {
    static CELL: Verdict = OnceLock::new();
    verified(&CELL, || {
        ProviderConformanceReport::MeshHiddenLines(Box::new(ConformanceSuite::run(
            &ScalarHiddenLines::new(),
        )))
    })
}

/// The hidden-line provider's verdict on its perspective claim.
fn perspective_hidden_lines_verdict() -> Result<(), &'static ProviderConformanceReport> {
    static CELL: Verdict = OnceLock::new();
    verified(&CELL, || {
        perspective_report(PerspectiveConformanceSuite::run(&ScalarHiddenLines::new()))
    })
}

/// A perspective suite report as a hidden-line report: a failure is the
/// case `perspective claim`, as `register_conformant` of the dispatch
/// registry reports it.
fn perspective_report(report: PerspectiveConformanceReport) -> ProviderConformanceReport {
    let mut folded = ConformanceReport::default();
    if !report.is_success() {
        folded.failures.push(ConformanceFailure::ReturnedError {
            case: "perspective claim",
            error: HiddenLineError::request(GeomError::InvalidInput(format!(
                "perspective suite: {report}"
            ))),
        });
    }
    ProviderConformanceReport::MeshHiddenLines(Box::new(folded))
}

/// The visible-face provider's verdict, tied to its suite.
fn visible_faces_verdict() -> Result<(), &'static ProviderConformanceReport> {
    static CELL: Verdict = OnceLock::new();
    verified(&CELL, || {
        ProviderConformanceReport::MeshVisibleFaces(Box::new(FaceConformanceSuite::run(
            &ScalarVisibleFaces::new(),
        )))
    })
}

/// Whether a request holds a see-through triangle (ADR 0111).
fn has_see_through(meshes: &[HiddenLineMesh<'_>]) -> bool {
    meshes.iter().any(|entry| entry.occlusion.has_see_through())
}

/// The hidden-line provider's verdict on its see-through claim.
fn see_through_hidden_lines_verdict() -> Result<(), &'static ProviderConformanceReport> {
    static CELL: Verdict = OnceLock::new();
    verified(&CELL, || {
        see_through_report(SeeThroughConformanceSuite::run(&ScalarHiddenLines::new()))
    })
}

/// The visible-face provider's verdict on its see-through claim.
fn see_through_visible_faces_verdict() -> Result<(), &'static ProviderConformanceReport> {
    static CELL: Verdict = OnceLock::new();
    verified(&CELL, || {
        see_through_faces_report(SeeThroughConformanceSuite::run_faces(
            &ScalarVisibleFaces::new(),
        ))
    })
}

/// The refusal a failed see-through suite is folded into: the case
/// `see-through claim`, as `register_conformant` reports it.
fn see_through_error(report: &SeeThroughConformanceReport) -> HiddenLineError {
    HiddenLineError::request(GeomError::InvalidInput(format!(
        "see-through suite: {report}"
    )))
}

/// A see-through suite report as a hidden-line report.
fn see_through_report(report: SeeThroughConformanceReport) -> ProviderConformanceReport {
    let mut folded = ConformanceReport::default();
    folded.failures.extend(
        (!report.is_success()).then(|| ConformanceFailure::ReturnedError {
            case: "see-through claim",
            error: see_through_error(&report),
        }),
    );
    ProviderConformanceReport::MeshHiddenLines(Box::new(folded))
}

/// A see-through suite report as a visible-face report.
fn see_through_faces_report(report: SeeThroughConformanceReport) -> ProviderConformanceReport {
    let mut folded = FaceConformanceReport::default();
    folded.failures.extend(
        (!report.is_success()).then(|| FaceConformanceFailure::ReturnedError {
            case: "see-through claim",
            error: see_through_error(&report),
        }),
    );
    ProviderConformanceReport::MeshVisibleFaces(Box::new(folded))
}

/// Run `run` the first time `cell` is asked and keep its verdict.
///
/// A report that is not a success is kept as the refusal reason; every later
/// call returns it without running the suite again.
fn verified(
    cell: &Verdict,
    run: impl FnOnce() -> ProviderConformanceReport,
) -> Result<(), &ProviderConformanceReport> {
    cell.get_or_init(|| {
        let report = run();
        if report.is_success() {
            Ok(())
        } else {
            Err(report)
        }
    })
    .as_ref()
    .map(|&()| ())
}

impl ApplicationBuilder {
    /// Select the exact scalar hidden-line provider.
    ///
    /// Its conformance suite runs on first use, not here (see the module
    /// docs).
    #[must_use]
    pub fn with_portable_hidden_lines(mut self) -> Self {
        self.hidden_lines.register(0, ScalarHiddenLines::new());
        self.hidden_lines_provider = Some(ScalarHiddenLines::ID);
        self
    }

    /// Select the exact scalar visible-face provider.
    ///
    /// Its conformance suite runs on first use, not here (see the module
    /// docs).
    #[must_use]
    pub fn with_portable_visible_faces(mut self) -> Self {
        self.visible_faces.register(0, ScalarVisibleFaces::new());
        self.visible_faces_provider = Some(ScalarVisibleFaces::ID);
        self
    }
}

/// Advertise the drawing providers the builder selected.
pub(super) fn push_capabilities(
    descriptor: &mut IntegrationDescriptor,
    hidden_lines: Option<BackendId>,
    visible_faces: Option<BackendId>,
) {
    for (id, operation, provider) in [
        (
            capability_ids::MESH_HIDDEN_LINES,
            Operation::MeshHiddenLines,
            hidden_lines,
        ),
        (
            capability_ids::MESH_VISIBLE_FACES,
            Operation::MeshVisibleFaces,
            visible_faces,
        ),
    ] {
        if let Some(provider) = provider {
            descriptor.capabilities.push(CapabilityDescriptor {
                id,
                operation,
                provider,
                required_feature: "application",
                inputs: MESH_INPUT,
                output: Representation::Profile2d,
                exactness: Exactness::ToleranceBounded,
                deterministic: true,
            });
        }
    }
}

impl Application {
    /// Draw the edges of `request.meshes` with hidden-line removal.
    ///
    /// A refusal that names a mesh keeps it in
    /// [`ApplicationErrorSource::HiddenLines`].
    pub fn hidden_lines(
        &self,
        request: &HiddenLineRequest<'_>,
        options: &ExecutionOptions,
    ) -> Result<HiddenLineOutcome, ApplicationError> {
        let provider = self.hidden_lines_provider.unwrap_or(APPLICATION_ID);
        options.check_cancelled().map_err(|source| {
            drawing_error(Operation::MeshHiddenLines, provider, options, source.into())
        })?;
        if self.hidden_lines_provider.is_some() {
            hidden_lines_verdict().map_err(|report| {
                conformance_error(Operation::MeshHiddenLines, provider, options, report)
            })?;
            if request.view.projection == ViewProjection::Perspective {
                perspective_hidden_lines_verdict().map_err(|report| {
                    conformance_error(Operation::MeshHiddenLines, provider, options, report)
                })?;
            }
            if has_see_through(request.meshes) {
                see_through_hidden_lines_verdict().map_err(|report| {
                    conformance_error(Operation::MeshHiddenLines, provider, options, report)
                })?;
            }
        }
        self.hidden_lines
            .hidden_lines(request, options)
            .map_err(|source| drawing_error(Operation::MeshHiddenLines, provider, options, source))
    }

    /// Every input mesh [`hidden_lines`](Self::hidden_lines) would refuse,
    /// in request order, without drawing.
    pub fn audit_hidden_lines(
        &self,
        request: &HiddenLineRequest<'_>,
        options: &ExecutionOptions,
    ) -> Result<Vec<RefusedMesh>, ApplicationError> {
        self.hidden_lines.audit(request, options).map_err(|source| {
            drawing_error(
                Operation::MeshHiddenLines,
                self.hidden_lines_provider.unwrap_or(APPLICATION_ID),
                options,
                source,
            )
        })
    }

    /// Visible face regions of `request.meshes` in the view's image plane
    /// (ADR 0094).
    pub fn visible_faces(
        &self,
        request: &VisibleFaceRequest<'_>,
        options: &ExecutionOptions,
    ) -> Result<VisibleFaceOutcome, ApplicationError> {
        let provider = self.visible_faces_provider.unwrap_or(APPLICATION_ID);
        options.check_cancelled().map_err(|source| {
            drawing_error(
                Operation::MeshVisibleFaces,
                provider,
                options,
                source.into(),
            )
        })?;
        if self.visible_faces_provider.is_some() {
            visible_faces_verdict().map_err(|report| {
                conformance_error(Operation::MeshVisibleFaces, provider, options, report)
            })?;
            if has_see_through(request.meshes) {
                see_through_visible_faces_verdict().map_err(|report| {
                    conformance_error(Operation::MeshVisibleFaces, provider, options, report)
                })?;
            }
        }
        self.visible_faces
            .visible_faces(request, options)
            .map_err(|source| drawing_error(Operation::MeshVisibleFaces, provider, options, source))
    }
}

fn drawing_error(
    operation: Operation,
    provider: BackendId,
    options: &ExecutionOptions,
    source: HiddenLineError,
) -> ApplicationError {
    ApplicationError {
        context: CallContext {
            operation,
            provider,
            tolerance: options.tolerance(),
        },
        source: Box::new(ApplicationErrorSource::HiddenLines(source)),
    }
}

fn conformance_error(
    operation: Operation,
    provider: BackendId,
    options: &ExecutionOptions,
    report: &ProviderConformanceReport,
) -> ApplicationError {
    ApplicationError {
        context: CallContext {
            operation,
            provider,
            tolerance: options.tolerance(),
        },
        source: Box::new(ApplicationErrorSource::ProviderConformance(report.clone())),
    }
}

#[cfg(test)]
mod tests {
    use std::cell::Cell;

    use axiolid_mesh_hidden_lines_contract::conformance::{ConformanceFailure, ConformanceReport};

    use super::*;

    fn failing() -> ProviderConformanceReport {
        let mut report = ConformanceReport::default();
        report
            .failures
            .push(ConformanceFailure::NonDeterministic { case: "memo" });
        ProviderConformanceReport::MeshHiddenLines(Box::new(report))
    }

    #[test]
    fn a_failing_verdict_is_returned_and_the_suite_runs_once() {
        static CELL: Verdict = OnceLock::new();
        let runs = Cell::new(0);
        for _ in 0..3 {
            let verdict = verified(&CELL, || {
                runs.set(runs.get() + 1);
                failing()
            });
            let report = verdict.expect_err("the failing report is the refusal");
            assert!(!report.is_success());
        }
        assert_eq!(runs.get(), 1);
    }

    #[test]
    fn a_failed_perspective_suite_is_folded_into_a_hidden_line_report() {
        use axiolid_mesh_hidden_lines_contract::conformance::perspective::PerspectiveConformanceFailure;

        let mut failed = PerspectiveConformanceSuite::run(&ScalarHiddenLines::new());
        assert!(failed.is_success(), "the shipped provider passes: {failed}");
        failed
            .failures
            .push(PerspectiveConformanceFailure::ClaimUnverifiable {
                operation: "hidden lines",
            });
        let folded = perspective_report(failed);
        assert!(!folded.is_success());
        let ProviderConformanceReport::MeshHiddenLines(report) = folded else {
            panic!("a hidden-line report");
        };
        assert!(matches!(
            report.failures[..],
            [ConformanceFailure::ReturnedError {
                case: "perspective claim",
                ..
            }]
        ));
        let passed =
            perspective_report(PerspectiveConformanceSuite::run(&ScalarHiddenLines::new()));
        assert!(passed.is_success());
    }

    #[test]
    fn a_failed_see_through_suite_is_folded_into_each_report() {
        use axiolid_mesh_hidden_lines_contract::conformance::see_through::SeeThroughConformanceFailure;

        let failure = SeeThroughConformanceFailure::NotNeutral { case: "probe" };
        let mut lines = SeeThroughConformanceSuite::run(&ScalarHiddenLines::new());
        assert!(lines.is_success(), "the shipped provider passes: {lines}");
        assert!(see_through_report(lines.clone()).is_success());
        lines.failures.push(failure.clone());
        let ProviderConformanceReport::MeshHiddenLines(report) = see_through_report(lines) else {
            panic!("a hidden-line report");
        };
        assert!(matches!(
            report.failures[..],
            [ConformanceFailure::ReturnedError {
                case: "see-through claim",
                ..
            }]
        ));

        let mut faces = SeeThroughConformanceSuite::run_faces(&ScalarVisibleFaces::new());
        assert!(faces.is_success(), "the shipped provider passes: {faces}");
        assert!(see_through_faces_report(faces.clone()).is_success());
        faces.failures.push(failure);
        let ProviderConformanceReport::MeshVisibleFaces(report) = see_through_faces_report(faces)
        else {
            panic!("a visible-face report");
        };
        assert!(matches!(
            report.failures[..],
            [FaceConformanceFailure::ReturnedError {
                case: "see-through claim",
                ..
            }]
        ));
    }

    #[test]
    fn a_passing_verdict_is_kept() {
        static CELL: Verdict = OnceLock::new();
        let runs = Cell::new(0);
        for _ in 0..2 {
            verified(&CELL, || {
                runs.set(runs.get() + 1);
                ProviderConformanceReport::MeshHiddenLines(Box::default())
            })
            .expect("an empty report passes");
        }
        assert_eq!(runs.get(), 1);
    }
}
