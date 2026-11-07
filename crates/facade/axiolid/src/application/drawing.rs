//! Hidden-line drawings and visible faces on [`Application`] (ADR 0105).
//!
//! The exact scalar providers are stateless and their conformance suites
//! deterministic, so each suite runs at most once per process, on the first
//! call that reaches the provider ([`verified`]). `Application::portable()`
//! therefore stays cheap (the face suite alone takes about half a second,
//! and the C ABI builds an application per context). A provider that fails
//! its suite refuses every call with
//! [`ApplicationErrorSource::ProviderConformance`]; there is no output.

use std::sync::OnceLock;

use axiolid_contracts::{
    capability_ids, BackendId, CapabilityDescriptor, Exactness, ExecutionOptions,
    IntegrationDescriptor, Operation, Representation,
};
use axiolid_hidden_lines::{ScalarHiddenLines, ScalarVisibleFaces};
use axiolid_mesh_hidden_lines_contract::conformance::{
    faces::FaceConformanceSuite, ConformanceSuite,
};
use axiolid_mesh_hidden_lines_contract::{
    HiddenLineError, HiddenLineOutcome, HiddenLineRequest, RefusedMesh, VisibleFaceOutcome,
    VisibleFaceRequest,
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

/// The visible-face provider's verdict, tied to its suite.
fn visible_faces_verdict() -> Result<(), &'static ProviderConformanceReport> {
    static CELL: Verdict = OnceLock::new();
    verified(&CELL, || {
        ProviderConformanceReport::MeshVisibleFaces(Box::new(FaceConformanceSuite::run(
            &ScalarVisibleFaces::new(),
        )))
    })
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
