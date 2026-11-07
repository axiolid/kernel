//! Perspective routing and refusals (ADR 0107). The hidden-line providers
//! claim perspective and pass its fixtures; the visible-face providers refuse
//! it by name when called directly, and the registries find no visible-face
//! provider that supports the projection.

mod common;

use axiolid_contracts::{ExecutionOptions, GeomError, Operation};
use axiolid_core::Tolerance;
use axiolid_dispatch::{MeshHiddenLinesRegistry, MeshVisibleFacesRegistry};
use axiolid_hidden_lines::{
    ScalarHiddenLines, ScalarVisibleFaces, ScanHiddenLines, ScanVisibleFaces, UnprunedVisibleFaces,
};
use axiolid_mesh_hidden_lines_contract::conformance::perspective::PerspectiveConformanceSuite;
use axiolid_mesh_hidden_lines_contract::{
    EdgeSelection, HiddenLineError, HiddenLineOutcome, HiddenLineRequest, MeshHiddenLines,
    MeshVisibleFaces, View, ViewProjection, VisibleFaceOutcome, VisibleFaceRequest,
};
use common::{box_mesh, plan_frame, solids, FACE_LIMITS, LIMITS};

fn options() -> ExecutionOptions {
    ExecutionOptions::new(Tolerance::new(1e-9, 1e-9).expect("tolerance"))
}

fn perspective() -> View {
    View::perspective(plan_frame(), -1.0).expect("valid")
}

fn is_named_refusal(error: &GeomError, operation: Operation) -> bool {
    matches!(
        error,
        GeomError::UnsupportedInput { operation: o, input: "perspective view", .. } if *o == operation
    )
}

#[test]
fn only_the_hidden_line_providers_support_perspective() {
    let lines: [&dyn MeshHiddenLines; 2] = [&ScalarHiddenLines, &ScanHiddenLines];
    let faces: [&dyn MeshVisibleFaces; 3] = [
        &ScalarVisibleFaces,
        &ScanVisibleFaces,
        &UnprunedVisibleFaces,
    ];
    for provider in lines {
        assert!(provider.supports_projection(ViewProjection::Orthographic));
        assert!(provider.supports_projection(ViewProjection::Perspective));
        let report = PerspectiveConformanceSuite::run(provider);
        assert!(report.claimed && report.is_success(), "{report}");
    }
    for provider in faces {
        assert!(provider.supports_projection(ViewProjection::Orthographic));
        assert!(!provider.supports_projection(ViewProjection::Perspective));
        let report = PerspectiveConformanceSuite::run_faces(provider);
        assert!(!report.claimed && report.is_success(), "{report}");
    }
}

/// A hidden-line provider that leaves the claim at its default.
#[derive(Debug)]
struct Plain;

impl axiolid_contracts::Backend for Plain {
    fn descriptor(&self) -> axiolid_contracts::BackendDescriptor {
        ScalarHiddenLines.descriptor()
    }
}

impl MeshHiddenLines for Plain {
    fn hidden_lines(
        &self,
        request: &HiddenLineRequest<'_>,
        options: &ExecutionOptions,
    ) -> Result<HiddenLineOutcome, HiddenLineError> {
        ScalarHiddenLines.hidden_lines(request, options)
    }
}

#[test]
fn a_provider_that_leaves_the_claim_at_its_default_claims_orthographic_only() {
    assert!(Plain.supports_projection(ViewProjection::Orthographic));
    assert!(!Plain.supports_projection(ViewProjection::Perspective));
}

#[test]
fn visible_face_providers_called_directly_refuse_a_perspective_view_by_name() {
    let boxes = [box_mesh([-1.0, -1.0, -4.0], [1.0, 1.0, -2.0])];
    let meshes = solids(&boxes);
    let faces = VisibleFaceRequest::new(&meshes, perspective(), FACE_LIMITS);
    for provider in [
        &ScalarVisibleFaces as &dyn MeshVisibleFaces,
        &ScanVisibleFaces,
        &UnprunedVisibleFaces,
    ] {
        let error = provider
            .visible_faces(&faces, &options())
            .expect_err("refused");
        assert!(
            is_named_refusal(&error.error, Operation::MeshVisibleFaces),
            "{error}"
        );
    }
}

#[test]
fn the_hidden_line_registry_draws_a_perspective_view_and_the_face_registry_declines() {
    let boxes = [box_mesh([-1.0, -1.0, -4.0], [1.0, 1.0, -2.0])];
    let meshes = solids(&boxes);

    let mut lines = MeshHiddenLinesRegistry::new();
    lines
        .register_perspective_conformant(0, ScalarHiddenLines::new())
        .unwrap_or_else(|report| panic!("{report}"));
    let request = HiddenLineRequest::new(&meshes, perspective(), EdgeSelection::new(None), LIMITS);
    let drawn = lines.hidden_lines(&request, &options()).expect("draws");
    assert_eq!(drawn.view, perspective());
    assert!(!drawn.segments.is_empty());

    let mut faces = MeshVisibleFacesRegistry::new();
    faces
        .register_perspective_conformant(0, ScalarVisibleFaces::new())
        .unwrap_or_else(|report| panic!("{report}"));
    let request = VisibleFaceRequest::new(&meshes, perspective(), FACE_LIMITS);
    let error = faces
        .visible_faces(&request, &options())
        .expect_err("unsupported");
    assert!(matches!(
        error.error,
        GeomError::Unsupported {
            operation: Operation::MeshVisibleFaces,
            ..
        }
    ));
}

#[test]
fn an_orthographic_view_is_unchanged_by_the_projection_field() {
    let boxes = [box_mesh([0.0, 0.0, 0.0], [1.0, 1.0, 1.0])];
    let meshes = solids(&boxes);
    let view = View::new(plan_frame());
    assert_eq!(view.projection, ViewProjection::Orthographic);
    let request = HiddenLineRequest::new(&meshes, view, EdgeSelection::new(None), LIMITS);
    ScalarHiddenLines
        .hidden_lines(&request, &options())
        .expect("draws");
}

/// An in-tree face provider that claims perspective, which it cannot draw.
#[derive(Debug)]
struct Claims;

impl axiolid_contracts::Backend for Claims {
    fn descriptor(&self) -> axiolid_contracts::BackendDescriptor {
        ScalarHiddenLines.descriptor()
    }
}

/// A hidden-line provider that claims perspective and ignores the
/// projection: it draws the orthographic picture.
#[derive(Debug)]
struct Ignores;

impl axiolid_contracts::Backend for Ignores {
    fn descriptor(&self) -> axiolid_contracts::BackendDescriptor {
        ScalarHiddenLines.descriptor()
    }
}

impl MeshHiddenLines for Ignores {
    fn supports_projection(&self, _projection: ViewProjection) -> bool {
        true
    }

    fn hidden_lines(
        &self,
        request: &HiddenLineRequest<'_>,
        options: &ExecutionOptions,
    ) -> Result<HiddenLineOutcome, HiddenLineError> {
        let mut orthographic = request.clone();
        orthographic.view = View::new(request.view.frame);
        ScalarHiddenLines.hidden_lines(&orthographic, options)
    }
}

impl MeshVisibleFaces for Claims {
    fn supports_projection(&self, _projection: ViewProjection) -> bool {
        true
    }

    fn visible_faces(
        &self,
        request: &VisibleFaceRequest<'_>,
        options: &ExecutionOptions,
    ) -> Result<VisibleFaceOutcome, HiddenLineError> {
        ScalarVisibleFaces.visible_faces(request, options)
    }
}

#[test]
fn register_conformant_refuses_an_unverified_perspective_claim() {
    use axiolid_mesh_hidden_lines_contract::conformance::{
        faces::FaceConformanceFailure, ConformanceFailure,
    };

    // A hidden-line provider that ignores the projection fails the fixtures.
    let mut lines = MeshHiddenLinesRegistry::new();
    let report = lines
        .register_conformant(0, Ignores)
        .expect_err("the claim is false");
    assert!(report.failures.iter().any(|f| matches!(
        f,
        ConformanceFailure::ReturnedError {
            case: "perspective claim",
            ..
        }
    )));
    assert_eq!(lines.providers().count(), 0);

    let mut faces = MeshVisibleFacesRegistry::new();
    let report = faces
        .register_conformant(0, Claims)
        .expect_err("the claim cannot be verified");
    assert!(report.failures.iter().any(|f| matches!(
        f,
        FaceConformanceFailure::ReturnedError {
            case: "perspective claim",
            ..
        }
    )));
    assert_eq!(faces.providers().count(), 0);

    // Without the claim the same provider registers as before.
    lines
        .register_conformant(0, ScalarHiddenLines::new())
        .expect("conformant");
    faces
        .register_conformant(0, ScalarVisibleFaces::new())
        .expect("conformant");
}
