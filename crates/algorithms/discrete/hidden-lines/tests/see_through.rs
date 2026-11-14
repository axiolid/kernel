//! See-through triangles (ADR 0111). A provider that does not honour them
//! keeps `supports_occlusion` at its default and refuses them by name when
//! called directly.

mod common;

use axiolid_contracts::{Backend, BackendDescriptor, ExecutionOptions, GeomError, Operation};
use axiolid_core::Tolerance;
use axiolid_dispatch::{MeshHiddenLinesRegistry, MeshVisibleFacesRegistry};
use axiolid_hidden_lines::{
    ScalarHiddenLines, ScalarVisibleFaces, ScanHiddenLines, ScanVisibleFaces, UnprunedVisibleFaces,
};
use axiolid_mesh_hidden_lines_contract::conformance::{
    faces::FaceConformanceFailure, ConformanceFailure,
};
use axiolid_mesh_hidden_lines_contract::{
    EdgeSelection, HiddenLineError, HiddenLineMesh, HiddenLineOutcome, HiddenLineRequest,
    MeshHiddenLines, MeshId, MeshRole, MeshVisibleFaces, Occlusion, TriangleOcclusion, View,
    VisibleFaceOutcome, VisibleFaceRequest,
};
use common::{box_mesh, plan_frame, FACE_LIMITS, LIMITS};

fn options() -> ExecutionOptions {
    ExecutionOptions::new(Tolerance::new(1e-9, 1e-9).expect("tolerance"))
}

/// A provider that claims see-through triangles and lets them hide: it draws
/// the request with every flag dropped.
#[derive(Debug)]
struct Ignores;

impl Backend for Ignores {
    fn descriptor(&self) -> BackendDescriptor {
        ScalarHiddenLines.descriptor()
    }
}

impl MeshHiddenLines for Ignores {
    fn supports_occlusion(&self, _occlusion: Occlusion) -> bool {
        true
    }

    fn hidden_lines(
        &self,
        request: &HiddenLineRequest<'_>,
        options: &ExecutionOptions,
    ) -> Result<HiddenLineOutcome, HiddenLineError> {
        let meshes: Vec<HiddenLineMesh<'_>> = request
            .meshes
            .iter()
            .map(|entry| {
                let mut entry = *entry;
                entry.occlusion = TriangleOcclusion::AllOpaque;
                entry
            })
            .collect();
        let mut opaque = request.clone();
        opaque.meshes = &meshes;
        ScalarHiddenLines.hidden_lines(&opaque, options)
    }
}

impl MeshVisibleFaces for Ignores {
    fn supports_occlusion(&self, _occlusion: Occlusion) -> bool {
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
fn register_conformant_refuses_an_unverified_see_through_claim() {
    let mut lines = MeshHiddenLinesRegistry::new();
    let report = lines
        .register_conformant(0, Ignores)
        .expect_err("the claim is false");
    assert!(report.failures.iter().any(|f| matches!(
        f,
        ConformanceFailure::ReturnedError {
            case: "see-through claim",
            ..
        }
    )));
    assert_eq!(lines.providers().count(), 0);

    let mut faces = MeshVisibleFacesRegistry::new();
    let report = faces
        .register_conformant(0, Ignores)
        .expect_err("the claim is not verified");
    assert!(report.failures.iter().any(|f| matches!(
        f,
        FaceConformanceFailure::ReturnedError {
            case: "see-through claim",
            ..
        }
    )));
    assert_eq!(faces.providers().count(), 0);
}

fn is_named_refusal(error: &GeomError, operation: Operation) -> bool {
    matches!(
        error,
        GeomError::UnsupportedInput { operation: o, input: "see-through triangles", .. } if *o == operation
    )
}

#[test]
fn providers_that_do_not_claim_see_through_refuse_it_by_name() {
    let mesh = box_mesh([0.0, 0.0, 0.0], [1.0, 1.0, 1.0]);
    let mut flags = [Occlusion::Opaque; 12];
    flags[3] = Occlusion::SeeThrough;
    let meshes = [HiddenLineMesh::new(MeshId(1), MeshRole::Solid, &mesh).with_occlusion(&flags)];
    let view = View::new(plan_frame());
    let lines: [&dyn MeshHiddenLines; 2] = [&ScalarHiddenLines, &ScanHiddenLines];
    for provider in lines {
        assert!(!provider.supports_occlusion(Occlusion::SeeThrough));
        let request = HiddenLineRequest::new(&meshes, view, EdgeSelection::new(None), LIMITS);
        let error = provider
            .hidden_lines(&request, &options())
            .expect_err("refused");
        assert!(
            is_named_refusal(&error.error, Operation::MeshHiddenLines),
            "{error:?}"
        );
    }
    let faces: [&dyn MeshVisibleFaces; 3] = [
        &ScalarVisibleFaces,
        &ScanVisibleFaces,
        &UnprunedVisibleFaces,
    ];
    for provider in faces {
        assert!(!provider.supports_occlusion(Occlusion::SeeThrough));
        let request = VisibleFaceRequest::new(&meshes, view, FACE_LIMITS);
        let error = provider
            .visible_faces(&request, &options())
            .expect_err("refused");
        assert!(
            is_named_refusal(&error.error, Operation::MeshVisibleFaces),
            "{error:?}"
        );
    }
}
