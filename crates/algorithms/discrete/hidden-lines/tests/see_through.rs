//! See-through triangles (ADR 0111). A provider that does not honour them
//! keeps `supports_occlusion` at its default and refuses them by name when
//! called directly.

mod common;

use axiolid_contracts::{ExecutionOptions, GeomError, Operation};
use axiolid_core::Tolerance;
use axiolid_hidden_lines::{
    ScalarHiddenLines, ScalarVisibleFaces, ScanHiddenLines, ScanVisibleFaces, UnprunedVisibleFaces,
};
use axiolid_mesh_hidden_lines_contract::{
    EdgeSelection, HiddenLineMesh, HiddenLineRequest, MeshHiddenLines, MeshId, MeshRole,
    MeshVisibleFaces, Occlusion, View, VisibleFaceRequest,
};
use common::{box_mesh, plan_frame, FACE_LIMITS, LIMITS};

fn options() -> ExecutionOptions {
    ExecutionOptions::new(Tolerance::new(1e-9, 1e-9).expect("tolerance"))
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
