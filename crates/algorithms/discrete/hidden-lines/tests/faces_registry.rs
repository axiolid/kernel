//! The provider behind `MeshVisibleFacesRegistry`.

mod common;

use axiolid_contracts::ExecutionOptions;
use axiolid_core::Tolerance;
use axiolid_dispatch::MeshVisibleFacesRegistry;
use axiolid_hidden_lines::ScalarVisibleFaces;
use axiolid_mesh_hidden_lines_contract::{FaceGroup, FaceSource, View, VisibleFaceRequest};
use common::{box_mesh, plan_frame, solids, source_area, FACE_LIMITS};

fn registry() -> MeshVisibleFacesRegistry {
    let mut registry = MeshVisibleFacesRegistry::new();
    registry
        .register_conformant(0, ScalarVisibleFaces::new())
        .unwrap_or_else(|report| panic!("{report}"));
    registry
}

#[test]
fn a_request_runs_end_to_end_through_the_registry() {
    // Box behind box in plan view (the hidden-line registry scene).
    let boxes = [
        box_mesh([0.0, 0.0, 2.0], [2.0, 2.0, 3.0]),
        box_mesh([1.0, 1.0, 0.0], [3.0, 3.0, 1.0]),
    ];
    let meshes = solids(&boxes);
    let request = VisibleFaceRequest::new(&meshes, View::new(plan_frame()), FACE_LIMITS);
    let options = ExecutionOptions::new(Tolerance::new(1e-9, 1e-9).expect("tolerance"));
    let outcome = registry().visible_faces(&request, &options).expect("faces");
    let top = FaceSource::Group(FaceGroup(0));
    assert_eq!(source_area(&outcome, 1, top), 4.0);
    assert_eq!(source_area(&outcome, 2, top), 4.0 - 1.0);
}
