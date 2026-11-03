#![cfg(all(
    feature = "dispatch-mesh-hidden-lines",
    feature = "hidden-lines-provider"
))]

use axiolid::contracts::ExecutionOptions;
use axiolid::core::{Frame3, Tolerance};
use axiolid::dispatch::MeshVisibleFacesRegistry;
use axiolid::mesh::TriMesh;
use axiolid::mesh_hidden_lines::{
    HiddenLineMesh, MeshId, MeshRole, OrthographicView, ScalarVisibleFaces, VisibleFaceLimits,
    VisibleFaceRequest,
};
use axiolid::{Point3, Vec3};

/// The facade re-export is the exact provider, and it answers through the
/// dispatch registry like any other.
#[test]
fn the_scalar_visible_faces_provider_answers_through_the_registry() {
    // A unit square sheet seen from above along +z.
    let sheet = TriMesh::new(
        vec![
            Point3::new(0.0, 0.0, 1.0),
            Point3::new(1.0, 0.0, 1.0),
            Point3::new(1.0, 1.0, 1.0),
            Point3::new(0.0, 1.0, 1.0),
        ],
        vec![0, 1, 2, 0, 2, 3],
    );
    let meshes = [HiddenLineMesh::new(MeshId(1), MeshRole::Sheet, &sheet)];
    let view = OrthographicView::new(Frame3 {
        origin: Point3::ZERO,
        x: Vec3::X,
        y: Vec3::Y,
        z: Vec3::Z,
    });
    let request = VisibleFaceRequest::new(
        &meshes,
        view,
        VisibleFaceLimits::new(4, 64, 64, 1024, 16, 256),
    );
    let mut registry = MeshVisibleFacesRegistry::new();
    registry
        .register_conformant(0, ScalarVisibleFaces::new())
        .expect("the exact provider passes the conformance suite");
    let options = ExecutionOptions::new(Tolerance::new(1.0e-9, 1.0e-9).unwrap());
    let outcome = registry
        .visible_faces(&request, &options)
        .expect("the sheet resolves");
    assert_eq!(outcome.regions.len(), 1);
    assert_eq!(outcome.regions[0].mesh, MeshId(1));
}
