//! The provider behind `MeshHiddenLinesRegistry`.

mod common;

use axiolid_contracts::ExecutionOptions;
use axiolid_core::{Point3, Tolerance};
use axiolid_dispatch::MeshHiddenLinesRegistry;
use axiolid_hidden_lines::ScalarHiddenLines;
use axiolid_mesh::TriMesh;
use axiolid_mesh_hidden_lines_contract::{
    EdgeKind, EdgeSelection, HiddenLineMesh, HiddenLineRequest, MeshId, MeshRole, OrthographicView,
    Visibility,
};
use common::{box_mesh, isometric_frame, plan_frame, LIMITS};

fn options() -> ExecutionOptions {
    ExecutionOptions::new(Tolerance::new(1e-9, 1e-9).expect("tolerance"))
}

fn registry() -> MeshHiddenLinesRegistry {
    let mut registry = MeshHiddenLinesRegistry::new();
    registry
        .register_conformant(0, ScalarHiddenLines::new())
        .unwrap_or_else(|report| panic!("{report}"));
    registry
}

#[test]
fn the_provider_registers_as_conformant() {
    assert_eq!(registry().providers().count(), 1);
}

#[test]
fn a_request_runs_end_to_end_through_the_registry() {
    // Box behind box in plan view (the conformance scene): visible length 14,
    // hidden length 18.
    let front = box_mesh([0.0, 0.0, 2.0], [2.0, 2.0, 3.0]);
    let back = box_mesh([1.0, 1.0, 0.0], [3.0, 3.0, 1.0]);
    let meshes = [
        HiddenLineMesh::new(MeshId(1), MeshRole::Solid, &front),
        HiddenLineMesh::new(MeshId(2), MeshRole::Solid, &back),
    ];
    let request = HiddenLineRequest::new(
        &meshes,
        OrthographicView::new(plan_frame()),
        EdgeSelection::new(Some(0.5)),
        LIMITS,
    );
    let outcome = registry()
        .hidden_lines(&request, &options())
        .expect("drawn");
    assert_eq!(outcome.view, request.view);
    assert_eq!(outcome.evidence.selected_edges, 24);
    assert_eq!(outcome.evidence.point_edges, 8);
    let total = |visibility| common::length(&outcome, visibility);
    assert!((total(Visibility::Visible) - 14.0).abs() < 1e-12);
    assert!((total(Visibility::Hidden) - 18.0).abs() < 1e-12);
}

#[test]
fn the_isometric_cube_through_the_registry_has_the_analytic_kinds() {
    let cube = box_mesh([0.0; 3], [1.0; 3]);
    let meshes = [HiddenLineMesh::new(MeshId(1), MeshRole::Solid, &cube)];
    let request = HiddenLineRequest::new(
        &meshes,
        OrthographicView::new(isometric_frame(Point3::new(0.1, 0.2, 0.3))),
        EdgeSelection::new(Some(0.5)),
        LIMITS,
    );
    let outcome = registry()
        .hidden_lines(&request, &options())
        .expect("drawn");
    let silhouettes = outcome
        .segments
        .iter()
        .filter(|s| s.kind == EdgeKind::Silhouette)
        .count();
    assert_eq!(silhouettes, 6);
    assert_eq!(outcome.evidence.visible_segments, 9);
    assert_eq!(outcome.evidence.hidden_segments, 3);
}

#[test]
fn a_near_plane_through_a_solid_caps_it() {
    // The cube is cut at 0.5: its top ring is removed and the cap hides the
    // bottom ring.
    let cube = box_mesh([0.0; 3], [1.0; 3]);
    let sheet = box_mesh([5.0; 3], [6.0; 3]);
    let meshes = [
        HiddenLineMesh::new(MeshId(4), MeshRole::Sheet, &sheet),
        HiddenLineMesh::new(MeshId(9), MeshRole::Solid, &cube),
    ];
    let outcome = registry()
        .hidden_lines(&plan_request(&meshes, 0.5), &options())
        .expect("drawn");
    assert_eq!(outcome.evidence.visible_segments, 0);
    assert_eq!(outcome.evidence.hidden_segments, 4);
    assert!(outcome.segments.iter().all(|s| s.source.mesh == MeshId(9)));
}

#[test]
fn a_face_in_the_near_plane_is_drawn() {
    // The base of the tetrahedron lies in the plane W = 1 and its apex is
    // beyond it: the cap and the base coincide, and a floor quad under the
    // base is hidden by it.
    let tetra = TriMesh::new(
        vec![
            Point3::new(0.0, 0.0, 1.0),
            Point3::new(1.0, 0.0, 1.0),
            Point3::new(0.0, 1.0, 1.0),
            Point3::new(0.0, 0.0, 3.0),
        ],
        vec![0, 2, 1, 0, 1, 3, 1, 2, 3, 2, 0, 3],
    );
    let floor = box_mesh([0.1, 0.1, -1.0], [0.4, 0.4, 0.0]);
    let meshes = [
        HiddenLineMesh::new(MeshId(3), MeshRole::Solid, &floor),
        HiddenLineMesh::new(MeshId(7), MeshRole::Solid, &tetra),
    ];
    let outcome = registry()
        .hidden_lines(&plan_request(&meshes, 1.0), &options())
        .expect("drawn");
    assert!(outcome
        .segments
        .iter()
        .filter(|s| s.source.mesh == MeshId(3))
        .all(|s| s.visibility == Visibility::Hidden));
}

#[test]
fn a_beam_whose_underside_is_at_the_cut_height_is_drawn() {
    // The beam [0,2] x [0,0.3] x [1,1.5] is cut at 1 on its underside; the
    // floor slab below it has its top edge y = 0 hidden under the beam over
    // x in [0, 2] and visible beyond.
    let beam = box_mesh([0.0, 0.0, 1.0], [2.0, 0.3, 1.5]);
    let floor = box_mesh([0.0, 0.0, 0.0], [3.0, 1.0, 0.5]);
    let meshes = [
        HiddenLineMesh::new(MeshId(1), MeshRole::Solid, &floor),
        HiddenLineMesh::new(MeshId(2), MeshRole::Solid, &beam),
    ];
    let outcome = registry()
        .hidden_lines(&plan_request(&meshes, 1.0), &options())
        .expect("drawn");
    let along = |visibility| {
        outcome
            .segments
            .iter()
            .filter(|s| s.source.mesh == MeshId(1) && s.source.vertices == [4, 5])
            .filter(|s| s.visibility == visibility)
            .map(|s| (s.end - s.start).length())
            .sum::<f64>()
    };
    assert!((along(Visibility::Hidden) - 2.0).abs() < 1e-12);
    assert!((along(Visibility::Visible) - 1.0).abs() < 1e-12);
}

fn plan_request<'a>(meshes: &'a [HiddenLineMesh<'a>], near: f64) -> HiddenLineRequest<'a> {
    HiddenLineRequest::new(
        meshes,
        OrthographicView::new(plan_frame()).with_near(near),
        EdgeSelection::new(Some(0.5)),
        LIMITS,
    )
}

#[test]
fn a_near_plane_beyond_a_solid_draws_it() {
    let cube = box_mesh([0.0; 3], [1.0; 3]);
    let meshes = [HiddenLineMesh::new(MeshId(1), MeshRole::Solid, &cube)];
    let outcome = registry()
        .hidden_lines(&plan_request(&meshes, 5.0), &options())
        .expect("drawn");
    assert_eq!(outcome.evidence.clipped_edges, 0);
    assert_eq!(outcome.evidence.visible_segments, 4);
    assert_eq!(outcome.evidence.hidden_segments, 4);
}

#[test]
fn only_the_solid_that_the_near_plane_cuts_is_capped() {
    let low = box_mesh([0.0; 3], [1.0; 3]);
    let tall = box_mesh([2.0, 0.0, 0.0], [3.0, 1.0, 3.0]);
    let meshes = [
        HiddenLineMesh::new(MeshId(1), MeshRole::Solid, &low),
        HiddenLineMesh::new(MeshId(2), MeshRole::Solid, &tall),
    ];
    let outcome = registry()
        .hidden_lines(&plan_request(&meshes, 2.0), &options())
        .expect("drawn");
    // The low cube is drawn whole (4 visible top, 4 hidden bottom); the tall
    // one only shows its bottom ring, hidden by its cap.
    assert_eq!(outcome.evidence.visible_segments, 4);
    assert_eq!(outcome.evidence.hidden_segments, 8);
    assert_eq!(outcome.evidence.clipped_edges, 4);
}

#[test]
fn a_solid_wholly_beyond_the_near_plane_is_removed() {
    let low = box_mesh([0.0; 3], [1.0; 3]);
    let beyond = box_mesh([0.0, 0.0, 3.0], [1.0, 1.0, 4.0]);
    let alone = [HiddenLineMesh::new(MeshId(1), MeshRole::Solid, &low)];
    let both = [
        HiddenLineMesh::new(MeshId(1), MeshRole::Solid, &low),
        HiddenLineMesh::new(MeshId(2), MeshRole::Solid, &beyond),
    ];
    let expected = registry()
        .hidden_lines(&plan_request(&alone, 2.0), &options())
        .expect("drawn");
    let outcome = registry()
        .hidden_lines(&plan_request(&both, 2.0), &options())
        .expect("drawn");
    // It does not occlude the cube below it and draws nothing itself.
    assert_eq!(outcome.segments, expected.segments);
    assert_eq!(outcome.evidence.clipped_edges, 8);
}
