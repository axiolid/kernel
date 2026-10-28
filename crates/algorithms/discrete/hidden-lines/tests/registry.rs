//! The provider behind `MeshHiddenLinesRegistry`.

mod common;

use axiolid_contracts::{ExecutionOptions, GeomError, Operation};
use axiolid_core::{Point3, Tolerance};
use axiolid_dispatch::MeshHiddenLinesRegistry;
use axiolid_hidden_lines::ScalarHiddenLines;
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
fn a_near_plane_through_a_solid_is_refused_and_names_the_mesh() {
    let cube = box_mesh([0.0; 3], [1.0; 3]);
    let sheet = box_mesh([5.0; 3], [6.0; 3]);
    let meshes = [
        HiddenLineMesh::new(MeshId(4), MeshRole::Sheet, &sheet),
        HiddenLineMesh::new(MeshId(9), MeshRole::Solid, &cube),
    ];
    let request = HiddenLineRequest::new(
        &meshes,
        OrthographicView::new(plan_frame()).with_near(0.5),
        EdgeSelection::new(Some(0.5)),
        LIMITS,
    );
    let error = registry()
        .hidden_lines(&request, &options())
        .expect_err("refused");
    assert_eq!(error.mesh, Some(MeshId(9)));
    assert!(matches!(
        error.error,
        GeomError::UnsupportedInput {
            operation: Operation::MeshHiddenLines,
            ..
        }
    ));
}
