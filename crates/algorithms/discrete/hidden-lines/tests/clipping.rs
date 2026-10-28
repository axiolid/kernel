//! A triangle cut by the near plane hides only the part that is kept.

mod common;

use axiolid_core::Point3;
use axiolid_mesh::TriMesh;
use axiolid_mesh_hidden_lines_contract::{
    EdgeSelection, HiddenLineMesh, HiddenLineRequest, MeshId, MeshRole, OrthographicView,
    SourceEdge, Visibility,
};
use common::{plan_frame, scalar, LIMITS};

fn triangle(corners: [[f64; 3]; 3]) -> TriMesh {
    TriMesh::new(
        corners
            .iter()
            .map(|c| Point3::new(c[0], c[1], c[2]))
            .collect(),
        vec![0, 1, 2],
    )
}

/// The pieces (visibility, start x, end x) of the edge from vertex 0 to 1 of
/// the lower triangle, drawn under a sloping one, with the given near plane.
fn pieces(near: Option<f64>) -> Vec<(Visibility, f64, f64)> {
    // The slope rises from W = 0 at x = 0 to W = 4 at x = 8, over y = 1.
    let slope = triangle([[0.0, 0.0, 0.0], [8.0, 0.0, 4.0], [0.0, 8.0, 0.0]]);
    // The edge from (1, 1) to (3, 1) lies at W = -1, below it.
    let floor = triangle([[1.0, 1.0, -1.0], [3.0, 1.0, -1.0], [1.0, -2.0, -1.0]]);
    let meshes = [
        HiddenLineMesh::new(MeshId(1), MeshRole::Sheet, &slope),
        HiddenLineMesh::new(MeshId(2), MeshRole::Sheet, &floor),
    ];
    let mut view = OrthographicView::new(plan_frame());
    if let Some(near) = near {
        view = view.with_near(near);
    }
    let request = HiddenLineRequest::new(&meshes, view, EdgeSelection::new(None), LIMITS);
    let options = axiolid_contracts::ExecutionOptions::new(
        axiolid_core::Tolerance::new(1e-9, 1e-9).expect("tolerance"),
    );
    let outcome = axiolid_mesh_hidden_lines_contract::MeshHiddenLines::hidden_lines(
        &scalar(),
        &request,
        &options,
    )
    .expect("drawn");
    let edge = SourceEdge {
        mesh: MeshId(2),
        vertices: [0, 1],
    };
    outcome
        .segments
        .iter()
        .filter(|s| s.source == edge)
        .map(|s| (s.visibility, s.start.x, s.end.x))
        .collect()
}

#[test]
fn without_a_near_plane_the_slope_hides_the_whole_edge() {
    assert_eq!(pieces(None), [(Visibility::Hidden, 1.0, 3.0)]);
}

#[test]
fn a_near_plane_through_the_slope_hides_only_the_kept_part() {
    // The slope is at W = x / 2, so with near = 1 it is kept up to x = 2.
    assert_eq!(
        pieces(Some(1.0)),
        [
            (Visibility::Hidden, 1.0, 2.0),
            (Visibility::Visible, 2.0, 3.0)
        ]
    );
}

#[test]
fn a_near_plane_above_the_slope_changes_nothing() {
    assert_eq!(pieces(Some(10.0)), [(Visibility::Hidden, 1.0, 3.0)]);
}
