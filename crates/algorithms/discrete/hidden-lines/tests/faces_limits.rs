//! The input limits and the memory budget of `ScalarVisibleFaces`, asked of
//! the provider directly.

mod common;

use axiolid_contracts::{ExecutionOptions, GeomError};
use axiolid_core::{Point3, Tolerance};
use axiolid_hidden_lines::ScalarVisibleFaces;
use axiolid_mesh::TriMesh;
use axiolid_mesh_hidden_lines_contract::{
    HiddenLineMesh, MeshId, MeshRole, MeshVisibleFaces, OrthographicView, VisibleFaceLimits,
    VisibleFaceRequest,
};
use common::{box_mesh, plan_frame, solids};

fn options() -> ExecutionOptions {
    ExecutionOptions::new(Tolerance::new(1e-9, 1e-9).expect("tolerance"))
}

fn refusal(limits: VisibleFaceLimits) -> Option<GeomError> {
    // Three boxes: 3 meshes, 24 vertices, 36 triangles.
    let boxes = [
        box_mesh([0.0, 0.0, 0.0], [1.0, 1.0, 1.0]),
        box_mesh([2.0, 0.0, 0.0], [3.0, 1.0, 1.0]),
        box_mesh([4.0, 0.0, 0.0], [5.0, 1.0, 1.0]),
    ];
    let meshes = solids(&boxes);
    let request = VisibleFaceRequest::new(&meshes, OrthographicView::new(plan_frame()), limits);
    ScalarVisibleFaces::new()
        .visible_faces(&request, &options())
        .err()
        .map(|e| e.error)
}

#[test]
fn the_input_limits_are_enforced_and_named() {
    let limits = |meshes, vertices, triangles| {
        VisibleFaceLimits::new(meshes, vertices, triangles, 1_000_000, 10_000, 100_000)
    };
    assert_eq!(refusal(limits(3, 24, 36)), None);
    assert_eq!(
        refusal(limits(2, 24, 36)),
        Some(GeomError::BudgetExceeded {
            resource: "visible-face meshes"
        })
    );
    assert_eq!(
        refusal(limits(3, 23, 36)),
        Some(GeomError::BudgetExceeded {
            resource: "visible-face source vertices"
        })
    );
    assert_eq!(
        refusal(limits(3, 24, 35)),
        Some(GeomError::BudgetExceeded {
            resource: "visible-face source triangles"
        })
    );
}

/// `n` sheets through the line y = 5, z = 0 at different slopes and
/// extents, over a flat sheet: every pair meets along that one line.
fn along_one_line(n: usize) -> Vec<TriMesh> {
    let mut sheets = vec![TriMesh::new(
        vec![
            Point3::new(-1.0, 0.0, 0.0),
            Point3::new(20.0, 0.0, 0.0),
            Point3::new(20.0, 10.0, 0.0),
            Point3::new(-1.0, 10.0, 0.0),
        ],
        vec![0, 1, 2, 0, 2, 3],
    )];
    for k in 0..n {
        let (x0, x1) = (k as f64 * 0.5, 10.0 + k as f64 * 0.25);
        let slope = 1.0 + k as f64;
        sheets.push(TriMesh::new(
            vec![
                Point3::new(x0, 4.0, -slope),
                Point3::new(x1, 4.0, -slope),
                Point3::new(x1, 6.0, slope),
                Point3::new(x0, 6.0, slope),
            ],
            vec![0, 1, 2, 0, 2, 3],
        ));
    }
    sheets
}

#[test]
fn intersection_segments_count_against_the_memory_budget() {
    let sheets = along_one_line(24);
    let meshes: Vec<HiddenLineMesh<'_>> = sheets
        .iter()
        .enumerate()
        .map(|(i, m)| HiddenLineMesh::new(MeshId(i as u64 + 1), MeshRole::Sheet, m))
        .collect();
    let limits = VisibleFaceLimits::new(100, 10_000, 10_000, 100_000_000, 100_000, 1_000_000);
    let request = VisibleFaceRequest::new(&meshes, OrthographicView::new(plan_frame()), limits);
    let free = ScalarVisibleFaces::new()
        .visible_faces(&request, &options())
        .expect("answers without a budget");
    let (triangles, segments) = (
        free.evidence.source_triangles,
        free.evidence.intersection_segments,
    );
    assert!(segments >= 24, "{segments} intersection segments");
    // 3,072 bytes per triangle are declared; 1,024 per intersection segment.
    let base = 3072 * triangles;
    let roomy = options().with_memory_budget(base + 1024 * segments + 256 * segments * segments);
    let answered = ScalarVisibleFaces::new().visible_faces(&request, &roomy);
    assert!(answered.is_ok(), "{:?}", answered.err());
    let tight = options().with_memory_budget(base + 1024 * (segments / 2));
    let error = ScalarVisibleFaces::new()
        .visible_faces(&request, &tight)
        .expect_err("over the memory budget");
    assert_eq!(
        error.error,
        GeomError::BudgetExceeded { resource: "memory" }
    );
}

/// Intersection segments that never meet each other cost no crossing and no
/// split, so only their own charge can refuse them.
#[test]
fn separate_intersection_segments_count_against_the_memory_budget() {
    let mut sheets = vec![TriMesh::new(
        vec![
            Point3::new(-1.0, 0.0, 0.0),
            Point3::new(60.0, 0.0, 0.0),
            Point3::new(60.0, 10.0, 0.0),
            Point3::new(-1.0, 10.0, 0.0),
        ],
        vec![0, 1, 2, 0, 2, 3],
    )];
    for k in 0..20 {
        let x = 3.0 * k as f64;
        sheets.push(TriMesh::new(
            vec![
                Point3::new(x, 4.0, -1.0),
                Point3::new(x + 1.0, 4.0, -1.0),
                Point3::new(x + 1.0, 6.0, 1.0),
                Point3::new(x, 6.0, 1.0),
            ],
            vec![0, 1, 2, 0, 2, 3],
        ));
    }
    let meshes: Vec<HiddenLineMesh<'_>> = sheets
        .iter()
        .enumerate()
        .map(|(i, m)| HiddenLineMesh::new(MeshId(i as u64 + 1), MeshRole::Sheet, m))
        .collect();
    let limits = VisibleFaceLimits::new(100, 10_000, 10_000, 100_000_000, 100_000, 1_000_000);
    let request = VisibleFaceRequest::new(&meshes, OrthographicView::new(plan_frame()), limits);
    let free = ScalarVisibleFaces::new()
        .visible_faces(&request, &options())
        .expect("answers without a budget");
    assert_eq!(free.evidence.crossings, 0);
    let (triangles, segments) = (
        free.evidence.source_triangles,
        free.evidence.intersection_segments,
    );
    assert!(segments >= 20, "{segments} intersection segments");
    let base = 3072 * triangles;
    assert!(ScalarVisibleFaces::new()
        .visible_faces(
            &request,
            &options().with_memory_budget(base + 1024 * segments + 256 * segments * segments)
        )
        .is_ok());
    let tight = options().with_memory_budget(base + 1024 * (segments / 2));
    let error = ScalarVisibleFaces::new()
        .visible_faces(&request, &tight)
        .expect_err("over the memory budget");
    assert_eq!(
        error.error,
        GeomError::BudgetExceeded { resource: "memory" }
    );
}
