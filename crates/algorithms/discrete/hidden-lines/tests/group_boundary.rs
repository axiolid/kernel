//! `EdgeKind::GroupBoundary` (ADR 0106): an edge between triangles of
//! different face groups that is neither a silhouette nor a selected crease.

mod common;

use axiolid_contracts::GeomError;
use axiolid_core::{Frame3, Point2, Point3, Vec3};
use axiolid_hidden_lines::{ScalarHiddenLines, ScalarVisibleFaces};
use axiolid_mesh::TriMesh;
use axiolid_mesh_hidden_lines_contract::{
    EdgeKind, EdgeSelection, FaceGroup, HiddenLineMesh, HiddenLineSegment, MeshId, MeshRole,
    OrthographicView, Visibility,
};
use common::{box_mesh, faces, frame, plan_frame};

fn square() -> TriMesh {
    TriMesh::new(
        vec![
            Point3::new(0.0, 0.0, 0.0),
            Point3::new(1.0, 0.0, 0.0),
            Point3::new(1.0, 1.0, 0.0),
            Point3::new(0.0, 1.0, 0.0),
        ],
        vec![0, 1, 2, 0, 2, 3],
    )
}

fn plan() -> OrthographicView {
    OrthographicView::new(plan_frame())
}

/// `(vertices, kind, visibility, start, end)` of each segment.
type Row = ([u32; 2], EdgeKind, Visibility, [f64; 2], [f64; 2]);

fn rows(segments: &[HiddenLineSegment]) -> Vec<Row> {
    segments
        .iter()
        .map(|s| {
            (
                s.source.vertices,
                s.kind,
                s.visibility,
                [s.start.x, s.start.y],
                [s.end.x, s.end.y],
            )
        })
        .collect()
}

#[test]
fn a_square_in_two_groups_draws_its_diagonal() {
    use EdgeKind::{Boundary, GroupBoundary};
    use Visibility::Visible;
    let mesh = square();
    let groups = [FaceGroup(0), FaceGroup(1)];
    let meshes = [HiddenLineMesh::new(MeshId(1), MeshRole::Sheet, &mesh).with_groups(&groups)];
    let selection = EdgeSelection::new(None).with_group_boundary(true);
    let outcome = draw_plan(&meshes, selection).expect("drawing");
    assert_eq!(
        rows(&outcome.segments),
        [
            ([0, 1], Boundary, Visible, [0.0, 0.0], [1.0, 0.0]),
            ([0, 2], GroupBoundary, Visible, [0.0, 0.0], [1.0, 1.0]),
            ([0, 3], Boundary, Visible, [0.0, 0.0], [0.0, 1.0]),
            ([1, 2], Boundary, Visible, [1.0, 0.0], [1.0, 1.0]),
            ([2, 3], Boundary, Visible, [1.0, 1.0], [0.0, 1.0]),
        ]
    );
    assert_eq!(outcome.evidence.selected_edges, 5);
}

#[test]
fn without_the_selection_the_groups_are_ignored_bit_for_bit() {
    let mesh = square();
    let groups = [FaceGroup(0), FaceGroup(1)];
    let grouped = [HiddenLineMesh::new(MeshId(1), MeshRole::Sheet, &mesh).with_groups(&groups)];
    let plain = [HiddenLineMesh::new(MeshId(1), MeshRole::Sheet, &mesh)];
    let selection = EdgeSelection::new(None);
    let (a, b) = (
        draw_plan(&grouped, selection).expect("grouped"),
        draw_plan(&plain, selection).expect("plain"),
    );
    assert_eq!(rows(&a.segments).len(), 4);
    assert!(a.segments.iter().all(|s| s.kind == EdgeKind::Boundary));
    assert_eq!(a, b);
    // Even groups of the wrong length are not read.
    let short = [FaceGroup(0)];
    let wrong = [HiddenLineMesh::new(MeshId(1), MeshRole::Sheet, &mesh).with_groups(&short)];
    assert_eq!(draw_plan(&wrong, selection).expect("wrong"), a);
}

#[test]
fn one_group_has_no_group_boundary() {
    let mesh = square();
    let single = [HiddenLineMesh::new(MeshId(1), MeshRole::Sheet, &mesh)];
    let selection = EdgeSelection::new(None).with_group_boundary(true);
    let outcome = draw_plan(&single, selection).expect("drawing");
    assert_eq!(outcome.segments.len(), 4);
    assert!(outcome
        .segments
        .iter()
        .all(|s| s.kind == EdgeKind::Boundary));
    // A planar square is one group as well.
    let planar = [HiddenLineMesh::new(MeshId(1), MeshRole::Sheet, &mesh).with_planar_groups()];
    assert_eq!(
        draw_plan(&planar, selection).expect("planar").segments,
        outcome.segments
    );
}

#[test]
fn a_group_slice_of_the_wrong_length_names_the_mesh_when_selected() {
    let mesh = square();
    let short = [FaceGroup(0)];
    let meshes = [
        HiddenLineMesh::new(MeshId(1), MeshRole::Sheet, &mesh),
        HiddenLineMesh::new(MeshId(2), MeshRole::Sheet, &mesh).with_groups(&short),
    ];
    let selection = EdgeSelection::new(None).with_group_boundary(true);
    let error = draw_plan(&meshes, selection).expect_err("length");
    assert_eq!(error.mesh, Some(MeshId(2)));
    assert!(matches!(error.error, GeomError::InvalidInput(_)));
}

#[test]
fn the_drawn_diagonal_is_the_edge_the_face_regions_share() {
    let mesh = square();
    let groups = [FaceGroup(0), FaceGroup(1)];
    let meshes = [HiddenLineMesh::new(MeshId(1), MeshRole::Sheet, &mesh).with_groups(&groups)];
    let outcome = faces(&ScalarVisibleFaces::new(), &meshes, plan());
    assert_eq!(outcome.regions.len(), 2);
    let edges = |ring: &[Point2]| -> Vec<[(u64, u64); 2]> {
        let key = |p: Point2| (p.x.to_bits(), p.y.to_bits());
        (0..ring.len())
            .map(|i| {
                let (a, b) = (key(ring[i]), key(ring[(i + 1) % ring.len()]));
                if a <= b {
                    [a, b]
                } else {
                    [b, a]
                }
            })
            .collect()
    };
    let (first, second) = (
        edges(&outcome.regions[0].outer),
        edges(&outcome.regions[1].outer),
    );
    let shared: Vec<_> = first.iter().filter(|e| second.contains(e)).collect();
    let key = |x: f64, y: f64| (x.to_bits(), y.to_bits());
    assert_eq!(shared, [&[key(0.0, 0.0), key(1.0, 1.0)]]);
}

/// The frame of `tests/c/smoke.c`: looking at the cube's corner `(1,1,1)`.
fn isometric() -> Frame3 {
    let (s2, s6, s3) = (
        std::f64::consts::FRAC_1_SQRT_2,
        1.0 / 6f64.sqrt(),
        1.0 / 3f64.sqrt(),
    );
    frame(
        Vec3::new(-s2, s2, 0.0),
        Vec3::new(-s6, -s6, 2.0 * s6),
        Vec3::new(s3, s3, s3),
        Point3::ZERO,
    )
}

fn cube_kinds(selection: EdgeSelection) -> Vec<([u32; 2], EdgeKind, Visibility)> {
    let cube = box_mesh([0.0; 3], [2.0; 3]);
    let meshes = [HiddenLineMesh::new(MeshId(1), MeshRole::Solid, &cube).with_planar_groups()];
    let outcome = common::try_draw_with(
        &ScalarHiddenLines,
        &meshes,
        OrthographicView::new(isometric()),
        selection,
    )
    .expect("drawing");
    outcome
        .segments
        .iter()
        .map(|s| (s.source.vertices, s.kind, s.visibility))
        .collect()
}

#[test]
fn the_planar_isometric_cube_has_group_boundaries_at_the_corners() {
    use EdgeKind::{GroupBoundary, Silhouette};
    use Visibility::{Hidden, Visible};
    let found = cube_kinds(EdgeSelection::new(None).with_group_boundary(true));
    let expected = [
        ([0, 1], GroupBoundary, Hidden),
        ([0, 3], GroupBoundary, Hidden),
        ([0, 4], GroupBoundary, Hidden),
        ([1, 2], Silhouette, Visible),
        ([1, 5], Silhouette, Visible),
        ([2, 3], Silhouette, Visible),
        ([2, 6], GroupBoundary, Visible),
        ([3, 7], Silhouette, Visible),
        ([4, 5], Silhouette, Visible),
        ([4, 7], Silhouette, Visible),
        ([5, 6], GroupBoundary, Visible),
        ([6, 7], GroupBoundary, Visible),
    ];
    assert_eq!(found, expected);
}

#[test]
fn a_selected_crease_takes_precedence_over_a_group_boundary() {
    let with_creases = cube_kinds(EdgeSelection::new(Some(0.5)).with_group_boundary(true));
    let creases_only = cube_kinds(EdgeSelection::new(Some(0.5)));
    assert_eq!(with_creases, creases_only);
    assert_eq!(
        with_creases
            .iter()
            .filter(|(_, kind, _)| *kind == EdgeKind::Crease)
            .count(),
        6
    );
    assert!(with_creases
        .iter()
        .all(|(_, kind, _)| *kind != EdgeKind::GroupBoundary));
}

#[test]
fn a_silhouette_is_not_a_group_boundary_even_when_silhouettes_are_off() {
    let found = cube_kinds(
        EdgeSelection::new(None)
            .with_silhouette(false)
            .with_group_boundary(true),
    );
    assert_eq!(found.len(), 6);
    assert!(found
        .iter()
        .all(|(_, kind, _)| *kind == EdgeKind::GroupBoundary));
}

fn draw_plan(
    meshes: &[HiddenLineMesh<'_>],
    selection: EdgeSelection,
) -> Result<
    axiolid_mesh_hidden_lines_contract::HiddenLineOutcome,
    axiolid_mesh_hidden_lines_contract::HiddenLineError,
> {
    common::try_draw_with(&ScalarHiddenLines, meshes, plan(), selection)
}
