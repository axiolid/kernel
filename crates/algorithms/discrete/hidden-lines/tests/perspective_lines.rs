//! Perspective hidden lines on hand-checked scenes (ADR 0107). The eye is the
//! origin of the view frame and looks along `-z`; the image of `(U, V, W)` is
//! `(U / -W, V / -W)`.

mod common;

use axiolid_contracts::ExecutionOptions;
use axiolid_core::{Aabb2, Point2, Point3, Tolerance, Vec3};
use axiolid_hidden_lines::{ScalarHiddenLines, ScanHiddenLines};
use axiolid_mesh::TriMesh;
use axiolid_mesh_hidden_lines_contract::{
    EdgeKind, EdgeSelection, HiddenLineMesh, HiddenLineRequest, HiddenLineSegment, MeshHiddenLines,
    MeshId, MeshRole, View, Visibility,
};
use common::{box_mesh, draw, plan_frame, solids};

fn options() -> ExecutionOptions {
    ExecutionOptions::new(Tolerance::new(1e-9, 1e-9).expect("tolerance"))
}

fn perspective(near: f64) -> View {
    View::perspective(plan_frame(), near).expect("valid")
}

fn segments(meshes: &[TriMesh], view: View) -> Vec<HiddenLineSegment> {
    let input = solids(meshes);
    let scalar = draw(&ScalarHiddenLines::new(), &input, view);
    // The linear-scan baseline draws the same.
    let scan = draw(&ScanHiddenLines, &input, view);
    assert_eq!(scalar.segments, scan.segments);
    scalar.segments
}

fn count(segments: &[HiddenLineSegment], kind: EdgeKind, visibility: Visibility) -> usize {
    segments
        .iter()
        .filter(|s| s.kind == kind && s.visibility == visibility)
        .count()
}

/// Whether the segment runs between the two image points, in either order.
fn joins(segment: &HiddenLineSegment, a: [f64; 2], b: [f64; 2]) -> bool {
    let (p, q) = (Point2::new(a[0], a[1]), Point2::new(b[0], b[1]));
    (segment.start == p && segment.end == q) || (segment.start == q && segment.end == p)
}

#[test]
fn a_box_in_front_of_the_eye_hides_its_back_edges_and_depth_edges() {
    // [-1, 1]^2 x [-4, -2], near -1: the front face is at depth 2 (image
    // +-1/2), the back face at depth 4 (image +-1/4).
    let boxes = [box_mesh([-1.0, -1.0, -4.0], [1.0, 1.0, -2.0])];
    let drawn = segments(&boxes, perspective(-1.0));
    assert_eq!(drawn.len(), 12);
    assert_eq!(
        count(&drawn, EdgeKind::Silhouette, Visibility::Visible),
        4,
        "{drawn:?}"
    );
    assert_eq!(count(&drawn, EdgeKind::Crease, Visibility::Hidden), 8);
    // Every endpoint is a corner of the front square (+-1/2) or of the back
    // square (+-1/4), exactly.
    for segment in &drawn {
        for p in [segment.start, segment.end] {
            assert!(
                [0.5, 0.25].contains(&p.x.abs()) && [0.5, 0.25].contains(&p.y.abs()),
                "{segment:?}"
            );
        }
    }
    // The front ring, edge by edge: corners (+-0.5, +-0.5).
    let ring: Vec<_> = drawn
        .iter()
        .filter(|s| s.visibility == Visibility::Visible)
        .collect();
    let expected = [
        ([-0.5, -0.5], [0.5, -0.5]),
        ([0.5, -0.5], [0.5, 0.5]),
        ([0.5, 0.5], [-0.5, 0.5]),
        ([-0.5, 0.5], [-0.5, -0.5]),
    ];
    for (a, b) in expected {
        assert_eq!(ring.iter().filter(|s| joins(s, a, b)).count(), 1);
    }
}

#[test]
fn an_endpoint_is_the_exact_division_rounded_once() {
    // Depth 3 and depth 7: the images 1/3 and 1/7 are not dyadic.
    let boxes = [box_mesh([-1.0, -1.0, -7.0], [1.0, 1.0, -3.0])];
    let drawn = segments(&boxes, perspective(-1.0));
    for segment in &drawn {
        for p in [segment.start, segment.end] {
            assert!(
                [1.0 / 3.0, -1.0 / 3.0, 1.0 / 7.0, -1.0 / 7.0].contains(&p.x),
                "{segment:?}"
            );
        }
    }
}

#[test]
fn the_cap_of_a_cut_box_hides_everything_behind_it() {
    // Cut at depth 3 the box keeps depths 3 to 4; its cap (image +-1/3)
    // covers the back ring (image +-1/4) and the kept depth edges.
    let boxes = [box_mesh([-1.0, -1.0, -4.0], [1.0, 1.0, -2.0])];
    let drawn = segments(&boxes, perspective(-3.0));
    assert_eq!(drawn.len(), 8, "{drawn:?}");
    assert!(drawn.iter().all(|s| s.visibility == Visibility::Hidden));
}

#[test]
fn an_eye_inside_a_box_sees_the_cap_and_nothing_of_the_far_wall() {
    // The box surrounds the eye; near -0.5 keeps depths 0.5 to 1.
    let boxes = [box_mesh([-1.0, -1.0, -1.0], [1.0, 1.0, 1.0])];
    let drawn = segments(&boxes, perspective(-0.5));
    assert_eq!(drawn.len(), 8, "{drawn:?}");
    assert!(drawn.iter().all(|s| s.visibility == Visibility::Hidden));
    // Whatever lies behind the eye is never drawn.
    for segment in &drawn {
        for p in [segment.start, segment.end] {
            assert!(p.x.abs() <= 2.0 && p.y.abs() <= 2.0, "{segment:?}");
        }
    }
}

#[test]
fn a_box_behind_the_eye_or_beyond_the_near_plane_draws_nothing() {
    let behind = [box_mesh([-1.0, -1.0, 2.0], [1.0, 1.0, 4.0])];
    assert!(segments(&behind, perspective(-1.0)).is_empty());
    let beyond = [box_mesh([-1.0, -1.0, -0.9], [1.0, 1.0, -0.2])];
    assert!(segments(&beyond, perspective(-1.0)).is_empty());
}

#[test]
fn a_nearer_box_hides_the_edge_it_covers_at_the_projected_boundary() {
    // A small cube at depth 1.5 to 2.5 in front of the big box's back wall:
    // the images are exact, so the covered part of the big box's front-face
    // edge ends where the small cube's silhouette crosses it.
    let big = box_mesh([-3.0, -3.0, -9.0], [3.0, 3.0, -6.0]);
    let small = box_mesh([-1.0, -1.0, -4.0], [1.0, 1.0, -3.0]);
    let drawn = segments(&[big, small], perspective(-1.0));
    let hidden = drawn
        .iter()
        .filter(|s| s.source.mesh.0 == 1 && s.visibility == Visibility::Hidden)
        .count();
    assert!(hidden > 0, "{drawn:?}");
    // No visible piece of the big box lies inside the small cube's image.
    let inside = |p: Point2| p.x.abs() < 1.0 / 3.0 && p.y.abs() < 1.0 / 3.0;
    for s in drawn
        .iter()
        .filter(|s| s.source.mesh.0 == 1 && s.visibility == Visibility::Visible)
    {
        let mid = Point2::new((s.start.x + s.end.x) / 2.0, (s.start.y + s.end.y) / 2.0);
        assert!(!inside(mid), "{s:?}");
    }
}

#[test]
fn a_view_that_is_not_in_front_of_the_eye_is_refused() {
    let boxes = [box_mesh([-1.0, -1.0, -4.0], [1.0, 1.0, -2.0])];
    let input = solids(&boxes);
    let mut view = perspective(-1.0);
    for near in [
        None,
        Some(0.0),
        Some(1.0),
        Some(f64::NAN),
        Some(f64::NEG_INFINITY),
    ] {
        view.near = near;
        let request =
            HiddenLineRequest::new(&input, view, EdgeSelection::new(None), common::LIMITS);
        let error = ScalarHiddenLines::new()
            .hidden_lines(&request, &options())
            .expect_err("refused");
        assert!(
            matches!(error.error, axiolid_contracts::GeomError::InvalidInput(_)),
            "{near:?}: {error}"
        );
    }
}

#[test]
fn a_rotated_eye_draws_the_same_picture_as_the_rotated_scene() {
    // Looking along -x from the origin at a box on the -x axis equals looking
    // along -z at the box rotated onto the -z axis.
    let along_z = [box_mesh([-1.0, -1.0, -4.0], [1.0, 1.0, -2.0])];
    let mut rotated = box_mesh([-1.0, -1.0, -4.0], [1.0, 1.0, -2.0]);
    for p in &mut rotated.positions {
        // (x, y, z) -> (z, x, y): a proper rotation taking -z to -x.
        *p = Point3::new(p.z, p.x, p.y);
    }
    let frame = common::frame(Vec3::Y, Vec3::Z, Vec3::X, Point3::ZERO);
    let view = View::perspective(frame, -1.0).expect("valid");
    let a = segments(&along_z, perspective(-1.0));
    let b = segments(&[rotated], view);
    assert_eq!(a.len(), b.len());
    for (x, y) in a.iter().zip(&b) {
        assert_eq!((x.visibility, x.kind), (y.visibility, y.kind));
        assert_eq!([x.start, x.end], [y.start, y.end], "{x:?} against {y:?}",);
    }
}

#[test]
fn a_far_depth_clips_the_edges_at_the_exact_projected_point() {
    // Far depth 3 keeps the front ring (depth 2) and the first unit of the
    // depth edges, which end at (+-1/3, +-1/3); the back ring is clipped.
    let boxes = [box_mesh([-1.0, -1.0, -4.0], [1.0, 1.0, -2.0])];
    let drawn = segments(&boxes, perspective(-1.0).with_far(-3.0));
    assert_eq!(drawn.len(), 8, "{drawn:?}");
    assert_eq!(count(&drawn, EdgeKind::Silhouette, Visibility::Visible), 4);
    assert_eq!(count(&drawn, EdgeKind::Crease, Visibility::Hidden), 4);
    for segment in drawn.iter().filter(|s| s.kind == EdgeKind::Crease) {
        let ends = [segment.start, segment.end];
        assert!(ends
            .iter()
            .any(|p| p.x.abs() == 1.0 / 3.0 && p.y.abs() == 1.0 / 3.0));
        assert!(ends.iter().any(|p| p.x.abs() == 0.5 && p.y.abs() == 0.5));
    }
}

#[test]
fn a_crop_region_is_in_image_coordinates() {
    let boxes = [box_mesh([-1.0, -1.0, -4.0], [1.0, 1.0, -2.0])];
    let input = solids(&boxes);
    let region = Aabb2 {
        min: Point2::new(-0.1, -1.0),
        max: Point2::new(0.3, 0.4),
    };
    let request = HiddenLineRequest::new(
        &input,
        perspective(-1.0),
        EdgeSelection::new(Some(0.5)),
        common::LIMITS,
    )
    .with_region(region);
    let drawn = ScalarHiddenLines::new()
        .hidden_lines(&request, &options())
        .expect("draws")
        .segments;
    // Visible: the lower front edge, y = -1/2, between x = -0.1 and 0.3.
    let visible: Vec<_> = drawn
        .iter()
        .filter(|s| s.visibility == Visibility::Visible)
        .collect();
    assert_eq!(visible.len(), 1, "{drawn:?}");
    assert!(joins(visible[0], [-0.1, -0.5], [0.3, -0.5]));
    // Hidden: the back ring's two horizontal edges clipped to x >= -0.1
    // (0.35 each), its right edge (0.5), and the two depth edges at x = +1/2
    // clipped to x <= 0.3, which keep 0.05 * sqrt(2) each.
    let hidden: f64 = drawn
        .iter()
        .filter(|s| s.visibility == Visibility::Hidden)
        .map(|s| (s.end - s.start).length())
        .sum();
    assert!(
        (hidden - (1.2 + 0.1 * 2f64.sqrt())).abs() < 1e-12,
        "{hidden}"
    );
    for segment in &drawn {
        for p in [segment.start, segment.end] {
            assert!(region.contains(p), "{segment:?}");
        }
    }
}

fn sheet(corners: [[f64; 3]; 4]) -> TriMesh {
    TriMesh::new(
        corners
            .iter()
            .map(|c| Point3::new(c[0], c[1], c[2]))
            .collect(),
        vec![0, 1, 2, 0, 2, 3],
    )
}

/// The pieces of the far sheet's edge (0, 1), from `y = -3` to `3` at depth
/// 9, as (visibility, start y, end y).
///
/// Both sheets lie in the plane `x = 0`, which contains the eye, so they are
/// seen edge-on. The near sheet spans depth 5 to 6 and `|y| <= 1`: it covers
/// the image height `|y| <= 1/5` (by its edge at depth 5) or, with a near
/// depth of 5.5 that cuts the sheet, `|y| <= 1/5.5` (its points at the plane).
fn far_edge(near: f64) -> Vec<(Visibility, f64, f64)> {
    let near_sheet = sheet([
        [0.0, -1.0, -6.0],
        [0.0, 1.0, -6.0],
        [0.0, 1.0, -5.0],
        [0.0, -1.0, -5.0],
    ]);
    let far_sheet = sheet([
        [0.0, -3.0, -9.0],
        [0.0, 3.0, -9.0],
        [0.0, 3.0, -9.5],
        [0.0, -3.0, -9.5],
    ]);
    let meshes = [
        HiddenLineMesh::new(MeshId(1), MeshRole::Sheet, &near_sheet),
        HiddenLineMesh::new(MeshId(2), MeshRole::Sheet, &far_sheet),
    ];
    let request = HiddenLineRequest::new(
        &meshes,
        perspective(near),
        EdgeSelection::new(None),
        common::LIMITS,
    );
    ScalarHiddenLines::new()
        .hidden_lines(&request, &options())
        .expect("draws")
        .segments
        .iter()
        .filter(|s| s.source.mesh.0 == 2 && s.source.vertices == [0, 1])
        .map(|s| {
            assert_eq!((s.start.x, s.end.x), (0.0, 0.0));
            (s.visibility, s.start.y, s.end.y)
        })
        .collect()
}

#[test]
fn an_edge_on_sheet_through_the_eye_hides_exactly_its_image() {
    use Visibility::{Hidden, Visible};
    let third = 1.0 / 3.0;
    assert_eq!(
        far_edge(-1.0),
        [
            (Visible, -third, -0.2),
            (Hidden, -0.2, 0.2),
            (Visible, 0.2, third)
        ]
    );
    // Near depth 5.5 removes the sheet's edge at depth 5; what it keeps of the
    // sheet reaches |y| = 1 at depth 5.5, so it hides up to the image 1/5.5.
    let reach = 1.0 / 5.5;
    assert_eq!(
        far_edge(-5.5),
        [
            (Visible, -third, -reach),
            (Hidden, -reach, reach),
            (Visible, reach, third)
        ]
    );
}
