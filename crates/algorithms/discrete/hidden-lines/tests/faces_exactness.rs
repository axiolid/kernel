//! Region vertices are the exact points rounded once, so they equal the
//! hidden-line endpoints of the same points bit for bit.

mod common;

use axiolid_core::{Point2, Point3, Vec3};
use axiolid_hidden_lines::ScalarVisibleFaces;
use axiolid_mesh_hidden_lines_contract::{OrthographicView, Visibility};
use common::{box_mesh, draw, faces, frame, isometric_frame, scalar, solids};

fn same(a: Point2, b: Point2) -> bool {
    a.x.to_bits() == b.x.to_bits() && a.y.to_bits() == b.y.to_bits()
}

/// Every ring vertex of every region is an endpoint of a visible hidden-line
/// segment of the same scene, in plan and isometric views.
#[test]
fn region_vertices_equal_visible_hidden_line_endpoints() {
    // Two boxes apart in space whose images overlap, so edge crossings appear.
    let boxes = [
        box_mesh([0.0, 0.0, 0.0], [2.0, 2.0, 2.0]),
        box_mesh([1.0, 1.0, 3.0], [3.0, 3.0, 5.0]),
        box_mesh([-1.0, 1.5, 0.5], [0.5, 2.5, 1.5]),
    ];
    let meshes = solids(&boxes);
    let o = Point3::new(0.1, 0.2, 0.3);
    for (name, frame) in [
        ("plan", common::plan_frame()),
        ("isometric", isometric_frame(o)),
        ("front", frame(Vec3::X, Vec3::Z, -Vec3::Y, o)),
    ] {
        let view = OrthographicView::new(frame);
        let lines = draw(&scalar(), &meshes, view);
        let endpoints: Vec<Point2> = lines
            .segments
            .iter()
            .filter(|s| s.visibility == Visibility::Visible)
            .flat_map(|s| [s.start, s.end])
            .collect();
        let outcome = faces(&ScalarVisibleFaces::new(), &meshes, view);
        assert!(!outcome.regions.is_empty(), "{name}");
        for region in &outcome.regions {
            for vertex in region.outer.iter().chain(region.holes.iter().flatten()) {
                assert!(
                    endpoints.iter().any(|e| same(*e, *vertex)),
                    "{name}: {vertex:?} is no visible endpoint"
                );
            }
        }
    }
}

/// Neighbouring regions share their common vertices bit for bit: every
/// vertex of a region is a vertex of a neighbour or on the outer outline.
#[test]
fn two_boxes_in_plan_share_the_vertices_of_their_common_boundary() {
    let boxes = [
        box_mesh([0.0, 0.0, 2.0], [2.0, 2.0, 3.0]),
        box_mesh([1.0, 1.0, 0.0], [3.0, 3.0, 1.0]),
    ];
    let outcome = faces(
        &ScalarVisibleFaces::new(),
        &solids(&boxes),
        OrthographicView::new(common::plan_frame()),
    );
    // The front box's top is whole; the back box's top is the L around it,
    // whose inner corner (2, 1) and (1, 2) are shared with the front box.
    let front = &outcome.regions[0];
    let back = &outcome.regions[1];
    for corner in [Point2::new(2.0, 1.0), Point2::new(1.0, 2.0)] {
        assert!(back.outer.iter().any(|p| same(*p, corner)));
    }
    assert!(front.outer.iter().any(|p| same(*p, Point2::new(2.0, 2.0))));
}
