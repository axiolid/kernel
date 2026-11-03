//! Depth clipping and cropping of visible faces: the near plane's cap, the far
//! plane, a cut sheet and the crop box.

mod common;

use axiolid_core::{Aabb2, Point2, Point3, Vec3};
use axiolid_hidden_lines::ScalarVisibleFaces;
use axiolid_mesh::TriMesh;
use axiolid_mesh_hidden_lines_contract::{
    FaceGroup, FaceSource, HiddenLineMesh, MeshId, MeshRole, OrthographicView,
};
use common::{box_mesh, faces, frame, plan_frame, source_area, try_faces};

const GROUP: FaceSource = FaceSource::Group(FaceGroup(0));

fn meshes<'a>(parts: &'a [(MeshRole, TriMesh)]) -> Vec<HiddenLineMesh<'a>> {
    parts
        .iter()
        .enumerate()
        .map(|(i, (role, mesh))| HiddenLineMesh::new(MeshId(i as u64 + 1), *role, mesh))
        .collect()
}

fn sheet(x0: f64, y0: f64, x1: f64, y1: f64, z: f64) -> TriMesh {
    TriMesh::new(
        vec![
            Point3::new(x0, y0, z),
            Point3::new(x1, y0, z),
            Point3::new(x1, y1, z),
            Point3::new(x0, y1, z),
        ],
        vec![0, 1, 2, 0, 2, 3],
    )
}

#[test]
fn a_near_plane_through_a_box_in_plan_gives_its_cap_alone() {
    let parts = [(MeshRole::Solid, box_mesh([0.0; 3], [4.0, 4.0, 2.0]))];
    let outcome = faces(
        &ScalarVisibleFaces::new(),
        &meshes(&parts),
        OrthographicView::new(plan_frame()).with_near(1.0),
    );
    assert_eq!(outcome.regions.len(), 1);
    assert_eq!(outcome.evidence.cap_regions, 1);
    assert_eq!(source_area(&outcome, 1, FaceSource::Cap), 16.0);
    assert!(outcome.regions[0].triangles.is_empty());
}

#[test]
fn a_near_plane_through_a_box_in_elevation_gives_a_cap_and_the_wall_behind_it() {
    // Looking along +y from -y: image x, up z, depth towards the viewer is -y.
    let view = OrthographicView::new(frame(Vec3::X, Vec3::Z, -Vec3::Y, Point3::ZERO));
    // Depth W = -y. The box spans y in [1, 3], so W in [-3, -1]; near = -2
    // keeps y >= 2 and cuts the box in half.
    let parts = [(MeshRole::Solid, box_mesh([0.0, 1.0, 0.0], [4.0, 3.0, 2.0]))];
    let outcome = faces(
        &ScalarVisibleFaces::new(),
        &meshes(&parts),
        view.with_near(-2.0),
    );
    assert_eq!(source_area(&outcome, 1, FaceSource::Cap), 8.0);
    // The far wall (y = 3) lies behind the cap and is hidden by it.
    assert_eq!(source_area(&outcome, 1, GROUP), 0.0);
}

#[test]
fn the_far_plane_removes_what_lies_behind_it() {
    let parts = [(MeshRole::Sheet, sheet(0.0, 0.0, 2.0, 2.0, 0.0))];
    let kept = faces(
        &ScalarVisibleFaces::new(),
        &meshes(&parts),
        OrthographicView::new(plan_frame()).with_far(0.0),
    );
    assert_eq!(source_area(&kept, 1, GROUP), 4.0);
    let gone = faces(
        &ScalarVisibleFaces::new(),
        &meshes(&parts),
        OrthographicView::new(plan_frame()).with_far(0.5),
    );
    assert!(gone.regions.is_empty());
    assert_eq!(gone.evidence.removed_triangles, 2);
}

#[test]
fn a_sheet_cut_by_the_near_plane_keeps_the_part_below_it() {
    // A ramp rising from z = 0 at y = 0 to z = 2 at y = 4; near = 1 keeps
    // y <= 2, half of its image.
    let ramp = TriMesh::new(
        vec![
            Point3::new(0.0, 0.0, 0.0),
            Point3::new(4.0, 0.0, 0.0),
            Point3::new(4.0, 4.0, 2.0),
            Point3::new(0.0, 4.0, 2.0),
        ],
        vec![0, 1, 2, 0, 2, 3],
    );
    let parts = [(MeshRole::Sheet, ramp)];
    let outcome = faces(
        &ScalarVisibleFaces::new(),
        &meshes(&parts),
        OrthographicView::new(plan_frame()).with_near(1.0),
    );
    assert_eq!(source_area(&outcome, 1, GROUP), 8.0);
    assert_eq!(outcome.evidence.cap_regions, 0, "a sheet has no cap");
}

#[test]
fn a_cut_solid_ties_with_a_sheet_in_the_near_plane() {
    let parts = [
        (MeshRole::Solid, box_mesh([0.0; 3], [4.0, 4.0, 2.0])),
        (MeshRole::Sheet, sheet(1.0, 1.0, 3.0, 3.0, 1.0)),
    ];
    let outcome = faces(
        &ScalarVisibleFaces::new(),
        &meshes(&parts),
        OrthographicView::new(plan_frame()).with_near(1.0),
    );
    // Equal depth never occludes: the sheet and the cap both own the square.
    assert_eq!(source_area(&outcome, 1, FaceSource::Cap), 16.0);
    assert_eq!(source_area(&outcome, 2, GROUP), 4.0);
    assert!(outcome.evidence.tied_sides > 0);
}

#[test]
fn the_crop_clips_every_region_to_the_box() {
    let parts = [(MeshRole::Sheet, sheet(0.0, 0.0, 4.0, 4.0, 0.0))];
    let region = Aabb2 {
        min: Point2::new(1.0, 1.0),
        max: Point2::new(3.0, 2.0),
    };
    let outcome = try_faces(
        &ScalarVisibleFaces::new(),
        &meshes(&parts),
        OrthographicView::new(plan_frame()),
        Some(region),
    )
    .expect("faces");
    assert_eq!(source_area(&outcome, 1, GROUP), 2.0);
    for point in outcome.regions.iter().flat_map(|r| &r.outer) {
        assert!((1.0..=3.0).contains(&point.x) && (1.0..=2.0).contains(&point.y));
    }
}
