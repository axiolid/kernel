//! Faces selected by normal: area and extent of a floor, a wall side and a
//! slope, with the tolerance's angular slack and the refusals.

use axiolid_core::{Interval, Point3, Tolerance, Vec3};
use axiolid_measure::{
    extent_along, select_by_normal, surface_properties, FaceSelectionError, NormalCone,
};
use axiolid_mesh::TriMesh;

const TOL: Tolerance = Tolerance::METRE;

/// Outward-wound faces of a box: bottom (triangles 0, 1), top (2, 3), front
/// (-Y, 4, 5), right (+X, 6, 7), back (+Y, 8, 9), left (-X, 10, 11).
const BOX: [u32; 36] = [
    0, 2, 1, 0, 3, 2, 4, 5, 6, 4, 6, 7, 0, 1, 5, 0, 5, 4, 1, 2, 6, 1, 6, 5, 2, 3, 7, 2, 7, 6, 3, 0,
    4, 3, 4, 7,
];

fn box_at(origin: [f64; 3], size: [f64; 3]) -> TriMesh {
    let corners = [
        [0.0, 0.0, 0.0],
        [1.0, 0.0, 0.0],
        [1.0, 1.0, 0.0],
        [0.0, 1.0, 0.0],
        [0.0, 0.0, 1.0],
        [1.0, 0.0, 1.0],
        [1.0, 1.0, 1.0],
        [0.0, 1.0, 1.0],
    ];
    let positions = corners
        .iter()
        .map(|c| {
            Point3::new(
                origin[0] + c[0] * size[0],
                origin[1] + c[1] * size[1],
                origin[2] + c[2] * size[2],
            )
        })
        .collect();
    TriMesh::new(positions, BOX.to_vec())
}

fn cone(axis: [f64; 3], half_angle: f64) -> NormalCone {
    NormalCone::new(Vec3::new(axis[0], axis[1], axis[2]), half_angle).unwrap()
}

fn down() -> NormalCone {
    cone([0.0, 0.0, -1.0], 0.0)
}

const X: Vec3 = Vec3::new(1.0, 0.0, 0.0);
const Y: Vec3 = Vec3::new(0.0, 1.0, 0.0);
const Z: Vec3 = Vec3::new(0.0, 0.0, 1.0);

#[test]
fn the_downward_faces_of_a_box_are_its_base() {
    let mesh = box_at([0.0; 3], [2.0, 3.0, 5.0]);
    let floor = select_by_normal(&mesh, down(), TOL).unwrap();
    assert_eq!(floor.triangles, vec![0, 1]);
    assert_eq!(floor.area, 6.0);
    assert_eq!(
        extent_along(&mesh, &floor.triangles, X).unwrap(),
        Some(Interval::new(0.0, 2.0))
    );
}

#[test]
fn the_axis_has_a_sign() {
    let mesh = box_at([0.0; 3], [2.0, 3.0, 5.0]);
    let up = select_by_normal(&mesh, cone([0.0, 0.0, 1.0], 0.0), TOL).unwrap();
    assert_eq!(up.triangles, vec![2, 3]);
}

#[test]
fn a_wall_side_has_its_area_and_extents() {
    let mesh = box_at([0.0; 3], [2.0, 3.0, 5.0]);
    let side = select_by_normal(&mesh, cone([1.0, 0.0, 0.0], 0.0), TOL).unwrap();
    assert_eq!(side.triangles, vec![6, 7]);
    assert_eq!(side.area, 15.0);
    let along_y = Some(Interval::new(0.0, 3.0));
    assert_eq!(extent_along(&mesh, &side.triangles, Y).unwrap(), along_y);
    assert_eq!(
        extent_along(&mesh, &side.triangles, Z).unwrap(),
        Some(Interval::new(0.0, 5.0))
    );
    // The direction is normalised, so its length does not matter.
    assert_eq!(
        extent_along(&mesh, &side.triangles, Vec3::new(0.0, 2.0, 0.0)).unwrap(),
        along_y
    );
}

#[test]
fn a_selection_with_no_face_has_no_area_and_no_extent() {
    let mesh = box_at([0.0; 3], [2.0, 3.0, 5.0]);
    let none = select_by_normal(&mesh, cone([1.0, 1.0, 0.0], 0.1), TOL).unwrap();
    assert!(none.triangles.is_empty());
    assert_eq!(none.area, 0.0);
    assert_eq!(extent_along(&mesh, &none.triangles, X).unwrap(), None);
}

/// One quad (two triangles) rising along +Y at `slope` radians over a run
/// of `run`, width 2: its normal is `slope` away from +Z.
fn ramp(slope: f64, run: f64) -> TriMesh {
    let (y, z) = (run, run * slope.tan());
    TriMesh::new(
        vec![
            Point3::new(0.0, 0.0, 0.0),
            Point3::new(2.0, 0.0, 0.0),
            Point3::new(2.0, y, z),
            Point3::new(0.0, y, z),
        ],
        vec![0, 1, 2, 0, 2, 3],
    )
}

#[test]
fn a_slope_is_selected_only_within_the_half_angle() {
    let slope = 30.0_f64.to_radians();
    let mesh = ramp(slope, 4.0);
    let up = [0.0, 0.0, 1.0];
    let miss = select_by_normal(&mesh, cone(up, 29.0_f64.to_radians()), TOL).unwrap();
    assert!(miss.triangles.is_empty());
    let hit = select_by_normal(&mesh, cone(up, 31.0_f64.to_radians()), TOL).unwrap();
    assert_eq!(hit.triangles, vec![0, 1]);
    // Width 2 times a slant length of run / cos(slope).
    assert!((hit.area - 2.0 * 4.0 / slope.cos()).abs() < 1e-12);
}

fn tilted_floor() -> TriMesh {
    TriMesh::new(
        vec![
            Point3::new(0.0, 0.0, 0.0),
            Point3::new(1.0, 0.0, 0.0),
            Point3::new(1.0, 1.0, 1e-12),
            Point3::new(0.0, 1.0, 1e-12),
        ],
        vec![0, 1, 2, 0, 2, 3],
    )
}

#[test]
fn the_angular_tolerance_is_slack_on_the_half_angle() {
    let up = cone([0.0, 0.0, 1.0], 0.0);
    let floor = tilted_floor();
    assert_eq!(
        select_by_normal(&floor, up, TOL).unwrap().triangles,
        vec![0, 1]
    );
    let no_slack = Tolerance::new(1e-6, 0.0).unwrap();
    assert!(select_by_normal(&floor, up, no_slack)
        .unwrap()
        .triangles
        .is_empty());
}

#[test]
fn an_exactly_level_floor_is_selected_even_without_slack() {
    let no_slack = Tolerance::new(1e-6, 0.0).unwrap();
    let mesh = box_at([0.0; 3], [2.0, 3.0, 5.0]);
    assert_eq!(
        select_by_normal(&mesh, down(), no_slack).unwrap().triangles,
        vec![0, 1]
    );
}

#[test]
fn a_half_angle_of_pi_reproduces_the_surface_area_bit_for_bit() {
    let all = cone([0.0, 0.0, 1.0], core::f64::consts::PI);
    for mesh in [
        box_at([0.0; 3], [2.0, 3.0, 5.0]),
        box_at([6.0e5, 5.6e6, 204.5], [2.0, 3.0, 5.0]),
        ramp(0.3, 4.0),
    ] {
        let selection = select_by_normal(&mesh, all, TOL).unwrap();
        assert_eq!(selection.triangles.len(), mesh.indices.len() / 3);
        let want = surface_properties(&mesh, TOL).unwrap().area;
        assert_eq!(selection.area.to_bits(), want.to_bits());
    }
}

#[test]
fn a_box_on_a_national_grid_has_its_closed_forms() {
    let mesh = box_at([6.0e5, 5.6e6, 204.5], [2.0, 3.0, 5.0]);
    let floor = select_by_normal(&mesh, down(), TOL).unwrap();
    assert_eq!(floor.triangles, vec![0, 1]);
    assert!((floor.area - 6.0).abs() < 1e-9);
    assert_eq!(
        extent_along(&mesh, &floor.triangles, X).unwrap(),
        Some(Interval::new(6.0e5, 6.0e5 + 2.0))
    );
    let side = select_by_normal(&mesh, cone([1.0, 0.0, 0.0], 0.0), TOL).unwrap();
    assert!((side.area - 15.0).abs() < 1e-9);
    assert_eq!(
        extent_along(&mesh, &side.triangles, Z).unwrap(),
        Some(Interval::new(204.5, 209.5))
    );
}

#[test]
fn a_degenerate_triangle_refuses_the_mesh() {
    let mut mesh = box_at([0.0; 3], [2.0, 3.0, 5.0]);
    mesh.indices.extend([0, 1, 1]);
    assert!(matches!(
        select_by_normal(&mesh, down(), TOL),
        Err(FaceSelectionError::MeshNotSurfaceUsable(_))
    ));
}

#[test]
fn a_cone_is_validated_and_its_axis_normalised() {
    for axis in [
        Vec3::ZERO,
        Vec3::new(f64::NAN, 0.0, 1.0),
        Vec3::new(f64::INFINITY, 0.0, 0.0),
    ] {
        assert_eq!(
            NormalCone::new(axis, 0.1),
            Err(FaceSelectionError::InvalidAxis)
        );
    }
    for half_angle in [-0.1, core::f64::consts::PI + 1e-9, f64::NAN, f64::INFINITY] {
        assert_eq!(
            NormalCone::new(Z, half_angle),
            Err(FaceSelectionError::InvalidHalfAngle)
        );
    }
    let cone = NormalCone::new(Vec3::new(0.0, 0.0, -10.0), 0.5).unwrap();
    assert_eq!(cone.axis(), Vec3::new(0.0, 0.0, -1.0));
    assert_eq!(cone.half_angle(), 0.5);
    assert!(NormalCone::new(Z, core::f64::consts::PI).is_ok());
}

#[test]
fn an_extent_refuses_a_bad_direction_and_bad_triangles() {
    let mesh = box_at([0.0; 3], [2.0, 3.0, 5.0]);
    for along in [
        Vec3::ZERO,
        Vec3::new(f64::NAN, 0.0, 0.0),
        Vec3::new(0.0, f64::INFINITY, 0.0),
    ] {
        assert_eq!(
            extent_along(&mesh, &[0], along),
            Err(FaceSelectionError::InvalidDirection)
        );
    }
    assert_eq!(
        extent_along(&mesh, &[0, 12], X),
        Err(FaceSelectionError::TriangleOutOfRange { triangle: 12 })
    );
    let mut bad = box_at([0.0; 3], [2.0, 3.0, 5.0]);
    bad.indices[3] = 99;
    assert_eq!(
        extent_along(&bad, &[0, 1], X),
        Err(FaceSelectionError::IndexOutOfRange { triangle: 1 })
    );
    let mut nan = box_at([0.0; 3], [2.0, 3.0, 5.0]);
    nan.positions[1] = Point3::new(f64::NAN, 0.0, 0.0);
    assert_eq!(
        extent_along(&nan, &[0], X),
        Err(FaceSelectionError::NonFinitePosition { triangle: 0 })
    );
    assert!(matches!(
        select_by_normal(&nan, down(), TOL),
        Err(FaceSelectionError::MeshNotSurfaceUsable(_))
    ));
}

#[test]
fn an_extent_ignores_order_and_duplicates() {
    let mesh = box_at([0.0; 3], [2.0, 3.0, 5.0]);
    assert_eq!(
        extent_along(&mesh, &[7, 6, 6, 7], Z).unwrap(),
        extent_along(&mesh, &[6, 7], Z).unwrap()
    );
}

fn huge_triangle() -> TriMesh {
    TriMesh::new(
        vec![
            Point3::new(0.0, 0.0, 0.0),
            Point3::new(0.0, 1e80, 0.0),
            Point3::new(0.0, 0.0, 1e80),
        ],
        vec![0, 1, 2],
    )
}

#[test]
fn an_overflowing_normal_or_area_is_refused_not_summed_or_dropped() {
    let mesh = huge_triangle();
    // The area would be infinite.
    assert_eq!(
        select_by_normal(&mesh, cone([1.0, 0.0, 0.0], 0.0), TOL),
        Err(FaceSelectionError::NonFiniteNormal { triangle: 0 })
    );
    // The face would be dropped by a wrong angle.
    assert_eq!(
        select_by_normal(
            &mesh,
            cone([1.0, 1.0, 0.0], core::f64::consts::FRAC_PI_3),
            TOL
        ),
        Err(FaceSelectionError::NonFiniteNormal { triangle: 0 })
    );
}

#[test]
fn directions_of_any_finite_magnitude_are_normalised() {
    for scale in [1e-170, 1e-300, 1e200, 1e300] {
        let cone = NormalCone::new(Vec3::new(scale, 0.0, 0.0), 0.0).unwrap();
        assert_eq!(cone.axis(), X);
        let mesh = box_at([0.0; 3], [2.0, 3.0, 5.0]);
        assert_eq!(
            extent_along(&mesh, &[6, 7], Vec3::new(0.0, scale, 0.0)).unwrap(),
            Some(Interval::new(0.0, 3.0))
        );
    }
}
