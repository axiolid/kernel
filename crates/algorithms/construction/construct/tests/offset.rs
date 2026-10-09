//! Solid offset and shelling (#78).
//!
//! The oracles are closed-form: a cube offset outward by `d` is a cube of
//! edge `a + 2d`, and a shell's volume is outer minus inner. Both are
//! computed from the request, not read back from the result.

use axiolid_construct::offset::{offset_solid, shell_solid, OffsetDirection};
use axiolid_construct::polyhedron::{triangulate, Polyhedron};
use axiolid_core::{Point3, Tolerance};
use axiolid_heal::mesh::MeshHealer;
use axiolid_heal::{self_intersections, Diagnose};
use axiolid_measure::volume_properties;

fn tol() -> Tolerance {
    Tolerance::new(1e-6, 1e-9).expect("tolerance")
}

fn box_solid(min: [f64; 3], max: [f64; 3]) -> Polyhedron {
    let [x0, y0, z0] = min;
    let [x1, y1, z1] = max;
    let p = |x: f64, y: f64, z: f64| Point3::new(x, y, z);
    Polyhedron::new(vec![
        vec![p(x0, y0, z0), p(x0, y1, z0), p(x1, y1, z0), p(x1, y0, z0)],
        vec![p(x0, y0, z1), p(x1, y0, z1), p(x1, y1, z1), p(x0, y1, z1)],
        vec![p(x0, y0, z0), p(x1, y0, z0), p(x1, y0, z1), p(x0, y0, z1)],
        vec![p(x0, y1, z0), p(x0, y1, z1), p(x1, y1, z1), p(x1, y1, z0)],
        vec![p(x0, y0, z0), p(x0, y0, z1), p(x0, y1, z1), p(x0, y1, z0)],
        vec![p(x1, y0, z0), p(x1, y1, z0), p(x1, y1, z1), p(x1, y0, z1)],
    ])
    .expect("box is a valid solid")
}

fn volume(solid: &Polyhedron) -> f64 {
    volume_properties(&triangulate(solid), tol())
        .expect("closed solid")
        .signed_volume
}

#[test]
fn a_cube_offset_outward_grows_by_twice_the_distance() {
    // The identity #78 names: this fails immediately if edges or corners are
    // mishandled, because a wrong corner changes the volume.
    let cube = box_solid([0.0, 0.0, 0.0], [2.0, 2.0, 2.0]);
    let d = 0.5;
    let grown = offset_solid(&cube, d, OffsetDirection::Outward).expect("offset");

    let expected = (2.0 + 2.0 * d).powi(3);
    let actual = volume(&grown);
    assert!(
        (actual - expected).abs() < 1e-9,
        "cube offset outward by {d}: expected {expected}, got {actual}"
    );
}

#[test]
fn a_cube_offset_inward_shrinks_by_twice_the_distance() {
    let cube = box_solid([0.0, 0.0, 0.0], [2.0, 2.0, 2.0]);
    let d = 0.25;
    let shrunk = offset_solid(&cube, d, OffsetDirection::Inward).expect("offset");

    let expected = (2.0 - 2.0 * d).powi(3);
    let actual = volume(&shrunk);
    assert!(
        (actual - expected).abs() < 1e-9,
        "cube offset inward by {d}: expected {expected}, got {actual}"
    );
}

#[test]
fn an_offset_result_is_closed_manifold_and_clean() {
    let cube = box_solid([0.0, 0.0, 0.0], [2.0, 2.0, 2.0]);
    for direction in [OffsetDirection::Outward, OffsetDirection::Inward] {
        let result = offset_solid(&cube, 0.3, direction).expect("offset");
        let mesh = triangulate(&result);
        let diagnosis = MeshHealer.diagnose(&mesh, tol()).expect("diagnose");
        assert!(
            diagnosis.is_clean(),
            "{direction:?} offset produced defects: {:?}",
            diagnosis.defects
        );
        assert!(
            self_intersections(&mesh).is_empty(),
            "{direction:?} offset self-intersects"
        );
    }
}

#[test]
fn an_inward_offset_past_the_half_thickness_is_refused() {
    // A 2-unit cube has half-thickness 1. Offsetting inward by 1.5 would
    // turn it inside out; the caller needs told, not handed a bad solid.
    let cube = box_solid([0.0, 0.0, 0.0], [2.0, 2.0, 2.0]);
    let error = offset_solid(&cube, 1.5, OffsetDirection::Inward)
        .expect_err("over-offset must be refused, not emitted");
    let text = format!("{error}");
    assert!(
        text.contains("half-thickness") || text.contains("passes through itself"),
        "the refusal must name the collapse, got: {text}"
    );
}

#[test]
fn an_exactly_collapsing_offset_is_refused() {
    // Exactly half-thickness: the cube shrinks to a point, zero volume. The
    // boundary is degenerate rather than reversed, which is the other
    // collapse mode.
    let cube = box_solid([0.0, 0.0, 0.0], [2.0, 2.0, 2.0]);
    assert!(
        offset_solid(&cube, 1.0, OffsetDirection::Inward).is_err(),
        "an offset that collapses the solid to nothing must be refused"
    );
}

#[test]
fn a_non_positive_distance_is_refused() {
    let cube = box_solid([0.0, 0.0, 0.0], [2.0, 2.0, 2.0]);
    assert!(offset_solid(&cube, 0.0, OffsetDirection::Outward).is_err());
    assert!(offset_solid(&cube, -1.0, OffsetDirection::Outward).is_err());
    assert!(offset_solid(&cube, f64::NAN, OffsetDirection::Outward).is_err());
}

#[test]
fn a_shelled_box_has_the_volume_of_outer_minus_inner() {
    let outer = box_solid([0.0, 0.0, 0.0], [4.0, 4.0, 4.0]);
    let thickness = 0.5;
    let shelled = shell_solid(&outer, thickness).expect("shell");

    // The cavity is the inward offset: a cube of edge 4 - 2*0.5 = 3.
    let expected = 4.0_f64.powi(3) - 3.0_f64.powi(3);
    let actual = volume(&shelled);
    assert!(
        (actual - expected).abs() < 1e-9,
        "shelled box: expected {expected}, got {actual}"
    );
}

#[test]
fn a_shelled_box_has_the_requested_wall_thickness() {
    let outer = box_solid([0.0, 0.0, 0.0], [4.0, 4.0, 4.0]);
    let thickness = 0.75;
    let shelled = shell_solid(&outer, thickness).expect("shell");

    // Measured, not assumed: the cavity's extent on each axis must be the
    // outer extent less two wall thicknesses. Reading it off the result
    // catches a shell that reported success while building the wrong wall.
    let xs: Vec<f64> = shelled
        .faces()
        .iter()
        .flat_map(|f| f.iter().map(|p| p.x))
        .collect();
    let inner_low = xs
        .iter()
        .copied()
        .filter(|&x| x > 1e-9)
        .fold(f64::INFINITY, f64::min);
    assert!(
        (inner_low - thickness).abs() < 1e-9,
        "wall thickness: expected {thickness}, measured {inner_low}"
    );
}

#[test]
fn a_thickness_that_collapses_the_cavity_is_refused() {
    let outer = box_solid([0.0, 0.0, 0.0], [2.0, 2.0, 2.0]);
    assert!(
        shell_solid(&outer, 1.5).is_err(),
        "a wall thicker than the half-thickness leaves no cavity and must refuse"
    );
}

/// An L-shaped prism: the concave-edge case a face-push offset gets wrong.
fn l_prism(z0: f64, z1: f64) -> Polyhedron {
    let p = |x: f64, y: f64, z: f64| Point3::new(x, y, z);
    let ring = [
        (0.0, 0.0),
        (4.0, 0.0),
        (4.0, 2.0),
        (2.0, 2.0),
        (2.0, 4.0),
        (0.0, 4.0),
    ];
    let mut faces: Vec<Vec<Point3>> = Vec::new();
    faces.push(ring.iter().rev().map(|&(x, y)| p(x, y, z0)).collect());
    faces.push(ring.iter().map(|&(x, y)| p(x, y, z1)).collect());
    for i in 0..ring.len() {
        let (x0, y0) = ring[i];
        let (x1, y1) = ring[(i + 1) % ring.len()];
        faces.push(vec![
            p(x0, y0, z0),
            p(x1, y1, z0),
            p(x1, y1, z1),
            p(x0, y0, z1),
        ]);
    }
    Polyhedron::new(faces).expect("L-prism is a valid solid")
}

/// Offsetting a non-convex solid handles its reflex edge correctly.
///
/// The L has one concave vertical edge at (2,2). Under an inward offset that
/// corner must move OUTWARD in plan -- away from the material -- which is the
/// opposite of what every convex corner does. A naive face-push leaves a gap
/// there; the miter solve gets it right because the incident planes still
/// meet in one point.
#[test]
fn a_non_convex_solid_offsets_its_reflex_edge_correctly() {
    let l = l_prism(0.0, 2.0);
    let d = 0.5;
    let shrunk = offset_solid(&l, d, OffsetDirection::Inward).expect("offset");

    // The L footprint offset inward by d is again an L, with every edge
    // pulled in by d: outer arms 4 -> 4 - 2d, and the notch corner pushed
    // out so the inner arms are 2 - d wide... giving footprint area
    // (4-2d)^2 - (4-2d-(2-d))^2 for this shape.
    let outer_side = (4.0 - d) - d;
    let notch_side = (4.0 - d) - (2.0 - d);
    let expected = (outer_side * outer_side - notch_side * notch_side) * (2.0 - 2.0 * d);
    let actual = volume(&shrunk);
    assert!(
        (actual - expected).abs() < 1e-9,
        "L-prism offset inward by {d}: expected {expected}, got {actual}"
    );
}

/// A cube with every zero coordinate of its first face written as `-0.0`,
/// as projecting onto plane axes can write it: the same corners, which
/// must weld to their twins on the other faces (#269).
fn box_with_signed_zeros(max: f64) -> Polyhedron {
    let base = box_solid([0.0, 0.0, 0.0], [max, max, max]);
    let mut faces = base.faces().to_vec();
    let flip = |v: f64| if v == 0.0 { -0.0 } else { v };
    for p in &mut faces[0] {
        *p = Point3::new(flip(p.x), flip(p.y), flip(p.z));
    }
    Polyhedron::new(faces).expect("the same cube")
}

#[test]
fn a_corner_written_as_negative_zero_welds_to_its_twin() {
    let cube = box_with_signed_zeros(2.0);
    assert_eq!(
        triangulate(&cube).positions.len(),
        8,
        "a cube has 8 corners"
    );
    assert!((volume(&cube) - 8.0).abs() < 1e-12);
    // Each corner moves once, to the meeting point of all three planes.
    let grown = offset_solid(&cube, 0.5, OffsetDirection::Outward).expect("offset");
    assert!((volume(&grown) - 27.0).abs() < 1e-9, "{}", volume(&grown));
}
