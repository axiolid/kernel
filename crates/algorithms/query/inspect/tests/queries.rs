//! Mesh queries: clearance, containment, ray casting, genus (#89).
//!
//! The oracles are hand-computed: two cubes a known distance apart, a ray
//! with a hand-derived parameter, and the L-prism notch from #77 -- the
//! non-convex point where a face-plane test gives the wrong answer.

use axiolid_core::{Point3, Vec3};
use axiolid_inspect::{contains, genus, min_gap, ray_cast, winding_number, GenusError};
use axiolid_mesh::TriMesh;

/// An axis-aligned box as a closed, outward-wound triangle mesh.
fn box_mesh(min: [f64; 3], max: [f64; 3]) -> TriMesh {
    let [x0, y0, z0] = min;
    let [x1, y1, z1] = max;
    let p = |x: f64, y: f64, z: f64| Point3::new(x, y, z);
    let positions = vec![
        p(x0, y0, z0),
        p(x1, y0, z0),
        p(x1, y1, z0),
        p(x0, y1, z0),
        p(x0, y0, z1),
        p(x1, y0, z1),
        p(x1, y1, z1),
        p(x0, y1, z1),
    ];
    let indices = vec![
        0, 2, 1, 0, 3, 2, // bottom, wound to face -z
        4, 5, 6, 4, 6, 7, // top
        0, 1, 5, 0, 5, 4, // -y
        1, 2, 6, 1, 6, 5, // +x
        2, 3, 7, 2, 7, 6, // +y
        3, 0, 4, 3, 4, 7, // -x
    ];
    TriMesh::new(positions, indices)
}

#[test]
fn two_cubes_half_a_unit_apart_report_that_gap() {
    let a = box_mesh([0.0, 0.0, 0.0], [1.0, 1.0, 1.0]);
    let b = box_mesh([1.5, 0.0, 0.0], [2.5, 1.0, 1.0]);

    let gap = min_gap(&a, &b, 2.0).expect("within the search length");
    assert!(
        (gap - 0.5).abs() < 1e-12,
        "expected a 0.5 gap between the facing walls, got {gap}"
    );
}

#[test]
fn overlapping_cubes_report_a_zero_gap() {
    // The clash case: surfaces interpenetrate, so the clearance is gone.
    let a = box_mesh([0.0, 0.0, 0.0], [1.0, 1.0, 1.0]);
    let b = box_mesh([0.5, 0.5, 0.5], [1.5, 1.5, 1.5]);

    assert_eq!(
        min_gap(&a, &b, 2.0),
        Some(0.0),
        "overlap must read as clash"
    );
}

#[test]
fn solids_beyond_the_search_length_report_nothing() {
    // "Far apart" must be distinguishable from "touching": a caller doing
    // clash detection treats these as opposite outcomes.
    let a = box_mesh([0.0, 0.0, 0.0], [1.0, 1.0, 1.0]);
    let b = box_mesh([50.0, 0.0, 0.0], [51.0, 1.0, 1.0]);

    assert_eq!(
        min_gap(&a, &b, 1.0),
        None,
        "nothing within the search radius"
    );
}

/// An L-shaped prism as a closed triangle mesh.
///
/// The notch is the point of this fixture: a point there is on the inner
/// side of every face plane while being outside the solid, so a convex
/// half-space test gets it wrong and ray parity gets it right.
fn l_prism() -> TriMesh {
    let p = |x: f64, y: f64, z: f64| Point3::new(x, y, z);
    // Footprint, counter-clockwise: (0,0) (4,0) (4,2) (2,2) (2,4) (0,4).
    let ring = [
        (0.0, 0.0),
        (4.0, 0.0),
        (4.0, 2.0),
        (2.0, 2.0),
        (2.0, 4.0),
        (0.0, 4.0),
    ];
    let mut positions = Vec::new();
    for (x, y) in ring {
        positions.push(p(x, y, 0.0));
    }
    for (x, y) in ring {
        positions.push(p(x, y, 1.0));
    }

    let n = ring.len() as u32;
    let mut indices = Vec::new();
    // Caps: fan from vertex 0 is safe here because the L is star-shaped
    // about its first corner.
    for i in 1..n - 1 {
        indices.extend_from_slice(&[0, i + 1, i]);
        indices.extend_from_slice(&[n, n + i, n + i + 1]);
    }
    // Walls.
    for i in 0..n {
        let j = (i + 1) % n;
        indices.extend_from_slice(&[i, j, n + j]);
        indices.extend_from_slice(&[i, n + j, n + i]);
    }
    TriMesh::new(positions, indices)
}

#[test]
fn a_point_in_the_l_notch_is_outside() {
    // (3, 3) is inside the bounding box and on the inner side of every face
    // plane, but outside the L. This is the case #89 names: a convex
    // half-space test answers "inside" here and is wrong.
    let mesh = l_prism();
    let notch = Point3::new(3.0, 3.0, 0.5);

    assert_eq!(winding_number(&mesh, notch), Some(0), "notch winds zero");
    assert_eq!(contains(&mesh, notch), Some(false), "notch is outside");
}

#[test]
fn points_in_both_arms_are_inside() {
    let mesh = l_prism();
    for point in [Point3::new(3.0, 1.0, 0.5), Point3::new(1.0, 3.0, 0.5)] {
        assert_eq!(
            winding_number(&mesh, point),
            Some(1),
            "{point:?} lies in an arm of the L"
        );
        assert_eq!(contains(&mesh, point), Some(true), "{point:?} is inside");
    }
}

#[test]
fn a_point_outside_the_bounding_box_is_outside() {
    let mesh = box_mesh([0.0, 0.0, 0.0], [1.0, 1.0, 1.0]);
    let far = Point3::new(9.0, 9.0, 9.0);

    assert_eq!(winding_number(&mesh, far), Some(0));
    assert_eq!(contains(&mesh, far), Some(false));
}

#[test]
fn a_ray_returns_the_face_it_strikes_and_where() {
    // Fire from outside along -x at the cube's +x wall (x = 1). From x = 3
    // the hit is at t = 2 exactly, hand-computed and independent of the
    // implementation.
    let cube = box_mesh([0.0, 0.0, 0.0], [1.0, 1.0, 1.0]);
    let origin = Point3::new(3.0, 0.4, 0.6);
    let hit = ray_cast(&cube, origin, Vec3::new(-1.0, 0.0, 0.0)).expect("the ray strikes the cube");

    assert!(
        (hit.t - 2.0).abs() < 1e-12,
        "expected t = 2.0 at the x = 1 wall, got {}",
        hit.t
    );
    assert!(
        (hit.point.x - 1.0).abs() < 1e-12,
        "hit point should lie on the x = 1 plane, got {:?}",
        hit.point
    );
    // The reported triangle must actually be one of the +x wall's two.
    let corners = &cube.indices[hit.triangle * 3..hit.triangle * 3 + 3];
    assert!(
        corners
            .iter()
            .all(|c| (cube.positions[*c as usize].x - 1.0).abs() < 1e-12),
        "the named triangle should lie in the x = 1 plane"
    );
}

#[test]
fn a_ray_that_misses_returns_nothing() {
    let cube = box_mesh([0.0, 0.0, 0.0], [1.0, 1.0, 1.0]);
    let origin = Point3::new(3.0, 9.0, 9.0);
    assert!(ray_cast(&cube, origin, Vec3::new(-1.0, 0.0, 0.0)).is_none());
}

#[test]
fn the_nearest_hit_is_returned_not_the_first_found() {
    // A ray crossing the whole cube meets two walls; the near one wins.
    let cube = box_mesh([0.0, 0.0, 0.0], [1.0, 1.0, 1.0]);
    let hit = ray_cast(&cube, Point3::new(3.0, 0.4, 0.6), Vec3::new(-1.0, 0.0, 0.0)).expect("hit");
    assert!(
        (hit.t - 2.0).abs() < 1e-12,
        "near wall at t = 2, far wall at t = 3; got {}",
        hit.t
    );
}

/// A torus as a closed quad-grid triangle mesh: genus 1 by construction.
fn torus(major_segments: u32, minor_segments: u32) -> TriMesh {
    let (major_radius, minor_radius) = (3.0, 1.0);
    let mut positions = Vec::new();
    for i in 0..major_segments {
        let theta = std::f64::consts::TAU * f64::from(i) / f64::from(major_segments);
        for j in 0..minor_segments {
            let phi = std::f64::consts::TAU * f64::from(j) / f64::from(minor_segments);
            let radius = major_radius + minor_radius * phi.cos();
            positions.push(Point3::new(
                radius * theta.cos(),
                radius * theta.sin(),
                minor_radius * phi.sin(),
            ));
        }
    }
    let mut indices = Vec::new();
    for i in 0..major_segments {
        for j in 0..minor_segments {
            let next_i = (i + 1) % major_segments;
            let next_j = (j + 1) % minor_segments;
            let a = i * minor_segments + j;
            let b = next_i * minor_segments + j;
            let c = next_i * minor_segments + next_j;
            let d = i * minor_segments + next_j;
            indices.extend_from_slice(&[a, b, c]);
            indices.extend_from_slice(&[a, c, d]);
        }
    }
    TriMesh::new(positions, indices)
}

#[test]
fn a_cube_has_genus_zero() {
    let cube = box_mesh([0.0, 0.0, 0.0], [1.0, 1.0, 1.0]);
    assert_eq!(genus(&cube), Ok(0));
}

#[test]
fn a_torus_has_genus_one() {
    // The case that distinguishes genus from "is it closed": a torus is a
    // perfectly good closed manifold with a hole through it.
    assert_eq!(genus(&torus(16, 8)), Ok(1));
}

#[test]
fn an_open_mesh_refuses_rather_than_guessing() {
    // A single triangle has three boundary edges. Euler's formula would
    // still produce a number; it would just be meaningless.
    let sheet = TriMesh::new(
        vec![
            Point3::new(0.0, 0.0, 0.0),
            Point3::new(1.0, 0.0, 0.0),
            Point3::new(0.0, 1.0, 0.0),
        ],
        vec![0, 1, 2],
    );
    assert!(matches!(
        genus(&sheet),
        Err(GenusError::NotClosedManifold { boundary: 3, .. })
    ));
}

/// Two disjoint cubes: the simplest multi-component case (#98).
///
/// `chi = 2c - 2g` gives `chi = 4` here, so the single-component formula
/// `g = (2 - chi) / 2` yields `-1`. Converting that to `u32` used to clamp
/// to `Ok(0)`, which is indistinguishable from a genuine sphere.
#[test]
fn two_disjoint_cubes_refuse_rather_than_reporting_genus_zero() {
    let mut mesh = box_mesh([0.0, 0.0, 0.0], [1.0, 1.0, 1.0]);
    let far = box_mesh([4.0, 0.0, 0.0], [5.0, 1.0, 1.0]);
    let base = mesh.positions.len() as u32;
    mesh.positions.extend(far.positions.iter().copied());
    mesh.indices.extend(far.indices.iter().map(|i| i + base));

    match genus(&mesh) {
        Err(GenusError::MultipleComponents {
            components,
            characteristic,
        }) => {
            assert_eq!(components, 2, "two cubes are two components");
            assert_eq!(characteristic, 4, "chi = 2c for two genus-0 solids");
        }
        other => panic!("expected MultipleComponents, got {other:?}"),
    }
}

/// A solid with interior cavities is the case that motivated #98.
///
/// A cube containing one hollow cube is 2 components with `chi = 4`. The
/// old clamp reported `Ok(0)` -- a hollow solid and a plain sphere gave the
/// same answer, so no caller could tell them apart.
#[test]
fn a_solid_with_a_cavity_is_not_reported_as_a_sphere() {
    let mut mesh = box_mesh([0.0, 0.0, 0.0], [4.0, 4.0, 4.0]);
    let cavity = box_mesh([1.0, 1.0, 1.0], [2.0, 2.0, 2.0]);
    let base = mesh.positions.len() as u32;
    mesh.positions.extend(cavity.positions.iter().copied());
    // Reversed winding: an interior void faces inward.
    for t in cavity.indices.chunks(3) {
        mesh.indices
            .extend_from_slice(&[t[0] + base, t[2] + base, t[1] + base]);
    }

    assert!(
        matches!(
            genus(&mesh),
            Err(GenusError::MultipleComponents { components: 2, .. })
        ),
        "a cavity-bearing solid must refuse, not report genus 0"
    );
}

/// The single-component path must keep working unchanged.
#[test]
fn one_component_still_reports_its_genus() {
    assert_eq!(genus(&box_mesh([0.0, 0.0, 0.0], [1.0, 1.0, 1.0])), Ok(0));
    assert_eq!(genus(&torus(16, 8)), Ok(1));
}

#[test]
fn containment_holds_at_scales_where_a_tolerance_test_fails() {
    // A tolerance-based orientation test treats anything within epsilon of a
    // face plane as "on" it. On a mesh whose features are SMALLER than that
    // epsilon, every crossing collapses to a tie and the parity is lost.
    //
    // This cube is 1e-7 across, so a 1e-9 absolute tolerance on the
    // orient3d determinant -- whose magnitude scales with the CUBE of the
    // edge length, here 1e-21 -- cannot resolve any crossing at all. The
    // exact predicate is scale-free and still answers correctly.
    let small = box_mesh([0.0, 0.0, 0.0], [1e-7, 1e-7, 1e-7]);
    let inside = Point3::new(5e-8, 5e-8, 5e-8);
    let outside = Point3::new(5e-7, 5e-8, 5e-8);

    assert_eq!(
        winding_number(&small, inside),
        Some(1),
        "a point at the centre of a small cube is inside it"
    );
    assert_eq!(
        winding_number(&small, outside),
        Some(0),
        "a point five diameters away is outside"
    );
    assert_eq!(contains(&small, inside), Some(true));
    assert_eq!(contains(&small, outside), Some(false));
}
