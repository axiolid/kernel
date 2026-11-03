//! Triangles that nearly share a plane, within a linear tolerance (#391).
//!
//! Exact contact and a crossing at an angle are decided cases and are not
//! reported; only a residue-thin offset is.

use axiolid_core::Point3;
use axiolid_inspect::near_coplanar;
use axiolid_mesh::TriMesh;

fn p(x: f64, y: f64, z: f64) -> Point3 {
    Point3::new(x, y, z)
}

/// A right triangle with legs of `1` at height `z`, its right angle at `origin`.
fn triangle(origin: [f64; 2], z: f64) -> TriMesh {
    TriMesh::new(
        vec![
            p(origin[0], origin[1], z),
            p(origin[0] + 1.0, origin[1], z),
            p(origin[0], origin[1] + 1.0, z),
        ],
        vec![0, 1, 2],
    )
}

const EPS: f64 = 1e-6;

#[test]
fn exact_face_contact_is_not_near() {
    assert_eq!(
        near_coplanar(&triangle([0.0, 0.0], 3.0), &triangle([0.25, 0.0], 3.0), EPS),
        None
    );
}

#[test]
fn a_rounding_residue_is_near_with_its_offset() {
    // The same plane at a site far from the origin, one ulp apart.
    let z: f64 = 5.6e6;
    let above = f64::from_bits(z.to_bits() + 1);
    let found = near_coplanar(&triangle([0.0, 0.0], z), &triangle([0.25, 0.0], above), EPS)
        .expect("a residue is near-coplanar");
    assert_eq!((found.first, found.second), (0, 0));
    assert!(found.offset > 0.0 && found.offset <= 1e-15 * z, "{found:?}");
}

#[test]
fn half_the_tolerance_is_near_and_twice_it_is_not() {
    let (a, near, far) = (
        triangle([0.0, 0.0], 0.0),
        triangle([0.25, 0.0], 0.5 * EPS),
        triangle([0.25, 0.0], 2.0 * EPS),
    );
    let found = near_coplanar(&a, &near, EPS).expect("within the tolerance");
    assert!((found.offset - 0.5 * EPS).abs() <= 1e-18, "{found:?}");
    assert_eq!(near_coplanar(&a, &far, EPS), None);
}

#[test]
fn a_crossing_at_an_angle_is_not_near() {
    // Tilted 30 degrees about the x axis through the first triangle.
    let (sin, cos) = (0.5_f64, 0.75_f64.sqrt());
    let tilted = TriMesh::new(
        vec![p(0.1, 0.1, 0.0), p(1.1, 0.1, 0.0), p(0.1, 0.1 + cos, sin)],
        vec![0, 1, 2],
    );
    assert_eq!(
        near_coplanar(&triangle([0.0, 0.0], 0.0), &tilted, EPS),
        None
    );
}

#[test]
fn exactly_coplanar_is_not_near_even_at_zero_linear() {
    let (a, b) = (triangle([0.0, 0.0], 1.0), triangle([0.25, 0.0], 1.0));
    assert_eq!(near_coplanar(&a, &b, 0.0), None);
    assert_eq!(near_coplanar(&a, &b, EPS), None);
}

#[test]
fn a_bad_linear_is_none() {
    let (a, b) = (triangle([0.0, 0.0], 0.0), triangle([0.25, 0.0], 0.5 * EPS));
    for linear in [f64::NAN, -EPS, f64::INFINITY] {
        assert_eq!(near_coplanar(&a, &b, linear), None, "{linear}");
    }
}

#[test]
fn a_degenerate_triangle_is_skipped() {
    // The first triangle of `b` has collinear corners, as far off the plane
    // as the second, which is real.
    let b = TriMesh::new(
        vec![
            p(0.0, 0.0, 0.5 * EPS),
            p(0.5, 0.0, 0.5 * EPS),
            p(1.0, 0.0, 0.5 * EPS),
            p(0.25, 0.0, 0.5 * EPS),
            p(1.25, 0.0, 0.5 * EPS),
            p(0.25, 1.0, 0.5 * EPS),
        ],
        vec![0, 1, 2, 3, 4, 5],
    );
    let found = near_coplanar(&triangle([0.0, 0.0], 0.0), &b, EPS).expect("the real one");
    assert_eq!(found.second, 1);
}

#[test]
fn the_first_pair_in_ascending_order_is_reported() {
    // `a`: a far triangle, then two that have a partner in `b`.
    let a = TriMesh::new(
        vec![
            p(10.0, 0.0, 0.0),
            p(11.0, 0.0, 0.0),
            p(10.0, 1.0, 0.0),
            p(0.0, 5.0, 0.0),
            p(1.0, 5.0, 0.0),
            p(0.0, 6.0, 0.0),
            p(0.0, 0.0, 0.0),
            p(1.0, 0.0, 0.0),
            p(0.0, 1.0, 0.0),
        ],
        vec![0, 1, 2, 3, 4, 5, 6, 7, 8],
    );
    let b = TriMesh::new(
        vec![
            p(0.25, 0.0, 0.5 * EPS),
            p(1.25, 0.0, 0.5 * EPS),
            p(0.25, 1.0, 0.5 * EPS),
            p(0.25, 5.0, 0.5 * EPS),
            p(1.25, 5.0, 0.5 * EPS),
            p(0.25, 6.0, 0.5 * EPS),
        ],
        vec![0, 1, 2, 3, 4, 5],
    );
    let found = near_coplanar(&a, &b, EPS).expect("near");
    assert_eq!((found.first, found.second), (1, 1));
}

#[test]
fn a_malformed_index_is_none() {
    let broken = TriMesh::new(vec![p(0.0, 0.0, 0.0)], vec![0, 1, 2]);
    assert_eq!(
        near_coplanar(&broken, &triangle([0.0, 0.0], 0.0), EPS),
        None
    );
    assert_eq!(
        near_coplanar(&triangle([0.0, 0.0], 0.0), &broken, EPS),
        None
    );
}

#[test]
fn a_tiny_tilted_triangle_on_the_plane_is_not_near() {
    // Its corners lie within the residue band of the big triangle's plane,
    // but its own plane stands at 45 degrees, so the big triangle's corners
    // are far from it.
    let tilted = TriMesh::new(
        vec![
            p(0.25, 0.25, 0.1 * EPS),
            p(0.25 + 1e-7, 0.25, 0.1 * EPS),
            p(0.25 + 1e-7, 0.25, 0.1 * EPS + 1e-7),
        ],
        vec![0, 1, 2],
    );
    let big = triangle([0.0, 0.0], 0.0);
    assert_eq!(near_coplanar(&big, &tilted, EPS), None);
    assert_eq!(near_coplanar(&tilted, &big, EPS), None);
}

#[test]
fn triangles_far_apart_in_their_plane_are_not_near() {
    // The boxes of the two overlap; the triangles are far apart.
    let a = TriMesh::new(
        vec![p(0.0, 0.0, 0.0), p(10.0, 0.0, 0.0), p(0.0, 10.0, 0.0)],
        vec![0, 1, 2],
    );
    let b = TriMesh::new(
        vec![
            p(9.0, 9.0, 0.5 * EPS),
            p(10.0, 9.0, 0.5 * EPS),
            p(9.0, 10.0, 0.5 * EPS),
        ],
        vec![0, 1, 2],
    );
    assert_eq!(near_coplanar(&a, &b, EPS), None);
    assert_eq!(near_coplanar(&b, &a, EPS), None);
}

#[test]
fn the_lowest_second_triangle_is_reported() {
    // One big triangle over two small ones: the lower index wins.
    let a = triangle([0.0, 0.0], 0.0);
    let small = |x: f64| {
        [
            p(x, 0.1, 0.5 * EPS),
            p(x + 0.1, 0.1, 0.5 * EPS),
            p(x, 0.2, 0.5 * EPS),
        ]
    };
    let b = TriMesh::new(
        [small(0.1), small(0.5)].into_iter().flatten().collect(),
        vec![0, 1, 2, 3, 4, 5],
    );
    let found = near_coplanar(&a, &b, EPS).expect("near");
    assert_eq!((found.first, found.second), (0, 0));
}
