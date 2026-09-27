//! Barycentric and mean-value coordinates (#143, ledger row H5).
//!
//! Oracles are closed forms: in the right triangle `(0,0), (4,0), (0,2)` the
//! point `(x, y)` has weights `(1 - x/4 - y/2, x/4, y/2)`; in the corner
//! tetrahedron `(1 - x - y - z, x, y, z)`; and on a triangle, mean-value
//! coordinates are the barycentric ones.

use axiolid_core::{Point2, Point3, Polygon2, Tolerance, Triangle2, Triangle3, Vec3};
use axiolid_spatial::{
    mean_value_coordinates2, tetrahedron_barycentric, triangle_barycentric2, triangle_barycentric3,
    BarycentricError,
};

const TOL: Tolerance = Tolerance::METRE;

fn p2(x: f64, y: f64) -> Point2 {
    Point2::new(x, y)
}

fn close(got: &[f64], want: &[f64], eps: f64) {
    assert_eq!(got.len(), want.len());
    for (g, w) in got.iter().zip(want) {
        assert!((g - w).abs() <= eps, "{got:?} vs {want:?}");
    }
}

/// The weights sum to one and reproduce `point` from `corners`.
fn reproduces2(weights: &[f64], corners: &[Point2], point: Point2) {
    let sum: f64 = weights.iter().sum();
    assert!((sum - 1.0).abs() < 1e-12, "sum {sum}");
    let back = corners
        .iter()
        .zip(weights)
        .fold(Point2::ZERO, |acc, (c, w)| acc + *c * *w);
    assert!((back - point).length() < 1e-12, "{back} vs {point}");
}

fn right_triangle() -> Triangle2 {
    Triangle2::new(p2(0.0, 0.0), p2(4.0, 0.0), p2(0.0, 2.0))
}

fn closed_form(x: f64, y: f64) -> [f64; 3] {
    [1.0 - x / 4.0 - y / 2.0, x / 4.0, y / 2.0]
}

#[test]
fn a_2d_triangle_has_the_closed_form_weights_inside_and_out() {
    for (x, y) in [(1.0, 0.5), (0.3, 1.1), (5.0, 1.0), (-2.0, 3.0), (2.0, 0.0)] {
        let got = triangle_barycentric2(&right_triangle(), p2(x, y), TOL).unwrap();
        close(&got, &closed_form(x, y), 1e-15);
    }
    // Winding does not matter: swap b and c, and their weights swap.
    let t = right_triangle();
    let cw = Triangle2::new(t.a, t.c, t.b);
    let got = triangle_barycentric2(&cw, p2(1.0, 0.5), TOL).unwrap();
    close(&got, &[0.5, 0.25, 0.25], 1e-15);
}

#[test]
fn a_triangle_corner_gets_exactly_one() {
    let t = right_triangle();
    for (k, corner) in [t.a, t.b, t.c].into_iter().enumerate() {
        let got = triangle_barycentric2(&t, corner, TOL).unwrap();
        let mut want = [0.0; 3];
        want[k] = 1.0;
        assert_eq!(got, want);
    }
}

#[test]
fn thin_or_non_finite_triangles_are_refused() {
    let thin = Triangle2::new(p2(0.0, 0.0), p2(1.0, 0.0), p2(2.0, 1e-7));
    match triangle_barycentric2(&thin, p2(0.5, 0.0), TOL) {
        Err(BarycentricError::Degenerate { thickness }) => assert!(thickness < 1e-6),
        other => panic!("{other:?}"),
    }
    // Above the tolerance it is answered.
    let slim = Triangle2::new(p2(0.0, 0.0), p2(1.0, 0.0), p2(2.0, 1e-5));
    assert!(triangle_barycentric2(&slim, p2(0.5, 0.0), TOL).is_ok());
    let nan = triangle_barycentric2(&right_triangle(), p2(f64::NAN, 0.0), TOL);
    assert_eq!(nan, Err(BarycentricError::NonFinite));
}

#[test]
fn a_3d_triangle_answers_for_the_projection() {
    // The right triangle, tilted and moved: weights are affine-invariant.
    let lift = |p: Point2| Point3::new(3.0 + p.x, -1.0 + 0.6 * p.y, 2.0 + 0.8 * p.y);
    let t = right_triangle();
    let t3 = Triangle3::new(lift(t.a), lift(t.b), lift(t.c));
    let normal = (t3.b - t3.a).cross(t3.c - t3.a).normalize();
    for (x, y) in [(1.0, 0.5), (0.3, 1.1), (5.0, 1.0)] {
        let on = lift(p2(x, y));
        close(
            &triangle_barycentric3(&t3, on, TOL).unwrap(),
            &closed_form(x, y),
            1e-14,
        );
        // Off the plane, the same weights: they locate the projection.
        let off = on + normal * 7.5;
        close(
            &triangle_barycentric3(&t3, off, TOL).unwrap(),
            &closed_form(x, y),
            1e-14,
        );
    }
    for (k, corner) in [t3.a, t3.b, t3.c].into_iter().enumerate() {
        let got = triangle_barycentric3(&t3, corner, TOL).unwrap();
        let mut want = [0.0; 3];
        want[k] = 1.0;
        assert_eq!(got, want);
    }
    let flat = Triangle3::new(Vec3::ZERO, Vec3::X, Vec3::X * 2.0 + Vec3::Y * 1e-8);
    assert!(matches!(
        triangle_barycentric3(&flat, Vec3::ZERO, TOL),
        Err(BarycentricError::Degenerate { .. })
    ));
}

#[test]
fn a_tetrahedron_has_the_closed_form_weights() {
    let corners = [Vec3::ZERO, Vec3::X, Vec3::Y, Vec3::Z];
    for [x, y, z] in [
        [0.25, 0.25, 0.25],
        [0.1, 0.2, 0.3],
        [1.0, 1.0, 1.0],
        [-0.5, 0.25, 2.0],
    ] {
        let got = tetrahedron_barycentric(corners, Point3::new(x, y, z), TOL).unwrap();
        close(&got, &[1.0 - x - y - z, x, y, z], 1e-15);
    }
    // Left-handed corner order: same weights, permuted with the corners.
    let swapped = [Vec3::ZERO, Vec3::Y, Vec3::X, Vec3::Z];
    let got = tetrahedron_barycentric(swapped, Point3::new(0.1, 0.2, 0.3), TOL).unwrap();
    close(&got, &[0.4, 0.2, 0.1, 0.3], 1e-15);
    for (k, corner) in corners.into_iter().enumerate() {
        let got = tetrahedron_barycentric(corners, corner, TOL).unwrap();
        let mut want = [0.0; 4];
        want[k] = 1.0;
        assert_eq!(got, want);
    }
    // A tetrahedron flattened to within the tolerance is refused.
    let flat = [Vec3::ZERO, Vec3::X, Vec3::Y, Vec3::new(0.3, 0.3, 1e-7)];
    assert!(matches!(
        tetrahedron_barycentric(flat, Vec3::ZERO, TOL),
        Err(BarycentricError::Degenerate { .. })
    ));
}

#[test]
fn mean_value_coordinates_on_a_triangle_are_barycentric() {
    let t = right_triangle();
    let polygon = Polygon2::new(vec![t.a, t.b, t.c]);
    for (x, y) in [
        (1.0, 0.5),
        (0.3, 1.1),
        (5.0, 1.0),
        (-2.0, 3.0),
        (0.01, 0.01),
    ] {
        let got = mean_value_coordinates2(&polygon, p2(x, y), TOL).unwrap();
        close(&got, &closed_form(x, y), 1e-13);
    }
}

fn l_shape() -> Vec<Point2> {
    vec![
        p2(0.0, 0.0),
        p2(3.0, 0.0),
        p2(3.0, 1.0),
        p2(1.0, 1.0),
        p2(1.0, 3.0),
        p2(0.0, 3.0),
    ]
}

#[test]
fn a_non_convex_polygon_reproduces_its_points() {
    let vertices = l_shape();
    let polygon = Polygon2::new(vertices.clone());
    let reversed = Polygon2::new(vertices.iter().rev().copied().collect());
    // Points in both arms, the reflex corner's neighbourhood, and the
    // corner square.
    for point in [
        p2(2.5, 0.5),
        p2(0.5, 2.5),
        p2(0.5, 0.5),
        p2(1.1, 0.9),
        p2(0.9, 1.1),
        p2(0.99, 0.99),
    ] {
        let got = mean_value_coordinates2(&polygon, point, TOL).unwrap();
        reproduces2(&got, &vertices, point);
        // Clockwise, the same weights in reverse order.
        let mut back = mean_value_coordinates2(&reversed, point, TOL).unwrap();
        back.reverse();
        close(&got, &back, 1e-14);
    }
    // In the polygon's kernel, the corner square, every weight is positive.
    for point in [p2(0.5, 0.5), p2(0.2, 0.9), p2(0.9, 0.2)] {
        let got = mean_value_coordinates2(&polygon, point, TOL).unwrap();
        assert!(got.iter().all(|&w| w > 0.0), "{got:?}");
    }
}

#[test]
fn a_square_centre_weighs_every_corner_equally() {
    let square = Polygon2::new(vec![p2(0.0, 0.0), p2(2.0, 0.0), p2(2.0, 2.0), p2(0.0, 2.0)]);
    let got = mean_value_coordinates2(&square, p2(1.0, 1.0), TOL).unwrap();
    close(&got, &[0.25; 4], 1e-15);
    // Outside a convex polygon the weights are still defined.
    let far = p2(5.0, -1.0);
    let got = mean_value_coordinates2(&square, far, TOL).unwrap();
    reproduces2(&got, &square.vertices, far);
}

#[test]
fn the_boundary_is_interpolated_linearly() {
    let vertices = l_shape();
    let polygon = Polygon2::new(vertices.clone());
    for (k, &corner) in vertices.iter().enumerate() {
        let got = mean_value_coordinates2(&polygon, corner, TOL).unwrap();
        let mut want = vec![0.0; vertices.len()];
        want[k] = 1.0;
        assert_eq!(got, want);
    }
    // Within the tolerance of a corner but past both its edges' ends: the
    // corner alone, not the formula's near-miss.
    let got = mean_value_coordinates2(&polygon, p2(-3e-7, -4e-7), TOL).unwrap();
    assert_eq!(got, vec![1.0, 0.0, 0.0, 0.0, 0.0, 0.0]);
    // A quarter of the way along edge 3, from (1, 1) to (1, 3).
    let got = mean_value_coordinates2(&polygon, p2(1.0, 1.5), TOL).unwrap();
    assert_eq!(got, vec![0.0, 0.0, 0.0, 0.75, 0.25, 0.0]);
    // Just inside that edge, the weights approach it continuously.
    let near = mean_value_coordinates2(&polygon, p2(0.999, 1.5), TOL).unwrap();
    close(&near, &[0.0, 0.0, 0.0, 0.75, 0.25, 0.0], 0.02);
    reproduces2(&near, &vertices, p2(0.999, 1.5));
}

#[test]
fn polygons_that_bound_nothing_are_refused() {
    let at = |pts: &[(f64, f64)]| Polygon2::new(pts.iter().map(|&(x, y)| p2(x, y)).collect());
    let origin = p2(0.25, 0.25);
    assert_eq!(
        mean_value_coordinates2(&at(&[(0.0, 0.0), (1.0, 0.0)]), origin, TOL),
        Err(BarycentricError::TooFewVertices { count: 2 })
    );
    assert_eq!(
        mean_value_coordinates2(
            &at(&[(0.0, 0.0), (1.0, 0.0), (1.0, 0.0), (0.0, 1.0)]),
            origin,
            TOL
        ),
        Err(BarycentricError::ShortEdge { index: 1 })
    );
    let bowtie = at(&[(0.0, 0.0), (1.0, 1.0), (1.0, 0.0), (0.0, 1.0)]);
    assert!(matches!(
        mean_value_coordinates2(&bowtie, origin, TOL),
        Err(BarycentricError::SelfIntersecting { .. })
    ));
    // A spike: edge 2 runs back down edge 1.
    let spike = at(&[(0.0, 0.0), (2.0, 0.0), (2.0, 2.0), (2.0, 1.0), (0.0, 1.0)]);
    assert!(matches!(
        mean_value_coordinates2(&spike, origin, TOL),
        Err(BarycentricError::SelfIntersecting { .. })
    ));
    // A vertex touching a non-adjacent edge.
    let touching = at(&[(0.0, 0.0), (4.0, 0.0), (4.0, 2.0), (2.0, 0.0), (0.0, 2.0)]);
    assert!(matches!(
        mean_value_coordinates2(&touching, origin, TOL),
        Err(BarycentricError::SelfIntersecting { .. })
    ));
    // A sliver thinner than the tolerance: its far side comes within the
    // tolerance of its near side, which the simplicity check sees first.
    let sliver = at(&[(0.0, 0.0), (1.0, 0.0), (2.0, 1e-7)]);
    assert!(matches!(
        mean_value_coordinates2(&sliver, origin, TOL),
        Err(BarycentricError::SelfIntersecting { .. } | BarycentricError::Degenerate { .. })
    ));
    let nan = at(&[(0.0, 0.0), (1.0, 0.0), (f64::NAN, 1.0)]);
    assert_eq!(
        mean_value_coordinates2(&nan, origin, TOL),
        Err(BarycentricError::NonFinite)
    );
}
