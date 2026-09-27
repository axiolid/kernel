//! Collinear but non-intersecting edges are valid input (#188).

use axiolid_core::{Point2, Tolerance};
use axiolid_overlay::{Polygon, Region, Ring};

fn tol() -> Tolerance {
    Tolerance::METRE
}

fn ring(points: &[(f64, f64)]) -> Ring {
    Ring {
        points: points.iter().map(|&(x, y)| Point2::new(x, y)).collect(),
    }
}

fn square(x0: f64, y0: f64, x1: f64, y1: f64) -> Region {
    Region::new(
        vec![Polygon {
            outer: ring(&[(x0, y0), (x1, y0), (x1, y1), (x0, y1)]),
            holes: vec![],
        }],
        tol(),
    )
    .unwrap()
}

#[test]
fn a_vertex_on_the_extension_of_a_non_adjacent_edge_is_valid() {
    // A convex pentagon with a straight vertex B: the vertex A, an end of
    // the edge E-A, lies on the extension of the edge B-C, which does not
    // touch E-A. Collinear, not intersecting.
    let pentagon = ring(&[(0.0, 0.0), (2.0, 0.0), (4.0, 0.0), (4.0, 3.0), (0.0, 3.0)]);
    let region = Region::new(
        vec![Polygon {
            outer: pentagon,
            holes: vec![],
        }],
        tol(),
    )
    .expect("a vertex on another edge's extension is not an intersection");
    assert!((region.area() - 12.0).abs() < 1e-12);
    // And it takes part in booleans.
    let union = region.union(&square(3.0, -1.0, 5.0, 1.0), tol()).unwrap();
    assert!(
        (union.area() - (12.0 + 4.0 - 1.0)).abs() < 1e-9,
        "{}",
        union.area()
    );
}

#[test]
fn overlay_output_with_collinear_edges_feeds_back_into_booleans() {
    // A room and a doorway strip on one wall: their union has collinear
    // consecutive edges along the wall on either side of the strip.
    let room = square(0.0, 0.0, 4.0, 3.0);
    let doorway = square(1.5, 3.0, 2.5, 3.5);
    let union = room.union(&doorway, tol()).unwrap();
    assert_eq!(union.component_count(), 1);
    let expected = 12.0 + 0.5;
    assert!((union.area() - expected).abs() < 1e-9, "{}", union.area());
    // Fed back: a difference, a union and an intersection with it.
    let column = square(0.5, 0.5, 1.0, 1.0);
    let difference = union.difference(&column, tol()).unwrap();
    assert!((difference.area() - (expected - 0.25)).abs() < 1e-9);
    let rebuilt =
        Region::new(difference.polygons().to_vec(), tol()).expect("overlay output is valid input");
    let again = rebuilt.union(&square(3.0, 3.0, 3.5, 3.2), tol()).unwrap();
    assert!((again.area() - (expected - 0.25 + 0.1)).abs() < 1e-9);
    let within = again.intersection(&room, tol()).unwrap();
    assert!((within.area() - (12.0 - 0.25)).abs() < 1e-9);
}
