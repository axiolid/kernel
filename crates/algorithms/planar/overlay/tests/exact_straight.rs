//! Straight-edge booleans are exact (#173).
//!
//! The polygon path used to map its operands onto an `i32` grid scaled to
//! their bounding box, so every output coordinate came back snapped to a
//! step of about 1.5e-8 of the extent, even input vertices the operation
//! does not move. Now an untouched input vertex comes back bit-identical,
//! and a crossing is the double nearest to the exact crossing point.

use axiolid_core::{Frame2, Point2, Tolerance, Vec2};
use axiolid_exact::{Arith, Dyadic};
use axiolid_overlay::{
    overlay, polygon_area, union_soup, FillRule, OverlayInput, OverlayOperation, Polygon, Region,
    Ring,
};

fn tol() -> Tolerance {
    Tolerance::new(1e-9, 1e-9).unwrap()
}

fn frame() -> Frame2 {
    Frame2 {
        origin: Point2::new(0.0, 0.0),
        x: Vec2::new(1.0, 0.0),
        y: Vec2::new(0.0, 1.0),
    }
}

fn rect(x0: f64, y0: f64, x1: f64, y1: f64) -> Polygon {
    Polygon {
        outer: Ring {
            points: vec![
                Point2::new(x0, y0),
                Point2::new(x1, y0),
                Point2::new(x1, y1),
                Point2::new(x0, y1),
            ],
        },
        holes: vec![],
    }
}

fn boolean(a: Polygon, b: Polygon, operation: OverlayOperation) -> Vec<Polygon> {
    let a = OverlayInput {
        frame: frame(),
        polygons: vec![a],
    };
    let b = OverlayInput {
        frame: frame(),
        polygons: vec![b],
    };
    overlay(&a, &b, operation, FillRule::NonZero, tol())
        .expect("valid operands")
        .polygons
}

/// The same points as a set: rings may start anywhere.
fn same_points(ring: &Ring, expected: &Polygon) {
    let mut got: Vec<(u64, u64)> = ring
        .points
        .iter()
        .map(|p| (p.x.to_bits(), p.y.to_bits()))
        .collect();
    let mut want: Vec<(u64, u64)> = expected
        .outer
        .points
        .iter()
        .map(|p| (p.x.to_bits(), p.y.to_bits()))
        .collect();
    got.sort_unstable();
    want.sort_unstable();
    assert_eq!(got, want, "{:?}", ring.points);
}

#[test]
fn a_rectangle_inside_the_clip_comes_back_bit_identical() {
    // The issue's cases: [0,4]x[0,0.2] came back with y = 0.20000000298...
    for wall in [
        rect(0.0, 0.0, 4.0, 0.2),
        rect(0.1, 0.7, 4.3, 0.9),
        rect(0.0, 0.0, 4.0, 0.25),
    ] {
        let result = boolean(
            wall.clone(),
            rect(-1.0, -1.0, 5.0, 5.0),
            OverlayOperation::Intersection,
        );
        assert_eq!(result.len(), 1);
        same_points(&result[0].outer, &wall);
        assert_eq!(polygon_area(&result[0]), polygon_area(&wall));
    }
    let result = boolean(
        rect(0.0, 0.0, 4.0, 0.2),
        rect(-1.0, -1.0, 5.0, 5.0),
        OverlayOperation::Intersection,
    );
    assert_eq!(polygon_area(&result[0]), 0.8);
}

#[test]
fn georeferenced_input_keeps_its_vertices() {
    let (x, y) = (500_000.25, 5_600_000.5);
    let wall = rect(x, y + 0.2, x + 4.0, y + 0.4);
    let site = rect(x - 1.0, y - 1.0, x + 5.0, y + 5.0);
    let result = boolean(wall.clone(), site.clone(), OverlayOperation::Intersection);
    same_points(&result[0].outer, &wall);
    assert_eq!(polygon_area(&result[0]), polygon_area(&wall));
    // Through Region, too, and a difference leaves the site's own corners.
    let site_region = Region::new(vec![site.clone()], tol()).unwrap();
    let wall_region = Region::new(vec![wall.clone()], tol()).unwrap();
    let inside = site_region.intersection(&wall_region, tol()).unwrap();
    same_points(&inside.polygons()[0].outer, &wall);
    let rest = site_region.difference(&wall_region, tol()).unwrap();
    same_points(&rest.polygons()[0].outer, &site);
    same_points(&rest.polygons()[0].holes[0], &wall);
}

/// Whether `value` is `num / den` (with `den > 0`) rounded to nearest:
/// the exact value lies between the midpoints to both neighbours.
fn correctly_rounded(value: f64, num: &Dyadic, den: &Dyadic) -> bool {
    use axiolid_guarantees::Sign;
    // Sign of `num / den - (a + b) / 2`; the midpoint of two neighbouring
    // doubles is no double, but it is a dyadic.
    let beyond = |a: f64, b: f64| {
        let mid = Dyadic::from_f64(a)
            .add(&Dyadic::from_f64(b))
            .mul(&Dyadic::from_f64(0.5));
        num.sub(&mid.mul(den))
            .sign()
            .expect("dyadic signs are decided")
    };
    beyond(value.next_down(), value) != Sign::Negative
        && beyond(value, value.next_up()) != Sign::Positive
}

#[test]
fn a_crossing_is_the_nearest_double_to_the_exact_point() {
    // A sliver from (-1, 0.3) to (2, 0.7) and back along (2, 0.9) cuts the
    // unit square's sides x = 0 and x = 1 at points no double holds.
    let (a, b) = (Point2::new(-1.0, 0.3), Point2::new(2.0, 0.7));
    let sliver = Polygon {
        outer: Ring {
            points: vec![a, b, Point2::new(2.0, 0.9), Point2::new(-1.0, 0.9)],
        },
        holes: vec![],
    };
    let result = boolean(
        rect(0.0, 0.0, 1.0, 1.0),
        sliver,
        OverlayOperation::Intersection,
    );
    assert_eq!(result.len(), 1);
    let d = |v: f64| Dyadic::from_f64(v);
    let mut checked = 0;
    for p in &result[0].outer.points {
        if p.x != 0.0 && p.x != 1.0 {
            continue;
        }
        if p.y == 0.9 || p.y == 1.0 {
            continue;
        }
        // On the line a -> b at x: y = a.y + (x - a.x) (b.y - a.y) / (b.x - a.x).
        let den = d(b.x).sub(&d(a.x));
        let num = d(a.y)
            .mul(&den)
            .add(&d(p.x).sub(&d(a.x)).mul(&d(b.y).sub(&d(a.y))));
        assert!(correctly_rounded(p.y, &num, &den), "{p:?}");
        checked += 1;
    }
    assert_eq!(checked, 2, "{:?}", result[0].outer.points);
}

#[test]
fn union_soup_keeps_shared_vertices() {
    // Two triangles of a projected quad, each given twice, as a mesh's
    // top and bottom project once turned the same way: the union is the
    // quad, its corners untouched.
    let (p, q, r, s) = (
        Point2::new(0.1, 0.1),
        Point2::new(4.3, 0.1),
        Point2::new(4.3, 0.7),
        Point2::new(0.1, 0.7),
    );
    let tri = |a, b, c| Ring {
        points: vec![a, b, c],
    };
    let soup = [tri(p, q, r), tri(p, r, s), tri(q, r, p), tri(s, p, r)];
    let result = union_soup(&soup, tol()).unwrap();
    assert_eq!(result.len(), 1);
    same_points(&result[0].outer, &rect(0.1, 0.1, 4.3, 0.7));
}

#[test]
fn rectangles_sharing_a_corner_and_two_half_sides() {
    // [0,1]x[0,1] inside [0,2]x[0,1]: the bottom and top edges share an
    // end and run along each other, so each must be split where the
    // other ends.
    let (small, big) = (rect(0.0, 0.0, 1.0, 1.0), rect(0.0, 0.0, 2.0, 1.0));
    let cases = [
        (OverlayOperation::Intersection, rect(0.0, 0.0, 1.0, 1.0)),
        (OverlayOperation::Union, rect(0.0, 0.0, 2.0, 1.0)),
        (OverlayOperation::Xor, rect(1.0, 0.0, 2.0, 1.0)),
    ];
    for (operation, expected) in cases {
        let result = boolean(small.clone(), big.clone(), operation);
        assert_eq!(result.len(), 1, "{operation:?}");
        same_points(&result[0].outer, &expected);
    }
    let result = boolean(big, small, OverlayOperation::Difference);
    same_points(&result[0].outer, &rect(1.0, 0.0, 2.0, 1.0));
}
