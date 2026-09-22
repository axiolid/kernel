use axiolid_core::Point2;
use axiolid_reference::{minimum_area_rectangle, side_lengths, strict_convex_hull};

fn p(x: f64, y: f64) -> Point2 {
    Point2::new(x, y)
}

#[test]
fn hull_is_stable_and_strict() {
    let points = [
        p(0., 0.),
        p(2., 0.),
        p(2., 2.),
        p(0., 2.),
        p(1., 0.),
        p(0., 0.),
    ];
    assert_eq!(strict_convex_hull(&points).unwrap(), vec![0, 1, 2, 3]);
}

#[test]
fn rectangle_has_expected_dimensions() {
    let points = [p(0., 0.), p(3., 0.), p(3., 2.), p(0., 2.)];
    let rectangle = minimum_area_rectangle(&points).unwrap();
    assert_eq!(rectangle.area(), 6.0);
    assert_eq!(side_lengths(&rectangle), [2.0, 3.0]);
}

/// The rectangle really encloses the input, corner for corner.
///
/// `minimum_area_rectangle` now returns the shared `Rectangle2` instead of a
/// local four-corner struct. Checking the area alone would not notice a
/// rectangle built from the right dimensions in the wrong place, which is the
/// mistake the origin-plus-edge-vectors form makes easy to introduce.
#[test]
fn rectangle_encloses_its_input_points() {
    // Deliberately offset from the origin AND rotated. An axis-aligned box at
    // the origin cannot distinguish a correctly placed rectangle from one
    // built at the origin by mistake, because there both are the same answer.
    let points = [p(11., 10.), p(14., 14.), p(10., 17.), p(7., 13.)];
    let rectangle = minimum_area_rectangle(&points).unwrap();

    let corners = rectangle.corners();
    for point in points {
        let hit = corners
            .iter()
            .any(|corner| (*corner - point).length() < 1e-9);
        assert!(hit, "{point:?} is not a corner of {corners:?}");
    }
}

/// A rotated square still yields its own area, not its axis-aligned bounds.
///
/// The diamond's axis-aligned box has twice the area, so this fails loudly if
/// the caliper ever collapses to an `Aabb2`.
#[test]
fn rotated_input_keeps_the_oriented_area() {
    let diamond = [p(1., 0.), p(2., 1.), p(1., 2.), p(0., 1.)];
    let rectangle = minimum_area_rectangle(&diamond).unwrap();
    assert!(
        (rectangle.area() - 2.0).abs() < 1e-9,
        "expected the rotated square's own area, got {}",
        rectangle.area()
    );
}

/// The caliper keeps the SMALLEST enclosing rectangle, not merely a valid one.
///
/// Every hull edge yields an enclosing rectangle, so a search that kept the
/// largest still returns something plausible-looking with the right corners.
/// Only comparing against the alternatives catches an inverted comparison.
#[test]
fn rectangle_is_minimal_across_every_hull_edge() {
    // An obtuse triangle, chosen because its three hull edges give genuinely
    // different rectangles: 12, 14.4 and 12.414. A right triangle is useless
    // here -- all three of its edges give exactly 12, so keeping the largest
    // would pass.
    let triangle = [p(0., 0.), p(6., 0.), p(5., 2.)];
    let rectangle = minimum_area_rectangle(&triangle).unwrap();
    assert!(
        (rectangle.area() - 12.0).abs() < 1e-9,
        "expected the minimal 12, got {}",
        rectangle.area()
    );
}

#[test]
fn non_finite_input_is_rejected() {
    assert!(strict_convex_hull(&[p(0., 0.), p(f64::NAN, 0.)]).is_err());
}

fn turn(a: (i64, i64), b: (i64, i64), c: (i64, i64)) -> i128 {
    i128::from(b.0 - a.0) * i128::from(c.1 - a.1) - i128::from(b.1 - a.1) * i128::from(c.0 - a.0)
}

/// Independent Jarvis march over exact integer coordinates. It deliberately
/// shares neither Axiolid's predicate implementation nor its monotone-chain
/// control flow.
fn integer_oracle(points: &[(i64, i64)]) -> Vec<usize> {
    let start = (0..points.len()).min_by_key(|&i| (points[i], i)).unwrap();
    let mut result = vec![start];
    loop {
        let current = *result.last().unwrap();
        let mut next = (0..points.len())
            .find(|&i| points[i] != points[current])
            .unwrap();
        for candidate in 0..points.len() {
            if candidate == current || points[candidate] == points[current] {
                continue;
            }
            let side = turn(points[current], points[next], points[candidate]);
            let farther = (points[candidate].0 - points[current].0).pow(2)
                + (points[candidate].1 - points[current].1).pow(2)
                > (points[next].0 - points[current].0).pow(2)
                    + (points[next].1 - points[current].1).pow(2);
            if side < 0 || (side == 0 && farther) {
                next = candidate;
            }
        }
        if next == start {
            break;
        }
        result.push(next);
    }
    result
}

#[test]
fn certified_hull_matches_independent_integer_oracle() {
    let integers = [
        (2, 0),
        (0, 2),
        (3, 3),
        (0, 0),
        (3, 0),
        (1, 1),
        (0, 0),
        (0, 3),
    ];
    let points: Vec<_> = integers
        .iter()
        .map(|&(x, y)| p(x as f64, y as f64))
        .collect();
    assert_eq!(
        strict_convex_hull(&points).unwrap(),
        integer_oracle(&integers)
    );
}
