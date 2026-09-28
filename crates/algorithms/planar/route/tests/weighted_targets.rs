//! Distance maps whose targets start at their own distances (#197).
//!
//! A floor with two ways out: an exit at (0, 2), weight 0, and a stair
//! landing at (10, 2) carrying 5 m of walk beyond it. The distance is the
//! least of `|p - exit|` and `|p - landing| + 5`.

use axiolid_core::Point2;
use axiolid_overlay::{Polygon, Ring};
use axiolid_route::{distance_map, distance_map_weighted, farthest_point, forced_walk, MapError};

fn p(x: f64, y: f64) -> Point2 {
    Point2::new(x, y)
}

fn rect(x0: f64, y0: f64, x1: f64, y1: f64) -> Polygon {
    Polygon {
        outer: Ring {
            points: vec![p(x0, y0), p(x1, y0), p(x1, y1), p(x0, y1)],
        },
        holes: Vec::new(),
    }
}

fn floor() -> [Polygon; 1] {
    [rect(0.0, 0.0, 10.0, 4.0)]
}

const EXIT: (f64, f64) = (0.0, 2.0);
const LANDING: (f64, f64) = (10.0, 2.0);

fn exact(q: Point2) -> (usize, f64) {
    let a = (q - p(EXIT.0, EXIT.1)).length();
    let b = (q - p(LANDING.0, LANDING.1)).length() + 5.0;
    if a <= b {
        (0, a)
    } else {
        (1, b)
    }
}

#[test]
fn the_nearest_target_switches_where_the_lengths_differ_by_the_weight() {
    let map = distance_map_weighted(
        &floor(),
        &[],
        &[(p(EXIT.0, EXIT.1), 0.0), (p(LANDING.0, LANDING.1), 5.0)],
    )
    .unwrap();
    // On the middle line the switch is where x = (10 - x) + 5: x = 7.5.
    for (x, target) in [(7.4, 0), (7.49, 0), (7.51, 1), (7.6, 1)] {
        let reach = map.nearest(p(x, 2.0)).unwrap().unwrap();
        assert_eq!(reach.target, target, "x = {x}");
        let (_, d) = exact(p(x, 2.0));
        assert!(
            (reach.distance - d).abs() < 1e-12,
            "{} vs {d}",
            reach.distance
        );
        // The route's own length leaves the weight out.
        let beyond = if target == 1 { 5.0 } else { 0.0 };
        assert!((reach.route.length + beyond - reach.distance).abs() < 1e-12);
    }
    // Everywhere else too.
    for i in 0..=20 {
        for j in 0..=8 {
            let q = p(f64::from(i) * 0.5, f64::from(j) * 0.5);
            let reach = map.nearest(q).unwrap().unwrap();
            let (target, d) = exact(q);
            assert!((reach.distance - d).abs() < 1e-12, "{q:?}");
            if (d - (q - p(EXIT.0, EXIT.1)).length()).abs() > 1e-9
                || (d - (q - p(LANDING.0, LANDING.1)).length() - 5.0).abs() > 1e-9
            {
                assert_eq!(reach.target, target, "{q:?}");
            }
        }
    }
}

#[test]
fn a_heavy_target_reached_from_a_lighter_one_names_the_lighter() {
    // A landing weighted 20 m, 10 m from the exit: from the landing itself
    // the exit is nearer.
    let map =
        distance_map_weighted(&floor(), &[], &[(p(0.0, 2.0), 0.0), (p(10.0, 2.0), 20.0)]).unwrap();
    let reach = map.nearest(p(10.0, 2.0)).unwrap().unwrap();
    assert_eq!(reach.target, 0);
    assert!((reach.distance - 10.0).abs() < 1e-12);
    assert_eq!(reach.route.polyline.last(), Some(&p(0.0, 2.0)));
}

#[test]
fn all_weights_zero_is_the_plain_map() {
    let targets = [p(0.0, 2.0), p(10.0, 2.0), p(4.0, 4.0)];
    let plain = distance_map(&floor(), &[], &targets).unwrap();
    let weighted: Vec<(Point2, f64)> = targets.iter().map(|t| (*t, 0.0)).collect();
    let zero = distance_map_weighted(&floor(), &[], &weighted).unwrap();
    for i in 0..=20 {
        for j in 0..=8 {
            let q = p(f64::from(i) * 0.5, f64::from(j) * 0.5);
            let a = plain.nearest(q).unwrap().unwrap();
            let b = zero.nearest(q).unwrap().unwrap();
            assert_eq!(a, b, "{q:?}");
            assert_eq!(a.distance, a.route.length);
        }
    }
    let room = rect(0.0, 0.0, 10.0, 4.0);
    assert_eq!(
        farthest_point(&plain, &room, 1e-3).unwrap(),
        farthest_point(&zero, &room, 1e-3).unwrap()
    );
}

#[test]
fn the_farthest_point_is_bracketed_with_weights() {
    let map = distance_map_weighted(
        &floor(),
        &[],
        &[(p(EXIT.0, EXIT.1), 0.0), (p(LANDING.0, LANDING.1), 5.0)],
    )
    .unwrap();
    let tol = 1e-3;
    let far = farthest_point(&map, &floor()[0], tol).unwrap();
    assert!(far.converged);
    // The greatest distance lies on the switching curve, farthest from
    // the axis: on the walls y = 0 and y = 4, where
    // sqrt(x^2 + 4) = sqrt((10 - x)^2 + 4) + 5.
    let g = |x: f64| (x * x + 4.0).sqrt() - (((10.0 - x) * (10.0 - x) + 4.0).sqrt() + 5.0);
    let (mut lo, mut hi) = (0.0, 10.0);
    for _ in 0..200 {
        let mid = 0.5 * (lo + hi);
        if g(mid) < 0.0 {
            lo = mid;
        } else {
            hi = mid;
        }
    }
    let x = 0.5 * (lo + hi);
    let max = (x * x + 4.0).sqrt();
    assert!(
        far.distance.lower <= max && max <= far.distance.upper,
        "{:?} misses {max}",
        far.distance
    );
    assert!(far.distance.upper - far.distance.lower <= tol * 1.01);
    let w = far.witness.unwrap();
    assert!(exact(w).1 >= far.distance.lower, "{w:?}");
}

#[test]
fn a_forced_walk_counts_the_origins_weight() {
    // Walks from a landing carrying 5 m to the exit.
    let from = distance_map_weighted(&floor(), &[], &[(p(10.0, 2.0), 5.0)]).unwrap();
    let to = distance_map(&floor(), &[], &[p(0.0, 2.0)]).unwrap();
    let walk = forced_walk(&from, &to, &rect(4.0, 1.0, 6.0, 3.0), 1e-3).unwrap();
    assert!((walk.shortest - 15.0).abs() < 1e-12);
    assert!(walk.length.lower <= 15.0 && 15.0 <= walk.length.upper);
}

#[test]
fn negative_and_non_finite_weights_are_refused() {
    for (w, index) in [(-1.0, 1), (f64::NAN, 1), (f64::INFINITY, 1)] {
        assert_eq!(
            distance_map_weighted(&floor(), &[], &[(p(0.0, 2.0), 0.0), (p(10.0, 2.0), w)])
                .unwrap_err(),
            MapError::InvalidWeight { index }
        );
    }
}

#[test]
fn of_two_targets_on_one_point_the_lighter_counts() {
    let map =
        distance_map_weighted(&floor(), &[], &[(p(0.0, 2.0), 5.0), (p(0.0, 2.0), 1.0)]).unwrap();
    let reach = map.nearest(p(3.0, 2.0)).unwrap().unwrap();
    assert_eq!(reach.target, 1);
    assert!((reach.distance - 4.0).abs() < 1e-12);
}

#[test]
fn a_route_bending_at_a_heavy_target_names_the_target_it_ends_at() {
    // A pillar [4, 6] x [1, 3]; a landing weighted 100 m at its corner
    // (6, 3). From (7, 2) the way to the exit bends round that corner,
    // through the landing's own graph vertex, and ends at the exit.
    let floor = [Polygon {
        outer: rect(0.0, 0.0, 10.0, 4.0).outer,
        holes: vec![rect(4.0, 1.0, 6.0, 3.0).outer],
    }];
    let map =
        distance_map_weighted(&floor, &[], &[(p(0.0, 2.0), 0.0), (p(6.0, 3.0), 100.0)]).unwrap();
    let reach = map.nearest(p(7.0, 2.0)).unwrap().unwrap();
    assert_eq!(reach.target, 0, "{reach:?}");
    assert!(reach.route.polyline.contains(&p(6.0, 3.0)));
    let exact = 2f64.sqrt() + 2.0 + 17f64.sqrt();
    assert!((reach.distance - exact).abs() < 1e-12, "{}", reach.distance);
}
