//! Weighted maps whose targets start at their own costs (#198).
//!
//! The room `[0, 10] x [0, 4]` with a stair `[4, 6] x [0, 4]` at factor 2
//! across it, an exit at `(0, 2)` and a landing at `(10, 2)` carrying 5 of
//! walk beyond it. On the middle line every walk crosses the stair square
//! to its edges, so the distances are exact: to the exit `x + 2` beyond the
//! stair, to the landing `15 - x`; they meet at `x = 6.5`.

use axiolid_core::Point2;
use axiolid_overlay::{Polygon, Ring};
use axiolid_route::{weighted_distance_map, weighted_distance_map_seeded, CostRegion, MapError};

const SPACING: f64 = 0.1;

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

fn room() -> [Polygon; 1] {
    [rect(0.0, 0.0, 10.0, 4.0)]
}

fn stair() -> [CostRegion; 1] {
    [CostRegion::new(rect(4.0, 0.0, 6.0, 4.0), 2.0)]
}

/// The cost to the exit and to the landing, with its weight, from a point
/// of the middle line.
fn exact(x: f64) -> (f64, f64) {
    let inside = |a: f64, b: f64| (b.min(6.0) - a.max(4.0)).max(0.0);
    let to_exit = x + inside(0.0, x);
    let to_landing = (10.0 - x) + inside(x, 10.0) + 5.0;
    (to_exit, to_landing)
}

#[test]
fn the_nearest_target_switches_where_the_costs_differ_by_the_weight() {
    let map = weighted_distance_map_seeded(
        &room(),
        &[],
        &[(p(0.0, 2.0), 0.0), (p(10.0, 2.0), 5.0)],
        &stair(),
        SPACING,
    )
    .unwrap();
    for x in [0.5, 3.0, 5.0, 6.4, 6.6, 8.0, 9.5] {
        let reach = map.nearest(p(x, 2.0)).unwrap().unwrap();
        let (exit, landing) = exact(x);
        let (target, cost) = if exit <= landing {
            (0, exit)
        } else {
            (1, landing)
        };
        assert_eq!(reach.target, target, "x = {x}");
        assert!(
            reach.cost.lower <= cost + 1e-9 && cost <= reach.cost.upper + 1e-9,
            "x = {x}: {:?} misses {cost}",
            reach.cost
        );
        assert!(reach.cost.upper - reach.cost.lower <= 1e-9, "{reach:?}");
        // The route runs straight along the line; its length is its own,
        // without the weight or the factor.
        let length = if target == 0 { x } else { 10.0 - x };
        assert!((reach.route.length - length).abs() < 1e-9, "{reach:?}");
    }
}

#[test]
fn a_heavy_target_reached_from_a_lighter_one_names_the_lighter() {
    // A landing 1 m from the exit, 5 behind: never the nearer.
    let map = weighted_distance_map_seeded(
        &room(),
        &[],
        &[(p(1.0, 2.0), 5.0), (p(0.0, 2.0), 0.0)],
        &stair(),
        SPACING,
    )
    .unwrap();
    for q in [p(1.0, 2.0), p(2.0, 3.0), p(9.0, 1.0)] {
        assert_eq!(map.nearest(q).unwrap().unwrap().target, 1, "{q:?}");
    }
}

#[test]
fn with_every_weight_zero_it_is_the_unseeded_map() {
    let targets = [p(0.0, 2.0), p(10.0, 1.0), p(5.0, 3.5)];
    let plain = weighted_distance_map(&room(), &[], &targets, &stair(), 0.25).unwrap();
    let zero: Vec<(Point2, f64)> = targets.iter().map(|t| (*t, 0.0)).collect();
    let seeded = weighted_distance_map_seeded(&room(), &[], &zero, &stair(), 0.25).unwrap();
    for i in 0..=20 {
        for j in 0..=8 {
            let q = p(f64::from(i) * 0.5, f64::from(j) * 0.5);
            assert_eq!(
                plain.nearest(q).unwrap(),
                seeded.nearest(q).unwrap(),
                "{q:?}"
            );
        }
    }
}

#[test]
fn a_weight_negative_or_not_finite_is_refused() {
    for w in [-1.0, f64::NAN, f64::INFINITY] {
        assert!(matches!(
            weighted_distance_map_seeded(
                &room(),
                &[],
                &[(p(0.0, 2.0), 0.0), (p(10.0, 2.0), w)],
                &stair(),
                SPACING,
            ),
            Err(MapError::InvalidWeight { index: 1 })
        ));
    }
}
