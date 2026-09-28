//! The cheapest walk forced through a region, over weighted maps (#198).
//!
//! A corridor `[0, 10] x [0, 1]` with an alcove `[8, 9] x [1, 3]` above
//! it, and a stair `[3, 7] x [0, 1]` at factor 2 across the corridor. From
//! the corner `(0, 1)` to the corner `(10, 1)` the cheapest walk runs along
//! the top wall, `3 + 2 * 4 + 3 = 14`. The cheapest one that reaches the
//! alcove's back half `[8, 9] x [2, 3]` leaves the corridor at `(8, 1)`,
//! touches `(8.5, 2)` and comes back by `(9, 1)`: `11 + 2 sqrt 1.25 + 1 =
//! 13 + sqrt 5`.

use axiolid_core::Point2;
use axiolid_overlay::{Polygon, Ring};
use axiolid_route::{
    distance_map, forced_walk, weighted_distance_map, weighted_distance_map_seeded,
    weighted_forced_walk, weighted_forced_walk_within, CostRegion, FarthestError, WeightedMap,
};

const TOL: f64 = 1e-3;
const SPACING: f64 = 0.1;

fn p(x: f64, y: f64) -> Point2 {
    Point2::new(x, y)
}

fn ring(points: &[(f64, f64)]) -> Ring {
    Ring {
        points: points.iter().map(|(x, y)| p(*x, *y)).collect(),
    }
}

fn rect(x0: f64, y0: f64, x1: f64, y1: f64) -> Polygon {
    Polygon {
        outer: ring(&[(x0, y0), (x1, y0), (x1, y1), (x0, y1)]),
        holes: Vec::new(),
    }
}

fn corridor() -> [Polygon; 1] {
    [Polygon {
        outer: ring(&[
            (0.0, 0.0),
            (10.0, 0.0),
            (10.0, 1.0),
            (9.0, 1.0),
            (9.0, 3.0),
            (8.0, 3.0),
            (8.0, 1.0),
            (0.0, 1.0),
        ]),
        holes: Vec::new(),
    }]
}

fn stair() -> Vec<CostRegion> {
    vec![CostRegion::new(rect(3.0, 0.0, 7.0, 1.0), 2.0)]
}

fn map(target: Point2) -> WeightedMap {
    weighted_distance_map(&corridor(), &[], &[target], &stair(), SPACING).unwrap()
}

#[test]
fn the_alcove_costs_more_than_the_cheapest_walk() {
    let (from, to) = (map(p(0.0, 1.0)), map(p(10.0, 1.0)));
    let back = weighted_forced_walk(&from, &to, &rect(8.0, 2.0, 9.0, 3.0), TOL).unwrap();
    let exact = 13.0 + 5f64.sqrt();
    assert!(back.converged, "{back:?}");
    assert!(
        back.cost.lower <= exact && exact <= back.cost.upper,
        "{back:?} misses {exact}"
    );
    assert!(back.cost.upper - back.cost.lower <= 2.0 * TOL, "{back:?}");
    assert!(back.shortest.lower <= 14.0 && 14.0 <= back.shortest.upper);
    // The decision: no cheapest walk reaches the alcove's back.
    assert!(back.cost.lower > back.shortest.upper);
    let w = back.witness.unwrap();
    assert!(
        (w.x - 8.5).abs() < 0.05 && (w.y - 2.0).abs() < 1e-2,
        "{w:?}"
    );

    // The stair lies on the cheapest walk: forcing through it costs
    // nothing more.
    let on = weighted_forced_walk(&from, &to, &rect(4.0, 0.0, 6.0, 1.0), TOL).unwrap();
    assert!(on.converged, "{on:?}");
    assert!(on.cost.lower <= 14.0 && 14.0 <= on.cost.upper, "{on:?}");
    assert!(on.cost.lower <= on.shortest.upper);
}

#[test]
fn without_costs_it_brackets_the_plain_forced_walk() {
    let region = corridor();
    let through = rect(8.0, 2.0, 9.0, 3.0);
    let plain = forced_walk(
        &distance_map(&region, &[], &[p(0.0, 1.0)]).unwrap(),
        &distance_map(&region, &[], &[p(10.0, 1.0)]).unwrap(),
        &through,
        TOL,
    )
    .unwrap();
    let from = weighted_distance_map(&region, &[], &[p(0.0, 1.0)], &[], SPACING).unwrap();
    let to = weighted_distance_map(&region, &[], &[p(10.0, 1.0)], &[], SPACING).unwrap();
    let weighted = weighted_forced_walk(&from, &to, &through, TOL).unwrap();
    assert!(weighted.converged);
    let exact = 9.0 + 5f64.sqrt();
    for (lower, upper) in [
        (plain.length.lower, plain.length.upper),
        (weighted.cost.lower, weighted.cost.upper),
    ] {
        assert!(lower <= exact && exact <= upper, "[{lower}, {upper}]");
        assert!(upper - lower <= 2.0 * TOL, "[{lower}, {upper}]");
    }
}

#[test]
fn the_bracket_holds_the_sampled_walks() {
    // A polygon across the stair's corner, where walks refract: every
    // sampled point's walk bounds the answer from above, and the sampled
    // lower bounds, less the most they can fall between samples, from
    // below.
    let (from, to) = (map(p(0.0, 0.0)), map(p(10.0, 1.0)));
    let through = rect(5.5, 0.25, 7.5, 0.75);
    let forced = weighted_forced_walk_within(&from, &to, &through, 1e-2, 4000).unwrap();
    let n = 40;
    let h = 2.0 / f64::from(n);
    let (mut least_upper, mut least_lower) = (f64::INFINITY, f64::INFINITY);
    for i in 0..=n {
        for j in 0..=n {
            let q = p(
                5.5 + 2.0 * f64::from(i) / f64::from(n),
                0.25 + 0.5 * f64::from(j) / f64::from(n),
            );
            let a = from.nearest(q).unwrap().unwrap().cost;
            let b = to.nearest(q).unwrap().unwrap().cost;
            least_upper = least_upper.min(a.upper + b.upper);
            least_lower = least_lower.min(a.lower + b.lower);
        }
    }
    assert!(
        forced.cost.lower <= least_upper,
        "{forced:?} vs {least_upper}"
    );
    // Two distances at factor at most 2, within h of a sample.
    assert!(
        forced.cost.upper >= least_lower - 2.0 * 2.0 * h,
        "{forced:?} vs {least_lower}"
    );
}

#[test]
fn origins_start_at_their_weights() {
    // The one origin carries 3 of walk behind it: every cost rises by 3.
    let seeded =
        weighted_distance_map_seeded(&corridor(), &[], &[(p(0.0, 1.0), 3.0)], &stair(), SPACING)
            .unwrap();
    let to = map(p(10.0, 1.0));
    let back = weighted_forced_walk(&seeded, &to, &rect(8.0, 2.0, 9.0, 3.0), TOL).unwrap();
    let exact = 3.0 + 13.0 + 5f64.sqrt();
    assert!(
        back.cost.lower <= exact && exact <= back.cost.upper,
        "{back:?}"
    );
    assert!(back.shortest.lower <= 17.0 && 17.0 <= back.shortest.upper);
}

#[test]
fn maps_over_different_costs_are_refused() {
    let from = map(p(0.0, 1.0));
    let to = weighted_distance_map(&corridor(), &[], &[p(10.0, 1.0)], &[], SPACING).unwrap();
    assert!(matches!(
        weighted_forced_walk(&from, &to, &rect(8.0, 2.0, 9.0, 3.0), TOL),
        Err(FarthestError::MismatchedMaps)
    ));
}
