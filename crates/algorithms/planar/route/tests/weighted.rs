//! Weighted distance maps (#195).
//!
//! Expected values are closed forms; every bracket must contain them. The
//! gap between the bounds is exact where every crossing of a cost edge is
//! square to it, and otherwise shrinks with the spacing of the points
//! along cost edges.

use axiolid_core::Point2;
use axiolid_overlay::{Polygon, Ring};
use axiolid_route::{
    distance_map, weighted_distance_map, weighted_farthest_point, CostRegion, MapError,
};

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

fn brackets(lower: f64, upper: f64, exact: f64, width: f64) {
    assert!(
        lower <= exact + 1e-12 && exact <= upper + 1e-12,
        "[{lower}, {upper}] misses {exact}"
    );
    assert!(
        upper - lower <= width,
        "[{lower}, {upper}] wider than {width}"
    );
}

#[test]
fn a_corridor_crossed_by_a_stair_costs_its_length_twice() {
    // The corridor [0, 10] x [0, 1]; the stair [3, 7] x [0, 1] at factor 2
    // spans it wall to wall. From (10, 0.5) to (0, 0.5): 10 + 4.
    let corridor = [rect(0.0, 0.0, 10.0, 1.0)];
    let stair = CostRegion::new(rect(3.0, 0.0, 7.0, 1.0), 2.0);
    let map = weighted_distance_map(&corridor, &[], &[p(0.0, 0.5)], &[stair], 0.1).unwrap();
    let reach = map.nearest(p(10.0, 0.5)).unwrap().unwrap();
    brackets(reach.cost.lower, reach.cost.upper, 14.0, 1e-9);
    assert_eq!(reach.target, 0);
    // The walk: straight along the corridor.
    assert!(
        (reach.route.length - 10.0).abs() < 1e-9,
        "{:?}",
        reach.route
    );
    // A point inside the stair.
    let reach = map.nearest(p(5.0, 0.5)).unwrap().unwrap();
    brackets(reach.cost.lower, reach.cost.upper, 3.0 + 4.0, 1e-9);
}

#[test]
fn a_detour_around_the_stair_is_taken_when_cheaper() {
    // The room [0, 10] x [0, 4]; the stair [4, 6] x [0, 3] leaves a gap
    // above it. From (0, 1) to (10, 1): straight through costs 10 + 2 at
    // factor 2, round the top over its corners 2 sqrt 20 + 2.
    let room = [rect(0.0, 0.0, 10.0, 4.0)];
    let detour = 2.0 * 20f64.sqrt() + 2.0;
    let stair = CostRegion::new(rect(4.0, 0.0, 6.0, 3.0), 2.0);
    let spacing = 0.05;
    let map = weighted_distance_map(&room, &[], &[p(0.0, 1.0)], &[stair], spacing).unwrap();
    let reach = map.nearest(p(10.0, 1.0)).unwrap().unwrap();
    // The walk is exact, round both corners; the lower bound may cut them
    // by about an interval.
    assert!((reach.cost.upper - detour).abs() < 1e-9);
    brackets(reach.cost.lower, reach.cost.upper, detour, 2.0 * spacing);
    assert!(
        reach.route.polyline.contains(&p(4.0, 3.0)) && reach.route.polyline.contains(&p(6.0, 3.0)),
        "{:?}",
        reach.route.polyline
    );
    // At factor 1.2 straight through is cheaper: 10 + 0.2 x 2.
    let light = CostRegion::new(rect(4.0, 0.0, 6.0, 3.0), 1.2);
    let map = weighted_distance_map(&room, &[], &[p(0.0, 1.0)], &[light], spacing).unwrap();
    let reach = map.nearest(p(10.0, 1.0)).unwrap().unwrap();
    brackets(reach.cost.lower, reach.cost.upper, 10.4, 1e-9);
}

#[test]
fn a_walk_refracts_across_a_cost_edge() {
    // The half-plane y < 2 of the room [0, 10] x [0, 4] costs 2. From
    // (0, 0) to (6, 4) the walk bends where it crosses y = 2, by Snell's
    // law: sin(a1) / sin(a2) = 1 / 2 ... solved numerically below.
    let room = [rect(0.0, 0.0, 10.0, 4.0)];
    let cheap = CostRegion::new(rect(0.0, 0.0, 10.0, 2.0), 2.0);
    let exact = {
        // Minimise 2 sqrt(x^2 + 4) + sqrt((6 - x)^2 + 4) over x by
        // bisection on the derivative.
        let g =
            |x: f64| 2.0 * x / (x * x + 4.0).sqrt() - (6.0 - x) / ((6.0 - x).powi(2) + 4.0).sqrt();
        let (mut lo, mut hi) = (0.0, 6.0);
        for _ in 0..200 {
            let mid = 0.5 * (lo + hi);
            if g(mid) > 0.0 {
                hi = mid;
            } else {
                lo = mid;
            }
        }
        let x = 0.5 * (lo + hi);
        2.0 * (x * x + 4.0).sqrt() + ((6.0 - x).powi(2) + 4.0).sqrt()
    };
    for spacing in [0.2, 0.05] {
        let map =
            weighted_distance_map(&room, &[], &[p(0.0, 0.0)], &[cheap.clone()], spacing).unwrap();
        let reach = map.nearest(p(6.0, 4.0)).unwrap().unwrap();
        // The gap shrinks with the spacing.
        brackets(reach.cost.lower, reach.cost.upper, exact, 2.0 * spacing);
    }
}

#[test]
fn without_costs_it_is_the_plain_distance() {
    let room = Polygon {
        outer: ring(&[(0.0, 0.0), (10.0, 0.0), (10.0, 6.0), (0.0, 6.0)]),
        holes: vec![ring(&[(2.0, 1.0), (2.0, 4.0), (8.0, 4.0), (8.0, 1.0)])],
    };
    let targets = [p(1.0, 1.0), p(9.0, 5.0)];
    let plain = distance_map(&[room.clone()], &[], &targets).unwrap();
    let weighted = weighted_distance_map(&[room.clone()], &[], &targets, &[], 1.0).unwrap();
    for q in [p(5.0, 0.5), p(5.0, 5.5), p(0.5, 5.5), p(9.5, 0.5)] {
        let d = plain.nearest(q).unwrap().unwrap().route.length;
        let reach = weighted.nearest(q).unwrap().unwrap();
        brackets(reach.cost.lower, reach.cost.upper, d, 1e-9);
    }
    // A factor-1 cost region changes nothing either.
    // Its points along the edges let the lower bound cut a corner by up to
    // an interval where a walk passes them.
    let one = CostRegion::new(rect(3.0, 4.0, 5.0, 6.0), 1.0);
    let weighted = weighted_distance_map(&[room], &[], &targets, &[one], 0.1).unwrap();
    for q in [p(5.0, 0.5), p(5.0, 5.5), p(0.5, 5.5)] {
        let d = plain.nearest(q).unwrap().unwrap().route.length;
        let reach = weighted.nearest(q).unwrap().unwrap();
        brackets(reach.cost.lower, reach.cost.upper, d, 0.1);
        assert!((reach.cost.upper - d).abs() < 1e-9);
    }
}

#[test]
fn overlapping_regions_take_the_greatest_factor() {
    // Nested: [3, 7] at 2 holds [4, 6] at 3 across the corridor. From
    // (10, 0.5) to (0, 0.5): 6 + 2 x 2 + 2 x 3 = 16.
    let corridor = [rect(0.0, 0.0, 10.0, 1.0)];
    let outer = CostRegion::new(rect(3.0, 0.0, 7.0, 1.0), 2.0);
    let inner = CostRegion::new(rect(4.0, 0.0, 6.0, 1.0), 3.0);
    let map = weighted_distance_map(&corridor, &[], &[p(0.0, 0.5)], &[outer, inner], 0.1).unwrap();
    let reach = map.nearest(p(10.0, 0.5)).unwrap().unwrap();
    brackets(reach.cost.lower, reach.cost.upper, 16.0, 1e-9);
    // The lesser inside the greater counts the greater: [4, 6] at 3 holds
    // [4.5, 5.5] at 2.
    let big = CostRegion::new(rect(4.0, 0.0, 6.0, 1.0), 3.0);
    let small = CostRegion::new(rect(4.5, 0.0, 5.5, 1.0), 2.0);
    let map = weighted_distance_map(&corridor, &[], &[p(0.0, 0.5)], &[big, small], 0.1).unwrap();
    let reach = map.nearest(p(10.0, 0.5)).unwrap().unwrap();
    brackets(reach.cost.lower, reach.cost.upper, 8.0 + 6.0, 1e-9);
}

#[test]
fn the_farthest_point_is_bracketed_by_weighted_distance() {
    // The corridor [0, 10] x [0, 1] with the stair [3, 7] across, factor
    // 2, and the exit at (0, 0.5): the far corners (10, 0) and (10, 1) are
    // 10 + 4 + ... sqrt(100.25) straight, weighted: the walk crosses the
    // stair on a slant. Bracket against sampled bounds.
    let corridor = [rect(0.0, 0.0, 10.0, 1.0)];
    let stair = CostRegion::new(rect(3.0, 0.0, 7.0, 1.0), 2.0);
    let map = weighted_distance_map(&corridor, &[], &[p(0.0, 0.5)], &[stair], 0.1).unwrap();
    let far = weighted_farthest_point(&map, &corridor[0], 0.05).unwrap();
    // No sampled point's lower bound exceeds the upper end; the best
    // sampled upper bound is no less than the lower end.
    let mut best_upper = 0.0f64;
    for i in 0..=100 {
        for j in 0..=10 {
            let q = p(f64::from(i) * 0.1, f64::from(j) * 0.1);
            let reach = map.nearest(q).unwrap().unwrap();
            assert!(
                reach.cost.lower <= far.distance.upper,
                "{q:?}: {:?} vs {far:?}",
                reach.cost
            );
            best_upper = best_upper.max(reach.cost.upper);
        }
    }
    assert!(far.distance.lower <= best_upper);
    assert!(far.witness.unwrap().x > 9.9);
}

#[test]
fn malformed_costs_are_refused() {
    let room = [rect(0.0, 0.0, 10.0, 4.0)];
    let t = [p(1.0, 1.0)];
    let bad = CostRegion::new(rect(4.0, 0.0, 6.0, 3.0), 0.5);
    assert_eq!(
        weighted_distance_map(&room, &[], &t, &[bad], 0.1).unwrap_err(),
        MapError::InvalidFactor { index: 0 }
    );
    let ok = CostRegion::new(rect(4.0, 0.0, 6.0, 3.0), 2.0);
    for spacing in [0.0, -1.0, f64::NAN, f64::INFINITY] {
        assert_eq!(
            weighted_distance_map(&room, &[], &t, &[ok.clone()], spacing).unwrap_err(),
            MapError::InvalidSpacing
        );
    }
    // Poking out through the wall.
    let out = CostRegion::new(rect(4.0, 3.0, 6.0, 5.0), 2.0);
    assert_eq!(
        weighted_distance_map(&room, &[], &t, &[ok.clone(), out], 0.1).unwrap_err(),
        MapError::CostCrossing { index: 1 }
    );
    // Crossing another cost region.
    let cross = CostRegion::new(rect(5.0, 1.0, 7.0, 2.0), 2.0);
    assert!(matches!(
        weighted_distance_map(&room, &[], &t, &[ok.clone(), cross], 0.1).unwrap_err(),
        MapError::CostCrossing { .. }
    ));
    // Along a barrier.
    let barrier = vec![vec![p(4.0, 3.0), p(6.0, 3.0)]];
    assert_eq!(
        weighted_distance_map(&room, &barrier, &t, &[ok], 0.1).unwrap_err(),
        MapError::CostCrossing { index: 0 }
    );
}

/// Two maps of the same scene at different spacings bracket the same
/// distance, so their brackets must overlap at every point; and no weighted
/// distance is below the plain one. A lower bound that cut a corner it may
/// not would show as a coarse lower end above the fine upper end.
fn cross_check(
    region: &[Polygon],
    barriers: &[Vec<Point2>],
    targets: &[Point2],
    costs: &[CostRegion],
    queries: &[Point2],
) {
    let coarse = weighted_distance_map(region, barriers, targets, costs, 0.4).unwrap();
    let fine = weighted_distance_map(region, barriers, targets, costs, 0.05).unwrap();
    let plain = distance_map(region, barriers, targets).unwrap();
    let mut checked = 0;
    for &q in queries {
        let (Ok(Ok(c)), Ok(Ok(f)), Ok(Ok(d))) =
            (coarse.nearest(q), fine.nearest(q), plain.nearest(q))
        else {
            continue;
        };
        assert!(c.cost.lower <= c.cost.upper && f.cost.lower <= f.cost.upper);
        assert!(
            c.cost.lower <= f.cost.upper + 1e-9 && f.cost.lower <= c.cost.upper + 1e-9,
            "{q:?}: coarse {:?}, fine {:?}",
            c.cost,
            f.cost
        );
        assert!(f.cost.upper >= d.route.length - 1e-9, "{q:?}");
        // The fine map is the tighter.
        assert!(f.cost.upper - f.cost.lower <= c.cost.upper - c.cost.lower + 1e-9);
        checked += 1;
    }
    assert!(
        checked > queries.len() / 2,
        "{checked} of {}",
        queries.len()
    );
}

fn grid(x0: f64, y0: f64, x1: f64, y1: f64, n: u32) -> Vec<Point2> {
    let mut out = Vec::new();
    for i in 0..=n {
        for j in 0..=n {
            out.push(p(
                x0 + (x1 - x0) * f64::from(i) / f64::from(n),
                y0 + (y1 - y0) * f64::from(j) / f64::from(n),
            ));
        }
    }
    out
}

#[test]
fn brackets_at_two_spacings_overlap() {
    // An L-shaped floor, a stair in its corner, a barrier.
    let l = Polygon {
        outer: ring(&[
            (0.0, 0.0),
            (8.0, 0.0),
            (8.0, 3.0),
            (3.0, 3.0),
            (3.0, 8.0),
            (0.0, 8.0),
        ]),
        holes: Vec::new(),
    };
    let stair = CostRegion::new(rect(0.5, 0.5, 2.5, 2.5), 2.0);
    let barrier = vec![vec![p(4.0, 0.0), p(4.0, 2.0)]];
    cross_check(
        &[l],
        &barrier,
        &[p(7.5, 0.5), p(1.5, 7.5)],
        &[stair],
        &grid(0.1, 0.1, 7.9, 7.9, 12),
    );
    // A room around a pillar; two cost regions sharing an edge, and one
    // with a hole.
    let room = Polygon {
        outer: ring(&[(0.0, 0.0), (10.0, 0.0), (10.0, 6.0), (0.0, 6.0)]),
        holes: vec![ring(&[(4.5, 2.5), (5.5, 2.5), (5.5, 3.5), (4.5, 3.5)])],
    };
    let a = CostRegion::new(rect(1.0, 1.0, 3.0, 5.0), 1.5);
    let b = CostRegion::new(rect(3.0, 1.0, 4.0, 5.0), 2.5);
    let ring_region = CostRegion::new(
        Polygon {
            outer: ring(&[(6.0, 1.0), (9.0, 1.0), (9.0, 5.0), (6.0, 5.0)]),
            holes: vec![ring(&[(7.0, 2.0), (8.0, 2.0), (8.0, 4.0), (7.0, 4.0)])],
        },
        3.0,
    );
    cross_check(
        &[room],
        &[],
        &[p(0.5, 3.0)],
        &[a, b, ring_region],
        &grid(0.2, 0.2, 9.8, 5.8, 12),
    );
    // A cost region along a wall.
    let hall = [rect(0.0, 0.0, 12.0, 3.0)];
    let bench = CostRegion::new(rect(2.0, 0.0, 10.0, 1.0), 2.0);
    cross_check(
        &hall,
        &[],
        &[p(0.0, 0.5), p(12.0, 2.5)],
        &[bench],
        &grid(0.1, 0.1, 11.9, 2.9, 12),
    );
}

#[test]
fn brackets_overlap_where_cost_regions_are_awkward() {
    let room = [rect(0.0, 0.0, 12.0, 8.0)];
    // A U-shaped region, open to the north: a hop across its notch is not
    // inside it though its corners are.
    let u = CostRegion::new(
        Polygon {
            outer: ring(&[
                (2.0, 1.0),
                (8.0, 1.0),
                (8.0, 6.0),
                (6.5, 6.0),
                (6.5, 2.5),
                (3.5, 2.5),
                (3.5, 6.0),
                (2.0, 6.0),
            ]),
            holes: Vec::new(),
        },
        3.0,
    );
    cross_check(
        &room,
        &[],
        // One target in the U's west arm: walks from the east arm cross
        // the notch.
        &[p(5.0, 0.5), p(11.5, 7.5), p(2.75, 4.0)],
        &[u],
        &grid(0.1, 0.1, 11.9, 7.9, 14),
    );
    // Regions touching corner to corner and corner to edge, a pillar's
    // corner on a cost edge, and a barrier ending on one.
    let pillar = Polygon {
        outer: ring(&[(0.0, 0.0), (12.0, 0.0), (12.0, 8.0), (0.0, 8.0)]),
        holes: vec![ring(&[(4.0, 4.0), (5.0, 4.0), (5.0, 5.0), (4.0, 5.0)])],
    };
    let a = CostRegion::new(rect(1.0, 1.0, 4.0, 4.0), 2.0);
    let b = CostRegion::new(rect(4.0, 0.5, 7.0, 3.0), 1.5);
    let c = CostRegion::new(rect(7.0, 2.0, 10.0, 6.0), 2.5);
    let barrier = vec![vec![p(8.5, 8.0), p(8.5, 6.0)]];
    cross_check(
        &[pillar],
        &barrier,
        &[p(0.5, 7.5), p(11.5, 0.5)],
        &[a, b, c],
        &grid(0.1, 0.1, 11.9, 7.9, 14),
    );
}

#[test]
fn a_coarse_farthest_point_still_holds_the_steepest_slope() {
    // The farthest point of a square inside a factor-3 stair, whose far
    // corners are no vertex of the free space: a coarse tolerance stops
    // with large cells, whose bounds must use the stair's factor.
    let corridor = [rect(0.0, 0.0, 10.0, 2.0)];
    let stair = CostRegion::new(rect(3.0, 0.0, 9.0, 2.0), 3.0);
    let map = weighted_distance_map(&corridor, &[], &[p(0.0, 1.0)], &[stair], 0.25).unwrap();
    let square = rect(5.0, 0.5, 7.5, 1.5);
    let far = weighted_farthest_point(&map, &square, 2.0).unwrap();
    for q in [p(7.5, 0.5), p(7.5, 1.5), p(7.5, 1.0)] {
        let reach = map.nearest(q).unwrap().unwrap();
        assert!(
            reach.cost.lower <= far.distance.upper,
            "{q:?}: {:?} above {:?}",
            reach.cost,
            far.distance
        );
    }
}

#[test]
fn turned_scenes_keep_their_brackets() {
    // The scenes above, turned: cost edges at an angle, whose points along
    // them are interpolated off their lines by rounding, and cost regions
    // meeting the walls only up to rounding (#198). Every bracket must
    // still hold the closed form.
    let refraction = {
        let g =
            |x: f64| 2.0 * x / (x * x + 4.0).sqrt() - (6.0 - x) / ((6.0 - x).powi(2) + 4.0).sqrt();
        let (mut lo, mut hi) = (0.0, 6.0);
        for _ in 0..200 {
            let mid = 0.5 * (lo + hi);
            if g(mid) > 0.0 {
                hi = mid;
            } else {
                lo = mid;
            }
        }
        let x = 0.5 * (lo + hi);
        2.0 * (x * x + 4.0).sqrt() + ((6.0 - x).powi(2) + 4.0).sqrt()
    };
    // (room, cost region, factor, target, query, exact, gap allowed)
    let scenes = [
        (
            (0.0, 0.0, 10.0, 1.0),
            (3.0, 0.0, 7.0, 1.0),
            2.0,
            (0.0, 0.5),
            (10.0, 0.5),
            14.0,
            0.1,
        ),
        (
            (0.0, 0.0, 10.0, 4.0),
            (4.0, 0.0, 6.0, 3.0),
            2.0,
            (0.0, 1.0),
            (10.0, 1.0),
            2.0 * 20f64.sqrt() + 2.0,
            0.2,
        ),
        (
            (0.0, 0.0, 10.0, 4.0),
            (0.0, 0.0, 10.0, 2.0),
            2.0,
            (0.0, 0.0),
            (6.0, 4.0),
            refraction,
            0.2,
        ),
        // Clear of the walls: round the square's corners.
        (
            (0.0, 0.0, 10.0, 4.0),
            (4.0, 1.0, 6.0, 3.0),
            2.0,
            (0.0, 2.0),
            (10.0, 2.0),
            2.0 * 17f64.sqrt() + 2.0,
            0.2,
        ),
    ];
    // Turned by the 3-4-5 angle and scaled by 5, every corner stays on
    // the integer grid: the walls and cost edges meet exactly, and only
    // the points along the cost edges are rounded. At other angles the
    // corners round too; the square clear of the walls, in a larger room
    // so that no query lies on a wall, is turned by those.
    let mut cases: Vec<(f64, f64, usize)> = (0..scenes.len()).map(|k| (4.0, 3.0, k)).collect();
    for degrees in [30.0f64, 17.0, 45.0] {
        let (s, c) = degrees.to_radians().sin_cos();
        cases.push((c, s, 3));
    }
    for (c, s, k) in cases {
        let scale = (c * c + s * s).sqrt();
        let turn = |x: f64, y: f64| p(c * x - s * y, s * x + c * y);
        let turned = |(x0, y0, x1, y1): (f64, f64, f64, f64)| Polygon {
            outer: Ring {
                points: vec![turn(x0, y0), turn(x1, y0), turn(x1, y1), turn(x0, y1)],
            },
            holes: Vec::new(),
        };
        let (mut room, cost, factor, from, to, exact, gap) = scenes[k];
        if k == 3 {
            room = (room.0 - 1.0, room.1 - 1.0, room.2 + 1.0, room.3 + 1.0);
        }
        let map = weighted_distance_map(
            &[turned(room)],
            &[],
            &[turn(from.0, from.1)],
            &[CostRegion::new(turned(cost), factor)],
            0.05 * scale,
        )
        .unwrap_or_else(|e| panic!("scene {k} turned ({c}, {s}): {e:?}"));
        let reach = map.nearest(turn(to.0, to.1)).unwrap().unwrap();
        let (lower, upper) = (reach.cost.lower / scale, reach.cost.upper / scale);
        assert!(
            lower <= exact + 1e-9 && exact <= upper + 1e-9,
            "scene {k} turned ({c}, {s}): [{lower}, {upper}] misses {exact}"
        );
        assert!(
            upper - lower <= gap,
            "scene {k} turned ({c}, {s}): [{lower}, {upper}]"
        );
    }
}
