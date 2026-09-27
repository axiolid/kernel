//! Many-target distance maps and the farthest point of a room (#186).
//!
//! Expected values are closed forms; the farthest-point interval must
//! contain them and be no wider than the tolerance asked for.

use axiolid_core::Point2;
use axiolid_overlay::{Polygon, Ring};
use axiolid_route::{
    distance_map, farthest_point, farthest_point_within, shortest_path, FarthestError, MapError,
    Unreachable,
};

fn p(x: f64, y: f64) -> Point2 {
    Point2::new(x, y)
}

fn ring(points: &[(f64, f64)]) -> Ring {
    Ring {
        points: points.iter().map(|(x, y)| p(*x, *y)).collect(),
    }
}

fn polygon(points: &[(f64, f64)]) -> Polygon {
    Polygon {
        outer: ring(points),
        holes: Vec::new(),
    }
}

fn rect(x0: f64, y0: f64, x1: f64, y1: f64) -> Polygon {
    polygon(&[(x0, y0), (x1, y0), (x1, y1), (x0, y1)])
}

const TOL: f64 = 1e-3;

fn assert_brackets(lower: f64, upper: f64, exact: f64) {
    assert!(
        lower <= exact && exact <= upper,
        "[{lower}, {upper}] misses {exact}"
    );
    assert!(
        upper - lower <= TOL * 1.01,
        "[{lower}, {upper}] wider than {TOL}"
    );
}

#[test]
fn one_door_the_far_corners() {
    let room = rect(0.0, 0.0, 10.0, 4.0);
    let map = distance_map(&[room.clone()], &[], &[p(0.0, 2.0)]).unwrap();
    let far = farthest_point(&map, &room, TOL).unwrap();
    assert!(far.converged);
    assert_brackets(far.distance.lower, far.distance.upper, 104f64.sqrt());
    let w = far.witness.unwrap();
    assert!(
        (w - p(10.0, 0.0)).length() < 1e-2 || (w - p(10.0, 4.0)).length() < 1e-2,
        "{w:?}"
    );
}

#[test]
fn two_doors_the_middle_of_the_long_walls() {
    // Nearest of two doors: the maximum is a crease, at (5, 0) and (5, 4).
    let room = rect(0.0, 0.0, 10.0, 4.0);
    let map = distance_map(&[room.clone()], &[], &[p(0.0, 2.0), p(10.0, 2.0)]).unwrap();
    let far = farthest_point(&map, &room, TOL).unwrap();
    assert!(far.converged);
    assert_brackets(far.distance.lower, far.distance.upper, 29f64.sqrt());
}

#[test]
fn an_l_shaped_room_is_measured_around_the_corner() {
    // The door is at the far end of the long arm; the farthest point is the
    // outer corner of the other arm, reached round the reflex corner.
    let room = polygon(&[
        (0.0, 0.0),
        (10.0, 0.0),
        (10.0, 10.0),
        (8.0, 10.0),
        (8.0, 2.0),
        (0.0, 2.0),
    ]);
    let map = distance_map(&[room.clone()], &[], &[p(0.0, 1.0)]).unwrap();
    let far = farthest_point(&map, &room, TOL).unwrap();
    assert!(far.converged);
    // Straight to (10, 10) would leave the room; round (8, 2) it is
    // sqrt(65) + sqrt(68). A straight-line bound would be sqrt(181).
    let exact = 65f64.sqrt() + 68f64.sqrt();
    assert_brackets(far.distance.lower, far.distance.upper, exact);
    assert!((far.witness.unwrap() - p(10.0, 10.0)).length() < 1e-2);
}

#[test]
fn a_subregion_is_measured_within_the_whole_floor() {
    // Floor with a pillar; the room is its right half; the exit is in the
    // left half, so the room's points are reached round the pillar.
    let floor = Polygon {
        outer: ring(&[(0.0, 0.0), (20.0, 0.0), (20.0, 10.0), (0.0, 10.0)]),
        holes: vec![ring(&[(9.0, 2.0), (9.0, 8.0), (11.0, 8.0), (11.0, 2.0)])],
    };
    let map = distance_map(&[floor], &[], &[p(0.0, 5.0)]).unwrap();
    let room = rect(10.0, 0.0, 20.0, 10.0);
    let far = farthest_point(&map, &room, TOL).unwrap();
    assert!(far.converged);
    // No point of the far wall is in sight of the door. Round the bottom
    // of the pillar or its top, the two ways tie at (20, 5):
    // sqrt(90) + 2 + sqrt(90).
    assert_brackets(
        far.distance.lower,
        far.distance.upper,
        2.0 * 90f64.sqrt() + 2.0,
    );
    assert!((far.witness.unwrap() - p(20.0, 5.0)).length() < 1e-2);
}

#[test]
fn a_barrier_lengthens_the_way() {
    // A free-standing wall at x = 5 from y = 1 to y = 9; the door at (0, 5).
    let room = rect(0.0, 0.0, 10.0, 10.0);
    let barrier = vec![vec![p(5.0, 1.0), p(5.0, 9.0)]];
    let map = distance_map(&[room.clone()], &barrier, &[p(0.0, 5.0)]).unwrap();
    let far = farthest_point(&map, &room, TOL).unwrap();
    assert!(far.converged);
    // (10, 5) is reached round either end of the wall: 2 sqrt(41).
    assert_brackets(far.distance.lower, far.distance.upper, 2.0 * 41f64.sqrt());
    assert!((far.witness.unwrap() - p(10.0, 5.0)).length() < 1e-2);
}

#[test]
fn the_far_side_of_a_barrier_is_not_bounded_from_the_near_side() {
    // A long wall close to the right of the room; the door on the left.
    // A point on the wall is in sight of the door, but the strip beyond it
    // is reached only round either end.
    let room = rect(0.0, 0.0, 6.0, 10.0);
    let barrier = vec![vec![p(5.0, 0.5), p(5.0, 8.5)]];
    let map = distance_map(&[room.clone()], &barrier, &[p(0.0, 5.0)]).unwrap();
    let far = farthest_point(&map, &room, TOL).unwrap();
    assert!(far.converged);
    // On the far wall x = 6 the ways round the two ends tie; bisect for
    // where (the maximum lies between sampled points, so the bracket's
    // upper end must come from the far side's own anchors).
    let (low_end, high_end) = (45.25f64.sqrt(), 37.25f64.sqrt());
    let below = |y: f64| low_end + (1.0 + (y - 0.5).powi(2)).sqrt();
    let above = |y: f64| high_end + (1.0 + (8.5 - y).powi(2)).sqrt();
    let (mut lo, mut hi) = (0.5, 8.5);
    for _ in 0..100 {
        let m = 0.5 * (lo + hi);
        if below(m) < above(m) {
            lo = m;
        } else {
            hi = m;
        }
    }
    let exact = below(lo);
    assert_brackets(far.distance.lower, far.distance.upper, exact);
}

#[test]
fn only_the_subregion_counts() {
    // The far end of the room is farther, but outside the subregion.
    let room = rect(0.0, 0.0, 10.0, 4.0);
    let map = distance_map(&[room], &[], &[p(0.0, 2.0)]).unwrap();
    let near = polygon(&[(0.0, 0.0), (5.0, 0.0), (5.0, 4.0), (0.0, 4.0)]);
    let far = farthest_point(&map, &near, TOL).unwrap();
    assert!(far.converged);
    assert_brackets(far.distance.lower, far.distance.upper, 29f64.sqrt());
}

#[test]
fn an_interior_maximum_is_found() {
    // Exits near the four corners: the farthest point is inside the room,
    // at the Voronoi vertex of the first three exits, where their
    // distances tie. The runner-up vertex is 0.0016 closer, so a bracket
    // to 0.001 must find this one. Value to 40 digits by mpmath.
    let room = rect(0.0, 0.0, 10.0, 7.0);
    let targets = [p(0.5, 1.0), p(9.0, 0.2), p(9.7, 6.1), p(1.3, 6.8)];
    let map = distance_map(&[room.clone()], &[], &targets).unwrap();
    let far = farthest_point(&map, &room, TOL).unwrap();
    assert!(far.converged);
    assert_brackets(
        far.distance.lower,
        far.distance.upper,
        5.075_490_229_455_75,
    );
    let w = far.witness.unwrap();
    assert!(
        (w - p(5.007_268_722_466_96, 3.333_480_176_211_45)).length() < 0.05,
        "{w:?}"
    );
}

#[test]
fn the_interval_contains_every_sampled_distance() {
    let room = polygon(&[
        (0.0, 0.0),
        (6.0, 0.0),
        (6.0, 3.0),
        (3.0, 3.0),
        (3.0, 6.0),
        (0.0, 6.0),
    ]);
    let map = distance_map(&[room.clone()], &[], &[p(6.0, 1.5), p(1.5, 6.0)]).unwrap();
    let far = farthest_point(&map, &room, TOL).unwrap();
    let mut sampled: f64 = 0.0;
    for i in 0..=120 {
        for j in 0..=120 {
            let q = p(6.0 * i as f64 / 120.0, 6.0 * j as f64 / 120.0);
            if let Ok(reach) = map.nearest(q).unwrap() {
                sampled = sampled.max(reach.route.length);
            }
        }
    }
    assert!(
        sampled <= far.distance.upper,
        "{sampled} > {:?}",
        far.distance
    );
    assert!(
        far.distance.lower >= sampled - 0.05,
        "{:?} vs {sampled}",
        far.distance
    );
}

#[test]
fn an_unreachable_pocket_is_refused_with_evidence() {
    let region = [rect(0.0, 0.0, 4.0, 4.0), rect(6.0, 0.0, 10.0, 4.0)];
    let map = distance_map(&region, &[], &[p(1.0, 1.0)]).unwrap();
    let whole = rect(0.0, 0.0, 10.0, 4.0);
    match farthest_point(&map, &whole, TOL) {
        Err(FarthestError::Unreachable { triangle }) => {
            assert!(triangle.iter().all(|c| c.x >= 6.0), "{triangle:?}");
        }
        other => panic!("{other:?}"),
    }
    // The reachable half alone is fine.
    let left = farthest_point(&map, &rect(0.0, 0.0, 4.0, 4.0), TOL).unwrap();
    assert_brackets(left.distance.lower, left.distance.upper, 18f64.sqrt());
    assert_eq!(
        map.nearest(p(7.0, 1.0)).unwrap(),
        Err(Unreachable::DisconnectedComponents)
    );
}

#[test]
fn nearest_agrees_with_shortest_path_to_each_target() {
    let floor = Polygon {
        outer: ring(&[(0.0, 0.0), (20.0, 0.0), (20.0, 10.0), (0.0, 10.0)]),
        holes: vec![ring(&[(9.0, 2.0), (9.0, 8.0), (11.0, 8.0), (11.0, 2.0)])],
    };
    let targets = [p(0.0, 5.0), p(20.0, 1.0), p(10.0, 9.0)];
    let region = [floor];
    let map = distance_map(&region, &[], &targets).unwrap();
    assert_eq!(map.targets(), 3);
    for q in [
        p(1.0, 1.0),
        p(15.0, 5.0),
        p(10.0, 1.0),
        p(12.0, 9.5),
        p(9.0, 5.0),
    ] {
        let reach = map.nearest(q).unwrap().unwrap();
        let best = targets
            .iter()
            .map(|t| shortest_path(&region, &[], q, *t).unwrap().unwrap().length)
            .fold(f64::INFINITY, f64::min);
        assert!(
            (reach.route.length - best).abs() < 1e-12,
            "{q:?}: {} vs {best}",
            reach.route.length
        );
        let own = shortest_path(&region, &[], q, targets[reach.target])
            .unwrap()
            .unwrap();
        assert!((own.length - best).abs() < 1e-12);
        assert_eq!(reach.route.polyline[0], q);
        assert_eq!(*reach.route.polyline.last().unwrap(), targets[reach.target]);
    }
    assert_eq!(
        map.nearest(p(10.0, 5.0)).unwrap(),
        Err(Unreachable::StartOutside)
    );
}

#[test]
fn malformed_maps_and_queries_are_refused() {
    let room = rect(0.0, 0.0, 4.0, 4.0);
    assert_eq!(
        distance_map(&[room.clone()], &[], &[]).unwrap_err(),
        MapError::NoTargets
    );
    assert_eq!(
        distance_map(&[room.clone()], &[], &[p(1.0, 1.0), p(5.0, 1.0)]).unwrap_err(),
        MapError::TargetOutside { index: 1 }
    );
    let map = distance_map(&[room.clone()], &[], &[p(1.0, 1.0)]).unwrap();
    assert_eq!(
        farthest_point(&map, &room, -1.0).unwrap_err(),
        FarthestError::InvalidTolerance
    );
    assert_eq!(
        farthest_point(&map, &rect(5.0, 5.0, 6.0, 6.0), TOL).unwrap_err(),
        FarthestError::Empty
    );
    let crossing = vec![
        vec![p(1.0, 0.5), p(3.0, 0.5)],
        vec![p(2.0, 0.0), p(2.0, 1.0)],
    ];
    let map = distance_map(&[room.clone()], &crossing, &[p(1.0, 1.0)]).unwrap();
    assert_eq!(
        farthest_point(&map, &room, TOL).unwrap_err(),
        FarthestError::CrossingObstacles
    );
}

#[test]
fn a_small_budget_still_brackets() {
    let room = rect(0.0, 0.0, 10.0, 4.0);
    let map = distance_map(&[room.clone()], &[], &[p(0.0, 2.0)]).unwrap();
    let far = farthest_point_within(&map, &room, 0.0, 12).unwrap();
    assert!(!far.converged);
    let exact = 104f64.sqrt();
    assert!(
        far.distance.lower <= exact && exact <= far.distance.upper,
        "{:?}",
        far.distance
    );
}
