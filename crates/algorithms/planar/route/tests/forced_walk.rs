//! The shortest walk forced through a region (#196).
//!
//! A ring corridor: the room [0, 10] x [0, 6] around the block
//! [2, 8] x [1, 4]. From (1, 1) to (9, 1) the short way runs along the
//! block's lower edge, length 8. The long way climbs to the block's upper
//! corners: sqrt(10) + 6 + sqrt(10).

use axiolid_core::Point2;
use axiolid_overlay::{Polygon, Ring};
use axiolid_route::{distance_map, forced_walk, forced_walk_within, FarthestError};

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

fn floor() -> Polygon {
    Polygon {
        outer: ring(&[(0.0, 0.0), (10.0, 0.0), (10.0, 6.0), (0.0, 6.0)]),
        holes: vec![ring(&[(2.0, 1.0), (2.0, 4.0), (8.0, 4.0), (8.0, 1.0)])],
    }
}

const TOL: f64 = 1e-3;

fn contains(lower: f64, upper: f64, exact: f64) {
    assert!(
        lower <= exact && exact <= upper,
        "[{lower}, {upper}] misses {exact}"
    );
    assert!(upper - lower <= TOL * 1.01, "[{lower}, {upper}] too wide");
}

#[test]
fn the_long_corridor_exceeds_the_shortest_walk_and_the_short_one_equals_it() {
    let region = [floor()];
    let from = distance_map(&region, &[], &[p(1.0, 1.0)]).unwrap();
    let to = distance_map(&region, &[], &[p(9.0, 1.0)]).unwrap();

    let short = forced_walk(&from, &to, &rect(2.0, 0.0, 8.0, 1.0), TOL).unwrap();
    assert!(short.converged);
    assert!((short.shortest - 8.0).abs() < 1e-12);
    contains(short.length.lower, short.length.upper, 8.0);

    let long = forced_walk(&from, &to, &rect(2.0, 4.0, 8.0, 6.0), TOL).unwrap();
    assert!(long.converged, "{long:?}");
    let exact = 6.0 + 2.0 * 10f64.sqrt();
    contains(long.length.lower, long.length.upper, exact);
    // The decision the bound exists for: no shortest walk goes this way.
    assert!(long.length.lower > long.shortest);
    // The witness is on the long way round, where the walk is shortest:
    // anywhere along the block's upper edge.
    let w = long.witness.unwrap();
    assert!(
        (w.y - 4.0).abs() < 1e-2 && (2.0..=8.0).contains(&w.x),
        "{w:?}"
    );
}

#[test]
fn a_polygon_holding_the_whole_shortest_walk_is_bounded_by_it() {
    let region = [floor()];
    let from = distance_map(&region, &[], &[p(1.0, 1.0)]).unwrap();
    let to = distance_map(&region, &[], &[p(9.0, 1.0)]).unwrap();
    let all = forced_walk(&from, &to, &rect(0.0, 0.0, 10.0, 6.0), TOL).unwrap();
    assert!(all.length.lower <= all.shortest);
    contains(all.length.lower, all.length.upper, 8.0);
}

#[test]
fn a_walk_forced_into_a_far_pocket_turns_back() {
    // An open room; origin and target at the south wall. Entering the
    // square [4, 6] x [4, 6] costs twice the way to its nearest point.
    let region = [rect(0.0, 0.0, 10.0, 10.0)];
    let from = distance_map(&region, &[], &[p(3.0, 0.0)]).unwrap();
    let to = distance_map(&region, &[], &[p(7.0, 0.0)]).unwrap();
    let walk = forced_walk(&from, &to, &rect(4.0, 4.0, 6.0, 6.0), TOL).unwrap();
    assert!(walk.converged);
    // Best point: the square's lower edge; by symmetry its middle (5, 4),
    // at sqrt(4 + 16) from each end.
    contains(walk.length.lower, walk.length.upper, 2.0 * 20f64.sqrt());
    assert!((walk.shortest - 4.0).abs() < 1e-12);
}

#[test]
fn several_origins_and_targets_take_the_best_pair() {
    let region = [rect(0.0, 0.0, 10.0, 10.0)];
    let from = distance_map(&region, &[], &[p(0.0, 5.0), p(10.0, 5.0)]).unwrap();
    let to = distance_map(&region, &[], &[p(5.0, 0.0)]).unwrap();
    // Through the square [8, 9] x [4, 6]: from (10, 5) to its corner
    // (9, 5)... the nearest point to both is on its edge x = 9 or corner
    // (8, 4); the sum is minimised on the segment from (10, 5) to (5, 0),
    // which passes (9, 4): the walk is the straight line, 5 sqrt 2.
    let walk = forced_walk(&from, &to, &rect(8.0, 4.0, 9.0, 6.0), TOL).unwrap();
    contains(walk.length.lower, walk.length.upper, 50f64.sqrt());
}

#[test]
fn a_polygon_no_walk_reaches_is_infinitely_far() {
    // Two rooms with no door; the polygon lies in the second.
    let region = [rect(0.0, 0.0, 4.0, 4.0), rect(6.0, 0.0, 10.0, 4.0)];
    let from = distance_map(&region, &[], &[p(1.0, 1.0)]).unwrap();
    let to = distance_map(&region, &[], &[p(3.0, 3.0)]).unwrap();
    let walk = forced_walk(&from, &to, &rect(7.0, 1.0, 8.0, 2.0), TOL).unwrap();
    assert!(walk.length.lower.is_infinite() && walk.witness.is_none());
    // Straddling both rooms: only the reachable part counts.
    let walk = forced_walk(&from, &to, &rect(3.0, 1.0, 7.0, 2.0), TOL).unwrap();
    assert!(walk.length.upper.is_finite());
}

#[test]
fn mismatched_and_malformed_queries_are_refused() {
    let region = [floor()];
    let from = distance_map(&region, &[], &[p(1.0, 1.0)]).unwrap();
    let other = distance_map(&[rect(0.0, 0.0, 10.0, 6.0)], &[], &[p(9.0, 1.0)]).unwrap();
    assert_eq!(
        forced_walk(&from, &other, &rect(0.0, 0.0, 1.0, 1.0), TOL).unwrap_err(),
        FarthestError::MismatchedMaps
    );
    let barrier = vec![vec![p(0.5, 5.0), p(1.5, 5.0)]];
    let walled = distance_map(&region, &barrier, &[p(9.0, 1.0)]).unwrap();
    assert_eq!(
        forced_walk(&from, &walled, &rect(0.0, 0.0, 1.0, 1.0), TOL).unwrap_err(),
        FarthestError::MismatchedMaps
    );
    let to = distance_map(&region, &[], &[p(9.0, 1.0)]).unwrap();
    assert_eq!(
        forced_walk(&from, &to, &rect(0.0, 0.0, 1.0, 1.0), -1.0).unwrap_err(),
        FarthestError::InvalidTolerance
    );
    // Inside the block only: no free space.
    assert_eq!(
        forced_walk(&from, &to, &rect(3.0, 2.0, 4.0, 3.0), TOL).unwrap_err(),
        FarthestError::Empty
    );
}

#[test]
fn a_small_budget_still_brackets() {
    let region = [floor()];
    let from = distance_map(&region, &[], &[p(1.0, 1.0)]).unwrap();
    let to = distance_map(&region, &[], &[p(9.0, 1.0)]).unwrap();
    let exact = 6.0 + 2.0 * 10f64.sqrt();
    let walk = forced_walk_within(&from, &to, &rect(2.0, 4.0, 8.0, 6.0), 0.0, 40).unwrap();
    assert!(!walk.converged);
    assert!(walk.length.lower <= exact && exact <= walk.length.upper);
}

/// The bracket against dense sampling: no sampled walk through the polygon
/// is shorter than the lower end, and the best one is near the upper end.
fn check_against_samples(
    region: &[Polygon],
    barriers: &[Vec<Point2>],
    origin: Point2,
    target: Point2,
    through: &Polygon,
) {
    let from = distance_map(region, barriers, &[origin]).unwrap();
    let to = distance_map(region, barriers, &[target]).unwrap();
    let walk = forced_walk(&from, &to, through, TOL).unwrap();
    assert!(walk.converged, "{walk:?}");
    let [x0, y0, x1, y1] = {
        let pts = &through.outer.points;
        let xs = pts.iter().map(|p| p.x);
        let ys = pts.iter().map(|p| p.y);
        [
            xs.clone().fold(f64::INFINITY, f64::min),
            ys.clone().fold(f64::INFINITY, f64::min),
            xs.fold(f64::NEG_INFINITY, f64::max),
            ys.fold(f64::NEG_INFINITY, f64::max),
        ]
    };
    let n = 60;
    let mut least = f64::INFINITY;
    for i in 0..=n {
        for j in 0..=n {
            let q = p(
                x0 + (x1 - x0) * f64::from(i) / f64::from(n),
                y0 + (y1 - y0) * f64::from(j) / f64::from(n),
            );
            // A point on a barrier takes each map's nearer side, so the
            // sum there is no walk: it would pass through the barrier.
            let on_barrier = barriers.iter().any(|b| {
                b.windows(2).any(|w| {
                    let (u, v) = (w[0], w[1]);
                    let cross = (v - u).perp_dot(q - u);
                    cross == 0.0
                        && q.x >= u.x.min(v.x)
                        && q.x <= u.x.max(v.x)
                        && q.y >= u.y.min(v.y)
                        && q.y <= u.y.max(v.y)
                })
            });
            if on_barrier {
                continue;
            }
            let (Ok(Ok(a)), Ok(Ok(b))) = (from.nearest(q), to.nearest(q)) else {
                continue;
            };
            let f = a.route.length + b.route.length;
            assert!(
                walk.length.lower <= f,
                "lower {} above a walk of {f} through {q:?}",
                walk.length.lower
            );
            least = least.min(f);
        }
    }
    // The grid's spacing bounds how far its best walk can miss the least.
    let spacing = ((x1 - x0).max(y1 - y0)) / f64::from(n);
    assert!(
        least <= walk.length.upper + 2.0 * spacing * std::f64::consts::SQRT_2,
        "sampled {least}, bracket {:?}",
        walk.length
    );
}

#[test]
fn the_bracket_holds_against_sampled_walks() {
    // An L-shaped room around a corner.
    let l = Polygon {
        outer: ring(&[
            (0.0, 0.0),
            (6.0, 0.0),
            (6.0, 2.0),
            (2.0, 2.0),
            (2.0, 6.0),
            (0.0, 6.0),
        ]),
        holes: Vec::new(),
    };
    check_against_samples(
        &[l.clone()],
        &[],
        p(5.5, 1.0),
        p(1.0, 5.5),
        &rect(0.0, 0.0, 2.0, 2.0),
    );
    check_against_samples(
        &[l],
        &[],
        p(5.5, 1.0),
        p(1.0, 5.5),
        &rect(3.0, 0.0, 4.0, 1.0),
    );
    // A barrier across most of an open room.
    let barrier = vec![vec![p(0.0, 5.0), p(7.0, 5.0)]];
    let room = [rect(0.0, 0.0, 10.0, 10.0)];
    check_against_samples(
        &room,
        &barrier,
        p(1.0, 1.0),
        p(1.0, 9.0),
        &rect(2.0, 6.0, 4.0, 8.0),
    );
    check_against_samples(
        &room,
        &barrier,
        p(1.0, 1.0),
        p(1.0, 9.0),
        &rect(0.0, 0.0, 10.0, 10.0),
    );
    // Two holes, the polygon between them.
    let holes = [Polygon {
        outer: ring(&[(0.0, 0.0), (12.0, 0.0), (12.0, 8.0), (0.0, 8.0)]),
        holes: vec![
            ring(&[(2.0, 2.0), (5.0, 2.0), (5.0, 6.0), (2.0, 6.0)]),
            ring(&[(7.0, 2.0), (10.0, 2.0), (10.0, 6.0), (7.0, 6.0)]),
        ],
    }];
    check_against_samples(
        &holes,
        &[],
        p(1.0, 4.0),
        p(11.0, 4.0),
        &rect(5.0, 3.0, 7.0, 5.0),
    );
    check_against_samples(
        &holes,
        &[],
        p(1.0, 1.0),
        p(11.0, 7.0),
        &rect(5.0, 0.0, 7.0, 8.0),
    );
}
