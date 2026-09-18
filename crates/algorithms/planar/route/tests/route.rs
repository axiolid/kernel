//! Shortest-path fixtures: known lengths, refusals, determinism (#46).
//!
//! Lengths are computed from the geometry by hand. Where a detour is forced,
//! the expected value is the exact sum of Euclidean legs, so a path that cuts
//! a corner it should not fails loudly.

use axiolid_core::Point2;
use axiolid_overlay::{Polygon, Ring};
use axiolid_route::{shortest_path, shortest_path_within, RouteError, Unreachable, MAX_VERTICES};

fn ring(points: &[(f64, f64)]) -> Ring {
    Ring {
        points: points.iter().map(|(x, y)| Point2::new(*x, *y)).collect(),
    }
}

/// A 10x10 open room.
fn room() -> Vec<Polygon> {
    vec![Polygon {
        outer: ring(&[(0.0, 0.0), (10.0, 0.0), (10.0, 10.0), (0.0, 10.0)]),
        holes: Vec::new(),
    }]
}

#[test]
fn an_unobstructed_path_is_the_straight_line() {
    // No obstacle, so the shortest path is the direct segment: length 10.
    let route = shortest_path(&room(), &[], Point2::new(1.0, 5.0), Point2::new(9.0, 5.0))
        .expect("valid query")
        .expect("both endpoints are inside an empty room");

    assert!(
        (route.length - 8.0).abs() < 1e-12,
        "expected 8, got {}",
        route.length
    );
    assert_eq!(route.polyline.len(), 2, "a clear line needs no waypoints");
    assert_eq!(route.polyline[0], Point2::new(1.0, 5.0));
    assert_eq!(route.polyline[1], Point2::new(9.0, 5.0));
}

#[test]
fn a_hole_forces_a_detour_of_known_length() {
    // Room with a 4x4 pillar centred on the straight line from (1,5) to (9,5).
    // The direct route is blocked, so the path bends around a pillar corner.
    let region = vec![Polygon {
        outer: ring(&[(0.0, 0.0), (10.0, 0.0), (10.0, 10.0), (0.0, 10.0)]),
        holes: vec![ring(&[(3.0, 3.0), (3.0, 7.0), (7.0, 7.0), (7.0, 3.0)])],
    }];

    let route = shortest_path(&region, &[], Point2::new(1.0, 5.0), Point2::new(9.0, 5.0))
        .expect("valid query")
        .expect("the pillar can be walked around");

    // Via a corner: (1,5)->(3,7)->(7,7)->(9,5) is 2*sqrt(8) + 4.
    let expected = 2.0 * 8f64.sqrt() + 4.0;
    assert!(
        (route.length - expected).abs() < 1e-9,
        "expected {expected}, got {}",
        route.length
    );
    // Strictly longer than the blocked straight line, which is 8.
    assert!(route.length > 8.0);
    assert!(route.polyline.len() >= 3, "a detour needs waypoints");
}

#[test]
fn a_zero_width_barrier_blocks_without_bounding_area() {
    // A wall from (5,0) to (5,8) leaves a 2-unit gap at the top. It has no
    // area at all, so it can only affect the route through visibility.
    let barrier = vec![vec![Point2::new(5.0, 0.0), Point2::new(5.0, 8.0)]];

    let start = Point2::new(1.0, 4.0);
    let goal = Point2::new(9.0, 4.0);
    let direct = shortest_path(&room(), &[], start, goal)
        .expect("valid")
        .expect("no barrier, clear line");
    assert!((direct.length - 8.0).abs() < 1e-12);

    let around = shortest_path(&room(), &barrier, start, goal)
        .expect("valid")
        .expect("the gap above the wall is passable");

    // Must route over the wall tip at (5,8): 2*sqrt(16+16) = 2*sqrt(32).
    let expected = 2.0 * 32f64.sqrt();
    assert!(
        (around.length - expected).abs() < 1e-9,
        "expected {expected} around the wall tip, got {}",
        around.length
    );
    assert!(
        around.length > direct.length,
        "a zero-width wall must still lengthen the route"
    );
}

#[test]
fn endpoints_outside_the_region_are_named_individually() {
    let outside = Point2::new(50.0, 50.0);
    let inside = Point2::new(5.0, 5.0);

    // Which endpoint is at fault is actionable information, so the two are
    // distinct variants rather than one "outside" answer.
    assert_eq!(
        shortest_path(&room(), &[], outside, inside).expect("valid query"),
        Err(Unreachable::StartOutside)
    );
    assert_eq!(
        shortest_path(&room(), &[], inside, outside).expect("valid query"),
        Err(Unreachable::GoalOutside)
    );
}

#[test]
fn disconnected_rooms_are_reported_as_disconnected() {
    // Two rooms with no doorway. Both endpoints are inside the region, so
    // this is genuinely a connectivity fact, not a containment one.
    let region = vec![
        Polygon {
            outer: ring(&[(0.0, 0.0), (4.0, 0.0), (4.0, 4.0), (0.0, 4.0)]),
            holes: Vec::new(),
        },
        Polygon {
            outer: ring(&[(20.0, 0.0), (24.0, 0.0), (24.0, 4.0), (20.0, 4.0)]),
            holes: Vec::new(),
        },
    ];

    let result = shortest_path(&region, &[], Point2::new(2.0, 2.0), Point2::new(22.0, 2.0))
        .expect("valid query");
    assert_eq!(result, Err(Unreachable::DisconnectedComponents));
}

#[test]
fn equal_length_paths_resolve_deterministically() {
    // A symmetric pillar offers two mirror-image detours of identical length.
    // Repeating the query must return the same one every time, otherwise the
    // answer depends on iteration order.
    let region = vec![Polygon {
        outer: ring(&[(0.0, 0.0), (10.0, 0.0), (10.0, 10.0), (0.0, 10.0)]),
        holes: vec![ring(&[(4.0, 3.0), (4.0, 7.0), (6.0, 7.0), (6.0, 3.0)])],
    }];
    let start = Point2::new(1.0, 5.0);
    let goal = Point2::new(9.0, 5.0);

    let first = shortest_path(&region, &[], start, goal)
        .expect("valid")
        .expect("routable");
    for _ in 0..8 {
        let again = shortest_path(&region, &[], start, goal)
            .expect("valid")
            .expect("routable");
        assert_eq!(first.polyline, again.polyline, "tie-breaking is unstable");
    }
}

/// A circular ring with `count` vertices, centred on (50, 50).
fn big_ring(count: usize) -> Ring {
    let mut points = Vec::new();
    for index in 0..count {
        let angle = (index as f64) * core::f64::consts::TAU / (count as f64);
        points.push((50.0 + 40.0 * angle.cos(), 50.0 + 40.0 * angle.sin()));
    }
    ring(&points)
}

#[test]
fn oversized_input_is_refused_rather_than_truncated() {
    // A ring with more vertices than the documented bound. Truncating would
    // answer a different question without telling the caller.
    let mut points = Vec::new();
    let count = MAX_VERTICES + 8;
    for index in 0..count {
        let angle = (index as f64) * core::f64::consts::TAU / (count as f64);
        points.push((50.0 + 40.0 * angle.cos(), 50.0 + 40.0 * angle.sin()));
    }
    let region = vec![Polygon {
        outer: ring(&points),
        holes: Vec::new(),
    }];

    let error = shortest_path(
        &region,
        &[],
        Point2::new(50.0, 50.0),
        Point2::new(51.0, 50.0),
    )
    .expect_err("the bound must be enforced");
    assert!(matches!(error, RouteError::TooManyVertices { .. }));
}

#[test]
fn no_shorter_path_exists_by_brute_force_enumeration() {
    // Independent check: enumerate every path of up to three legs through the
    // pillar corners and confirm none beats the reported length. This does not
    // reuse the visibility graph -- it re-derives legality from the geometry.
    let hole = [(3.0, 3.0), (3.0, 7.0), (7.0, 7.0), (7.0, 3.0)];
    let region = vec![Polygon {
        outer: ring(&[(0.0, 0.0), (10.0, 0.0), (10.0, 10.0), (0.0, 10.0)]),
        holes: vec![ring(&hole)],
    }];
    let start = Point2::new(1.0, 5.0);
    let goal = Point2::new(9.0, 5.0);

    let route = shortest_path(&region, &[], start, goal)
        .expect("valid")
        .expect("a route exists");

    let corners: Vec<Point2> = hole.iter().map(|(x, y)| Point2::new(*x, *y)).collect();

    // A leg is legal when it stays clear of the pillar interior. Sampled
    // densely and independently of the visibility predicate under test.
    let legal = |a: Point2, b: Point2| -> bool {
        (1..200).all(|step| {
            let t = f64::from(step) / 200.0;
            let x = a.x + (b.x - a.x) * t;
            let y = a.y + (b.y - a.y) * t;
            let inside_pillar =
                x > 3.0 + 1e-9 && x < 7.0 - 1e-9 && y > 3.0 + 1e-9 && y < 7.0 - 1e-9;
            let inside_room = (0.0..=10.0).contains(&x) && (0.0..=10.0).contains(&y);
            !inside_pillar && inside_room
        })
    };

    let dist = |a: Point2, b: Point2| ((b.x - a.x).powi(2) + (b.y - a.y).powi(2)).sqrt();
    let mut best = f64::INFINITY;
    if legal(start, goal) {
        best = dist(start, goal);
    }
    for a in &corners {
        if legal(start, *a) && legal(*a, goal) {
            best = best.min(dist(start, *a) + dist(*a, goal));
        }
        for b in &corners {
            if legal(start, *a) && legal(*a, *b) && legal(*b, goal) {
                best = best.min(dist(start, *a) + dist(*a, *b) + dist(*b, goal));
            }
        }
    }

    assert!(
        route.length <= best + 1e-9,
        "brute force found a shorter path than the graph"
    );
    assert!(
        best.is_finite(),
        "brute force must find at least one legal path"
    );
}

/// A refusal must carry the bound and the budget, not just a complaint.
///
/// The endpoints are 1.0 apart, so no route between them can be
/// shorter than 1.0 however the region is shaped (kernel#92).
#[test]
fn an_oversized_refusal_carries_a_proven_lower_bound() {
    let region = vec![Polygon {
        outer: big_ring(MAX_VERTICES + 8),
        holes: Vec::new(),
    }];
    let error = shortest_path(
        &region,
        &[],
        Point2::new(50.0, 50.0),
        Point2::new(51.0, 50.0),
    )
    .expect_err("the budget must be enforced");

    let RouteError::TooManyVertices {
        supplied,
        budget,
        lower_bound,
    } = error
    else {
        panic!("expected a budget refusal, got {error:?}");
    };
    assert!(supplied > budget, "{supplied} must exceed {budget}");
    assert_eq!(budget, MAX_VERTICES);
    assert_eq!(lower_bound, 1.0);
}

/// A raised budget accepts input the default refuses.
///
/// The consumer case from kernel#92: a caller with a larger time
/// budget should not have to shrink its supported model size.
#[test]
fn a_raised_budget_accepts_what_the_default_refuses() {
    let region = vec![Polygon {
        outer: big_ring(MAX_VERTICES + 8),
        holes: Vec::new(),
    }];
    let start = Point2::new(50.0, 50.0);
    let goal = Point2::new(51.0, 50.0);

    assert!(
        shortest_path(&region, &[], start, goal).is_err(),
        "the default budget must still refuse"
    );

    let route = shortest_path_within(&region, &[], start, goal, MAX_VERTICES * 2)
        .expect("a raised budget admits this input")
        .expect("both endpoints are inside a convex region");
    // Unobstructed inside a convex ring, so the route is the straight
    // line and the earlier lower bound is exactly attained.
    assert_eq!(route.length, 1.0);
}

/// The bound stays valid when obstacles force a detour.
///
/// This is what makes it a BOUND rather than an estimate. A barrier
/// between the endpoints lengthens the real route; the straight-line
/// distance must still be less than or equal to it, never above.
#[test]
fn the_lower_bound_never_exceeds_the_real_route() {
    let region = vec![Polygon {
        outer: big_ring(MAX_VERTICES + 8),
        holes: Vec::new(),
    }];
    // A wall across the direct line, open at one end only.
    let barriers = vec![vec![Point2::new(50.5, 20.0), Point2::new(50.5, 70.0)]];
    let start = Point2::new(40.0, 50.0);
    let goal = Point2::new(60.0, 50.0);

    let error = shortest_path(&region, &barriers, start, goal)
        .expect_err("the default budget must refuse this input");
    let RouteError::TooManyVertices { lower_bound, .. } = error else {
        panic!("expected a budget refusal, got {error:?}");
    };

    let route = shortest_path_within(&region, &barriers, start, goal, MAX_VERTICES * 2)
        .expect("a raised budget admits this input")
        .expect("the wall is open at one end, so a route exists");

    assert!(
        lower_bound <= route.length,
        "bound {lower_bound} must not exceed the real route {}",
        route.length
    );
    assert!(
        route.length > lower_bound,
        "the detour must actually be longer than the straight line"
    );
}
