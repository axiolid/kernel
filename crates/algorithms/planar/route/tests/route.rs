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

/// Two 4 x 4 rooms joined by a 1 m corridor (#187).
fn h_rooms() -> Vec<Polygon> {
    vec![Polygon {
        outer: ring(&[
            (0.0, 0.0),
            (4.0, 0.0),
            (4.0, 1.0),
            (5.0, 1.0),
            (5.0, 0.0),
            (9.0, 0.0),
            (9.0, 4.0),
            (5.0, 4.0),
            (5.0, 2.0),
            (4.0, 2.0),
            (4.0, 4.0),
            (0.0, 4.0),
        ]),
        holes: vec![],
    }]
}

#[test]
fn a_cut_corridor_disconnects_even_where_walls_line_up() {
    // The rooms' top walls lie on one line, y = 4, with a gap outside the
    // region between x = 4 and 5. With the corridor cut, no route runs
    // along the walls across that gap.
    let p = Point2::new;
    let barriers = vec![vec![p(4.5, 0.9), p(4.5, 2.1)]];
    assert_eq!(
        shortest_path(&h_rooms(), &barriers, p(1.0, 3.0), p(8.0, 3.0)),
        Ok(Err(Unreachable::DisconnectedComponents))
    );
    // Without the cut, the route through the corridor.
    let route = shortest_path(&h_rooms(), &[], p(1.0, 3.0), p(8.0, 3.0))
        .unwrap()
        .unwrap();
    assert_eq!(
        route.polyline,
        vec![p(1.0, 3.0), p(4.0, 2.0), p(5.0, 2.0), p(8.0, 3.0)]
    );
    assert!((route.length - (2.0 * 10f64.sqrt() + 1.0)).abs() < 1e-12);
}

/// Whether `q` lies in the closed region: inside an outer ring and outside
/// every hole, or within a hair of a boundary.
fn in_closed(region: &[Polygon], q: Point2) -> bool {
    let inside = |r: &Ring| {
        let mut c = false;
        let n = r.points.len();
        for i in 0..n {
            let (a, b) = (r.points[i], r.points[(i + 1) % n]);
            if (a.y > q.y) != (b.y > q.y) && q.x < (b.x - a.x) * (q.y - a.y) / (b.y - a.y) + a.x {
                c = !c;
            }
        }
        c
    };
    let near = |r: &Ring| {
        let n = r.points.len();
        (0..n).any(|i| {
            let (a, b) = (r.points[i], r.points[(i + 1) % n]);
            let d = b - a;
            let t = ((q - a).dot(d) / d.dot(d)).clamp(0.0, 1.0);
            (a + d * t - q).length() <= 1e-9
        })
    };
    region.iter().any(|poly| {
        let rings = std::iter::once(&poly.outer).chain(&poly.holes);
        rings.clone().any(&near) || (inside(&poly.outer) && !poly.holes.iter().any(&inside))
    })
}

#[test]
fn every_route_segment_lies_in_the_closed_region() {
    let p = Point2::new;
    // Regions whose walls line up across gaps, run along each other and
    // around holes: the H with and without its cut, a comb, and a room
    // with a hole whose side lies on the line of a wall.
    let comb = vec![Polygon {
        outer: ring(&[
            (0.0, 0.0),
            (7.0, 0.0),
            (7.0, 3.0),
            (6.0, 3.0),
            (6.0, 1.0),
            (5.0, 1.0),
            (5.0, 3.0),
            (4.0, 3.0),
            (4.0, 1.0),
            (3.0, 1.0),
            (3.0, 3.0),
            (2.0, 3.0),
            (2.0, 1.0),
            (1.0, 1.0),
            (1.0, 3.0),
            (0.0, 3.0),
        ]),
        holes: vec![],
    }];
    let holed = vec![Polygon {
        outer: ring(&[(0.0, 0.0), (6.0, 0.0), (6.0, 4.0), (0.0, 4.0)]),
        holes: vec![ring(&[
            (2.0, 1.0),
            (4.0, 1.0),
            (4.0, 4.0 - 1e-3),
            (2.0, 4.0 - 1e-3),
        ])],
    }];
    let cases: Vec<(Vec<Polygon>, Vec<Vec<Point2>>)> = vec![
        (h_rooms(), vec![]),
        (h_rooms(), vec![vec![p(4.5, 0.9), p(4.5, 2.1)]]),
        (comb, vec![]),
        (holed, vec![]),
    ];
    let mut routes = 0;
    for (region, barriers) in &cases {
        // Deterministic start and goal pairs over the region's box.
        let mut state = 0x2545_f491_4f6c_dd1d_u64;
        let mut next = || {
            state ^= state << 13;
            state ^= state >> 7;
            state ^= state << 17;
            (state >> 11) as f64 / (1u64 << 53) as f64
        };
        for _ in 0..60 {
            let (s, g) = (p(9.0 * next(), 4.0 * next()), p(9.0 * next(), 4.0 * next()));
            let Ok(Ok(route)) = shortest_path(region, barriers, s, g) else {
                continue;
            };
            routes += 1;
            for w in route.polyline.windows(2) {
                for k in 0..=200 {
                    let q = w[0] + (w[1] - w[0]) * (k as f64 / 200.0);
                    assert!(
                        in_closed(region, q),
                        "{q:?} on {:?} leaves the region",
                        route.polyline
                    );
                }
            }
        }
    }
    assert!(routes > 60, "{routes} routes");
}

fn length_of(
    region: &[Polygon],
    barriers: &[Vec<Point2>],
    from: (f64, f64),
    to: (f64, f64),
) -> f64 {
    let (a, b) = (Point2::new(from.0, from.1), Point2::new(to.0, to.1));
    let route = shortest_path(region, barriers, a, b).unwrap().unwrap();
    // Whatever the length, the polyline must be it.
    let walked: f64 = route
        .polyline
        .windows(2)
        .map(|w| (w[1] - w[0]).length())
        .sum();
    assert!((walked - route.length).abs() < 1e-9);
    route.length
}

#[test]
fn a_barrier_touching_a_wall_is_not_squeezed_past() {
    // #189: the barrier's foot touches the bottom wall. Along the wall
    // through the foot is a gap of zero width; the way is over the top.
    let barrier = vec![vec![Point2::new(5.0, 0.0), Point2::new(5.0, 8.0)]];
    let length = length_of(&room(), &barrier, (0.0, 0.0), (10.0, 0.0));
    assert!((length - 2.0 * 89f64.sqrt()).abs() < 1e-9, "{length}");
    // From just above the wall, too.
    let length = length_of(&room(), &barrier, (1.0, 1.0), (9.0, 1.0));
    assert!(
        (length - 2.0 * (16.0f64 + 49.0).sqrt()).abs() < 1e-9,
        "{length}"
    );
}

#[test]
fn holes_touching_at_a_corner_are_not_squeezed_between() {
    // Two pillars meet corner to corner at (5, 5); the straight line from
    // (3, 7) to (7, 3) runs through that point. Round either pillar is
    // 2 sqrt(5) + 6.
    let region = vec![Polygon {
        outer: ring(&[(0.0, 0.0), (10.0, 0.0), (10.0, 10.0), (0.0, 10.0)]),
        holes: vec![
            ring(&[(2.0, 2.0), (2.0, 5.0), (5.0, 5.0), (5.0, 2.0)]),
            ring(&[(5.0, 5.0), (5.0, 8.0), (8.0, 8.0), (8.0, 5.0)]),
        ],
    }];
    let length = length_of(&region, &[], (3.0, 7.0), (7.0, 3.0));
    assert!(
        (length - (2.0 * 5f64.sqrt() + 6.0)).abs() < 1e-9,
        "{length}"
    );
    // The other diagonal runs through the same point, between the two
    // open quarters: round a pillar again.
    let length = length_of(&region, &[], (6.0, 4.0), (4.0, 6.0));
    assert!(
        (length - (2.0 * 5f64.sqrt() + 6.0)).abs() < 1e-9,
        "{length}"
    );
}

#[test]
fn a_bent_barrier_is_one_wall() {
    // A wall drawn in two collinear pieces: its middle vertex is no gap.
    let region = vec![Polygon {
        outer: ring(&[(0.0, 0.0), (6.0, 0.0), (6.0, 10.0), (0.0, 10.0)]),
        holes: Vec::new(),
    }];
    let wall = vec![vec![
        Point2::new(5.0, 0.5),
        Point2::new(5.0, 4.0),
        Point2::new(5.0, 8.5),
    ]];
    let length = length_of(&region, &wall, (0.0, 5.0), (6.0, 4.0));
    let expected = 45.25f64.sqrt() + 13.25f64.sqrt();
    assert!((length - expected).abs() < 1e-9, "{length} vs {expected}");
    // Bent at a right angle, the corner is no gap either.
    let corner = vec![vec![
        Point2::new(2.0, 2.0),
        Point2::new(2.0, 8.0),
        Point2::new(8.0, 8.0),
    ]];
    // From outside the corner to inside it, the straight line runs through
    // the bend; round either end is sqrt(50) + sqrt(26).
    let length = length_of(&room(), &corner, (1.0, 9.0), (3.0, 7.0));
    let expected = 50f64.sqrt() + 26f64.sqrt();
    assert!((length - expected).abs() < 1e-9, "{length} vs {expected}");
}

#[test]
fn a_barrier_along_a_wall_leaves_the_wall_walkable() {
    let barrier = vec![vec![Point2::new(2.0, 0.0), Point2::new(4.0, 0.0)]];
    let length = length_of(&room(), &barrier, (0.0, 0.0), (10.0, 0.0));
    assert!((length - 10.0).abs() < 1e-12, "{length}");
}

#[test]
fn two_rooms_sharing_a_wall_are_walkable_along_it() {
    // Two polygons share the edge y = 5: a wall between rooms, walkable
    // on either side, closed to crossing.
    let region = vec![
        Polygon {
            outer: ring(&[(0.0, 0.0), (10.0, 0.0), (10.0, 5.0), (0.0, 5.0)]),
            holes: Vec::new(),
        },
        Polygon {
            outer: ring(&[(0.0, 5.0), (10.0, 5.0), (10.0, 10.0), (0.0, 10.0)]),
            holes: Vec::new(),
        },
    ];
    let length = length_of(&region, &[], (1.0, 5.0), (9.0, 5.0));
    assert!((length - 8.0).abs() < 1e-12, "{length}");
    let length = length_of(&region, &[], (0.0, 5.0), (3.0, 9.0));
    assert!((length - 5.0).abs() < 1e-12, "{length}");
}
