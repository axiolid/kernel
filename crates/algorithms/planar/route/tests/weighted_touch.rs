//! Cost regions touching a wall within rounding (#198).
//!
//! A costed footprint clipped to the free region lands on the wall only up
//! to rounding: a vertex a few nanometres inside leaves a sliver along the
//! wall a walk could slip through at factor 1, one a few nanometres beyond
//! makes the edge cross the wall. Either way the region touches the wall
//! and is taken as touching it.
//!
//! The room `[0, 10] x [0, 4]` and a factor-2 corridor `[0, 10] x [0, 1.5]`
//! along its bottom wall, turned 30 degrees. From the room's corner at
//! `(10, 0)` to the one at `(0, 0)`: along the wall only the corridor's side
//! is free, so the cheapest walk leaves the corridor, crosses the room
//! above it and comes back, refracting at the corridor's edge: each leg
//! through the corridor runs `sqrt 0.75` along it while rising 1.5, for
//! `2 * 2 sqrt 3 + 10 - 2 sqrt 0.75 = 10 + 3 sqrt 3`. A sliver along the
//! wall would allow 10.

use axiolid_core::{Frame2, Point2, Tolerance, Vec2};
use axiolid_overlay::{overlay, FillRule, OverlayInput, OverlayOperation, Polygon, Ring};
use axiolid_route::{weighted_distance_map, CostRegion, MapError};

const SPACING: f64 = 0.25;

/// Corner to corner: `10 + 2 * 1.5 * sqrt 3`.
const EXACT: f64 = 10.0 + 3.0 * 1.732_050_807_568_877_2;

fn turn(x: f64, y: f64) -> Point2 {
    let (s, c) = 30f64.to_radians().sin_cos();
    Point2::new(c * x - s * y, s * x + c * y)
}

fn turned(x0: f64, y0: f64, x1: f64, y1: f64) -> Polygon {
    Polygon {
        outer: Ring {
            points: vec![turn(x0, y0), turn(x1, y0), turn(x1, y1), turn(x0, y1)],
        },
        holes: Vec::new(),
    }
}

fn room() -> [Polygon; 1] {
    [turned(0.0, 0.0, 10.0, 4.0)]
}

/// The cost from corner to corner, bracketed.
fn corner_to_corner(corridor: Polygon) -> Result<(f64, f64), MapError> {
    let map = weighted_distance_map(
        &room(),
        &[],
        &[turn(0.0, 0.0)],
        &[CostRegion::new(corridor, 2.0)],
        SPACING,
    )?;
    let reach = map.nearest(turn(10.0, 0.0)).unwrap().unwrap();
    Ok((reach.cost.lower, reach.cost.upper))
}

fn holds_exact((lower, upper): (f64, f64)) {
    assert!(
        lower <= EXACT + 1e-9 && EXACT - 1e-9 <= upper,
        "[{lower}, {upper}] misses {EXACT}"
    );
    // As wide as the same map unturned: first order in the spacing.
    assert!(upper - lower <= 2.0 * SPACING, "[{lower}, {upper}]");
}

#[test]
fn a_corridor_clipped_to_the_room_along_a_turned_wall_builds_a_map() {
    // The footprint overshoots the room on three sides; clipped to it, its
    // new corners on the side walls are rounded.
    let frame = Frame2 {
        origin: Point2::new(0.0, 0.0),
        x: Vec2::new(1.0, 0.0),
        y: Vec2::new(0.0, 1.0),
    };
    let clipped = overlay(
        &OverlayInput {
            frame,
            polygons: vec![turned(-1.0, -1.0, 11.0, 1.5)],
        },
        &OverlayInput {
            frame,
            polygons: room().to_vec(),
        },
        OverlayOperation::Intersection,
        FillRule::NonZero,
        Tolerance::new(1e-9, 1e-9).unwrap(),
    )
    .unwrap();
    assert_eq!(clipped.polygons.len(), 1);
    holds_exact(corner_to_corner(clipped.polygons[0].clone()).expect("a map"));
}

#[test]
fn a_corridor_a_hair_inside_or_beyond_its_walls_touches_them() {
    // Each corner off the walls by a few nanometres, inside or beyond,
    // in every combination: none leaves a sliver, none is refused.
    let (s, c) = 30f64.to_radians().sin_cos();
    let along = Vec2::new(c, s);
    let up = Vec2::new(-s, c);
    for pattern in 0..16u32 {
        let off = |bit: u32| if pattern >> bit & 1 == 1 { 3e-9 } else { -3e-9 };
        let corridor = Polygon {
            outer: Ring {
                points: vec![
                    turn(0.0, 0.0) + up * off(0) + along * off(1),
                    turn(10.0, 0.0) + up * off(2) - along * off(3),
                    turn(10.0, 1.5) - along * off(1),
                    turn(0.0, 1.5) + along * off(3),
                ],
            },
            holes: Vec::new(),
        };
        let bracket = corner_to_corner(corridor)
            .unwrap_or_else(|e| panic!("pattern {pattern:04b} refused: {e:?}"));
        eprintln!("pattern {pattern:04b}: {bracket:?}");
        holds_exact(bracket);
    }
}

#[test]
fn a_cost_edge_crossing_a_wall_far_from_its_ends_is_still_refused() {
    // A corridor poking a metre through the far wall is no rounding.
    let corridor = turned(2.0, -1.0, 8.0, 1.5);
    assert!(matches!(
        corner_to_corner(corridor),
        Err(MapError::CostCrossing { index: 0 })
    ));
}

/// A factor-2 square standing on a wall turned by the 3-4-5 angle, cut to
/// the room by the exact overlay: its corners on the wall lie a rounding
/// step beyond it, so no wall runs exactly along the edge between them.
/// Along that edge only the region's side is free, as along the wall
/// itself; counted at the cheaper side, the walk under the square cost 1
/// a metre and the lower bound fell to 9.28 (#222).
#[test]
fn a_cost_edge_just_beyond_a_turned_wall_counts_its_free_side() {
    let p = Point2::new;
    let room = [Polygon {
        outer: Ring {
            points: vec![p(-2.4, 3.2), p(0.0, 0.0), p(8.0, 6.0), p(5.6, 9.2)],
        },
        holes: Vec::new(),
    }];
    let square = CostRegion::new(
        Polygon {
            outer: Ring {
                points: vec![
                    p(1.4000000000000004, 4.800000000000001),
                    p(3.2, 2.4),
                    p(4.800000000000001, 3.6),
                    p(3.000000000000001, 6.0),
                ],
            },
            holes: Vec::new(),
        },
        2.0,
    );
    let map = weighted_distance_map(
        &room,
        &[],
        &[p(7.000000000000001, 6.5)],
        &[square],
        0.009765625,
    )
    .expect("a map");
    let reach = map
        .nearest(p(0.8 * 0.5 - 0.6, 0.6 * 0.5 + 0.8))
        .unwrap()
        .unwrap();
    let exact = 2.0 * 3.5f64.hypot(2.0) + 2.0;
    assert!(
        reach.cost.lower <= exact + 1e-9 && exact <= reach.cost.upper + 1e-9,
        "{:?} misses {exact}",
        reach.cost
    );
    assert!(
        reach.cost.upper - reach.cost.lower <= 0.05,
        "{:?}",
        reach.cost
    );
}
