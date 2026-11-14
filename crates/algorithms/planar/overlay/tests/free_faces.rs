//! Bounded free faces between polygons and open segments, closing gaps with
//! reported bridges (#490, ADR 0109).
//!
//! `scripts/probe_free_faces_mutants.py` lists the faults these tests must
//! catch.

use axiolid_contracts::{CancellationToken, ExecutionOptions};
use axiolid_core::{Point2, Tolerance};
use axiolid_overlay::{
    free_faces, free_faces_within, BarrierDefect, FaceLocation, FreeFace, FreeFaces,
    FreeFacesError, FreeRing, FreeSpaceInput, OverlayError, PieceSource, Polygon, Ring,
};

fn t() -> Tolerance {
    Tolerance::new(1e-9, 1e-9).unwrap()
}

fn ring(points: &[(f64, f64)]) -> Ring {
    Ring {
        points: points.iter().map(|&(x, y)| Point2::new(x, y)).collect(),
    }
}

fn polygon(points: &[(f64, f64)]) -> Polygon {
    Polygon {
        outer: ring(points),
        holes: Vec::new(),
    }
}

/// An axis-aligned rectangle, counter-clockwise.
fn rect(x0: f64, y0: f64, x1: f64, y1: f64) -> Polygon {
    polygon(&[(x0, y0), (x1, y0), (x1, y1), (x0, y1)])
}

fn point(x: f64, y: f64) -> Point2 {
    Point2::new(x, y)
}

fn barrier(a: (f64, f64), b: (f64, f64)) -> [Point2; 2] {
    [point(a.0, a.1), point(b.0, b.1)]
}

fn input(obstacles: Vec<Polygon>, barriers: Vec<[Point2; 2]>) -> FreeSpaceInput {
    FreeSpaceInput {
        obstacles,
        barriers,
    }
}

const GAP: f64 = 1e-4;
/// 2^-17, about 7.6e-6.
const D: f64 = 1.0 / 131_072.0;

/// Four wall pieces round the rectangle `[x0, x0 + w] x [y0, y0 + h]`,
/// `k` thick, with mitred ends each pulled back along its wall by the entry
/// of `pulls`: south (west, east end), west (south, north), north (east,
/// west) and east (north, south). Where two pulls add up to something
/// positive the walls miss each other at that corner.
fn mitred_room(x0: f64, y0: f64, (w, h): (f64, f64), k: f64, p: [f64; 8]) -> Vec<Polygon> {
    vec![
        polygon(&[
            (x0 + p[0], y0),
            (x0 + w - p[1], y0),
            (x0 + w + k - p[1], y0 - k),
            (x0 - k + p[0], y0 - k),
        ]),
        polygon(&[
            (x0, y0 + h - p[3]),
            (x0, y0 + p[2]),
            (x0 - k, y0 - k + p[2]),
            (x0 - k, y0 + h + k - p[3]),
        ]),
        polygon(&[
            (x0 + w - p[4], y0 + h),
            (x0 + p[5], y0 + h),
            (x0 - k + p[5], y0 + h + k),
            (x0 + w + k - p[4], y0 + h + k),
        ]),
        polygon(&[
            (x0 + w, y0 + p[7]),
            (x0 + w, y0 + h - p[6]),
            (x0 + w + k, y0 + h + k - p[6]),
            (x0 + w + k, y0 - k + p[7]),
        ]),
    ]
}

/// The room `[0, 4] x [0, 3]`, 0.25 thick walls, every end pulled back by
/// `d`: the mitre faces at a corner are `d` times root two apart.
fn mitred_walls(d: f64) -> Vec<Polygon> {
    mitred_room(0.0, 0.0, (4.0, 3.0), 0.25, [d; 8])
}

fn area(points: &[Point2]) -> f64 {
    let n = points.len();
    (0..n)
        .map(|i| points[i].x * points[(i + 1) % n].y - points[(i + 1) % n].x * points[i].y)
        .sum::<f64>()
        / 2.0
}

fn is_obstacle(sources: &[PieceSource]) -> bool {
    matches!(sources, [PieceSource::Obstacle { .. }])
}

fn is_bridge(sources: &[PieceSource]) -> bool {
    matches!(sources, [PieceSource::Bridge { .. }])
}

fn count_edges(face: &FreeFace, kind: fn(&[PieceSource]) -> bool) -> usize {
    face.outer
        .sources
        .iter()
        .filter(|sources| kind(sources))
        .count()
}

/// Walls round `[0, 4] x [0, 3]` as four rectangles that overlap at the
/// corners: a closed room with no gap at all.
fn closed_room() -> Vec<Polygon> {
    vec![
        rect(-0.25, -0.25, 4.25, 0.0),
        rect(-0.25, 3.0, 4.25, 3.25),
        rect(-0.25, 0.0, 0.0, 3.0),
        rect(4.0, 0.0, 4.25, 3.0),
    ]
}

/// A wall round the hole `room`: the free face is the hole.
fn walled(room: &[(f64, f64)]) -> Polygon {
    Polygon {
        outer: ring(&[(-10.0, -10.0), (10.0, -10.0), (10.0, 10.0), (-10.0, 10.0)]),
        holes: vec![ring(room)],
    }
}

fn faces_of(obstacles: Vec<Polygon>, barriers: Vec<[Point2; 2]>, gap: f64) -> FreeFaces {
    free_faces(&input(obstacles, barriers), gap, t()).unwrap()
}

// --- the issue's cases -----------------------------------------------------

#[test]
fn four_mitred_walls_with_gaps_close_into_one_room() {
    let faces = faces_of(mitred_walls(D), Vec::new(), GAP);
    assert_eq!(faces.faces().len(), 1);
    let room = &faces.faces()[0];
    assert_eq!(room.outer.points.len(), 8);
    assert!(room.holes.is_empty());
    assert!((area(&room.outer.points) - 12.0).abs() < 4.0 * D);
    assert_eq!(faces.closed_gaps().len(), 4);
    assert_eq!(faces.bridges().len(), 8);
    assert!(faces.bridges().iter().all(|bridge| bridge.length < GAP));
    assert_eq!(count_edges(room, is_bridge), 4);
    assert_eq!(count_edges(room, is_obstacle), 4);
    assert!(faces.evidence().candidate_pairs > 0);
    assert_eq!(faces.evidence().collapsed_faces, 0);
}

#[test]
fn mitred_walls_pulled_back_by_more_than_the_gap_leave_the_room_open() {
    let faces = faces_of(mitred_walls(1.0 / 1024.0), Vec::new(), GAP);
    assert!(faces.faces().is_empty());
    assert!(faces.bridges().is_empty());
    assert!(faces.closed_gaps().is_empty());
}

#[test]
fn three_walls_and_a_barrier_close_a_room_whose_fourth_side_is_the_barrier() {
    let walls = vec![
        rect(-0.25, -0.25, 4.25, 0.0),
        rect(-0.25, 3.0, 4.25, 3.25),
        rect(4.0, 0.0, 4.25, 3.0),
    ];
    let faces = faces_of(walls, vec![barrier((0.0, 0.0), (0.0, 3.0))], GAP);
    assert_eq!(faces.faces().len(), 1);
    let room = &faces.faces()[0];
    assert_eq!(
        room.outer.points,
        vec![
            point(0.0, 0.0),
            point(4.0, 0.0),
            point(4.0, 3.0),
            point(0.0, 3.0)
        ]
    );
    assert_eq!(count_edges(room, is_obstacle), 3);
    assert_eq!(
        room.outer.sources[3],
        vec![PieceSource::Barrier { barrier: 0 }]
    );
    assert!(faces.bridges().is_empty());
}

#[test]
fn a_dangling_barrier_is_in_no_ring_and_its_points_are_inside() {
    let faces = faces_of(closed_room(), vec![barrier((0.0, 1.5), (2.0, 1.5))], GAP);
    assert_eq!(faces.faces().len(), 1);
    let room = &faces.faces()[0];
    assert_eq!(room.outer.points.len(), 4);
    assert!(room
        .outer
        .sources
        .iter()
        .flatten()
        .all(|source| !matches!(source, PieceSource::Barrier { .. })));
    assert_eq!(
        faces.locate(point(1.0, 1.5)).unwrap(),
        FaceLocation::Inside(0)
    );
    assert_eq!(
        faces.locate(point(2.0, 1.5)).unwrap(),
        FaceLocation::Inside(0)
    );
}

#[test]
fn an_island_is_a_hole_with_its_own_sources() {
    let mut obstacles = closed_room();
    obstacles.push(rect(1.0, 1.0, 2.0, 2.0));
    let faces = faces_of(obstacles, Vec::new(), GAP);
    assert_eq!(faces.faces().len(), 1);
    let room = &faces.faces()[0];
    assert_eq!(room.holes.len(), 1);
    let island = &room.holes[0];
    assert_eq!(island.points.len(), 4);
    assert!(island.sources.iter().all(|sources| matches!(
        sources[..],
        [PieceSource::Obstacle {
            polygon: 4,
            ring: 0,
            ..
        }]
    )));
    assert!(area(&island.points) < 0.0, "a hole is clockwise");
}

#[test]
fn locate_names_every_variant() {
    let faces = faces_of(mitred_walls(D), Vec::new(), GAP);
    assert_eq!(
        faces.locate(point(2.0, 1.5)).unwrap(),
        FaceLocation::Inside(0)
    );
    let room = &faces.faces()[0];
    let middle = |edge: usize| {
        let (a, b) = (
            room.outer.points[edge],
            room.outer.points[(edge + 1) % room.outer.points.len()],
        );
        point((a.x + b.x) / 2.0, (a.y + b.y) / 2.0)
    };
    // On a wall edge of the room, and on a bridge edge: the sources say which.
    for kind in [is_obstacle as fn(&[PieceSource]) -> bool, is_bridge] {
        let edge = room.outer.sources.iter().position(|s| kind(s)).unwrap();
        assert_eq!(
            faces.locate(middle(edge)).unwrap(),
            FaceLocation::OnBoundary {
                face: 0,
                ring: 0,
                edge
            }
        );
    }
    // In a closed gap: the centre of its four corners.
    let gap = &faces.closed_gaps()[0].outer.points;
    let centre = point(
        gap.iter().map(|p| p.x).sum::<f64>() / 4.0,
        gap.iter().map(|p| p.y).sum::<f64>() / 4.0,
    );
    assert!(matches!(
        faces.locate(centre).unwrap(),
        FaceLocation::InClosedGap(_)
    ));
    assert_eq!(
        faces.locate(point(2.0, -0.1)).unwrap(),
        FaceLocation::InObstacle(0)
    );
    assert_eq!(
        faces.locate(point(-0.1, 1.5)).unwrap(),
        FaceLocation::InObstacle(1)
    );
    assert_eq!(
        faces.locate(point(-5.0, 1.5)).unwrap(),
        FaceLocation::Unbounded
    );
    assert_eq!(
        faces.locate(point(2.0, 100.0)).unwrap(),
        FaceLocation::Unbounded
    );
    assert_eq!(
        faces.locate(point(f64::NAN, 0.0)),
        Err(FreeFacesError::NonFinitePoint)
    );
}

#[test]
fn locate_puts_a_point_on_a_hole_boundary_on_the_boundary_and_inside_an_island_in_the_island() {
    let mut obstacles = closed_room();
    obstacles.push(rect(1.0, 1.0, 2.0, 2.0));
    let faces = faces_of(obstacles, Vec::new(), GAP);
    assert!(matches!(
        faces.locate(point(1.5, 1.0)).unwrap(),
        FaceLocation::OnBoundary { ring: 1, .. }
    ));
    assert_eq!(
        faces.locate(point(1.5, 1.5)).unwrap(),
        FaceLocation::InObstacle(4)
    );
    assert_eq!(
        faces.locate(point(3.0, 1.5)).unwrap(),
        FaceLocation::Inside(0)
    );
}

// --- refusals --------------------------------------------------------------

#[test]
fn a_gap_that_is_not_finite_or_is_negative_is_refused() {
    let bad = input(closed_room(), Vec::new());
    for gap in [f64::NAN, f64::INFINITY, f64::NEG_INFINITY, -1e-9, -0.1] {
        assert_eq!(
            free_faces(&bad, gap, t()).unwrap_err(),
            FreeFacesError::InvalidGap
        );
    }
    assert!(free_faces(&bad, 0.0, t()).is_ok());
    assert!(free_faces(&bad, -0.0, t()).is_ok());
}

#[test]
fn every_obstacle_defect_is_named_with_the_index_of_the_polygon() {
    let reason = |bad: Polygon| -> (usize, OverlayError) {
        let obstacles = vec![rect(0.0, 0.0, 1.0, 1.0), rect(2.0, 0.0, 3.0, 1.0), bad];
        match free_faces(&input(obstacles, Vec::new()), GAP, t()).unwrap_err() {
            FreeFacesError::InvalidObstacle { polygon, reason } => (polygon, reason),
            other => panic!("{other:?}"),
        }
    };
    assert_eq!(
        reason(polygon(&[(0.0, 0.0), (1.0, 0.0)])),
        (2, OverlayError::RingTooShort)
    );
    assert_eq!(
        reason(polygon(&[(0.0, 0.0), (1.0, f64::NAN), (1.0, 1.0)])),
        (2, OverlayError::NonFinitePoint)
    );
    assert_eq!(
        reason(polygon(&[(0.0, 0.0), (1.0, 0.0), (1.0, 0.0), (1.0, 1.0)])),
        (2, OverlayError::RepeatedVertex)
    );
    assert_eq!(
        reason(polygon(&[(0.0, 0.0), (2.0, 2.0), (2.0, 0.0), (0.0, 2.0)])),
        (2, OverlayError::SelfIntersection)
    );
    // On one line: no wall at all, which a room must not leak through.
    assert_eq!(
        reason(polygon(&[(0.0, 0.0), (1.0, 0.0), (2.0, 0.0)])),
        (2, OverlayError::ZeroArea)
    );
    let square = [(0.0, 0.0), (4.0, 0.0), (4.0, 4.0), (0.0, 4.0)];
    let with_holes = |holes: &[&[(f64, f64)]]| Polygon {
        outer: ring(&square),
        holes: holes.iter().map(|hole| ring(hole)).collect(),
    };
    assert_eq!(
        reason(with_holes(&[&[(1.0, 1.0), (2.0, 1.0), (3.0, 1.0)]])),
        (2, OverlayError::ZeroArea)
    );
    assert_eq!(
        reason(with_holes(&[&[(5.0, 5.0), (6.0, 5.0), (6.0, 6.0)]])),
        (2, OverlayError::HoleOutsideOuter)
    );
    assert_eq!(
        reason(with_holes(&[
            &[(1.0, 1.0), (3.0, 1.0), (3.0, 3.0), (1.0, 3.0)],
            &[(2.0, 2.0), (2.5, 2.0), (2.5, 2.5), (2.0, 2.5)],
        ])),
        (2, OverlayError::HolesOverlap)
    );
}

#[test]
fn the_first_defect_wins_in_the_order_gap_obstacles_barriers() {
    let bad_obstacle = polygon(&[(0.0, 0.0), (1.0, 0.0)]);
    let bad_barrier = barrier((0.0, 0.0), (0.0, 0.0));
    let both = input(vec![bad_obstacle.clone()], vec![bad_barrier]);
    assert_eq!(
        free_faces(&both, f64::NAN, t()).unwrap_err(),
        FreeFacesError::InvalidGap
    );
    assert_eq!(
        free_faces(&both, GAP, t()).unwrap_err(),
        FreeFacesError::InvalidObstacle {
            polygon: 0,
            reason: OverlayError::RingTooShort
        }
    );
    // Of two bad obstacles the lower index.
    let two = input(vec![bad_obstacle.clone(), bad_obstacle], Vec::new());
    assert!(matches!(
        free_faces(&two, GAP, t()).unwrap_err(),
        FreeFacesError::InvalidObstacle { polygon: 0, .. }
    ));
}

#[test]
fn every_barrier_defect_is_named_with_the_index_of_the_barrier() {
    let tolerance = Tolerance::new(1e-3, 1e-3).unwrap();
    let with = |bad: [Point2; 2]| {
        let barriers = vec![
            barrier((0.0, 0.0), (1.0, 0.0)),
            barrier((0.0, 1.0), (1.0, 1.0)),
            bad,
        ];
        free_faces(&input(Vec::new(), barriers), GAP, tolerance).unwrap_err()
    };
    for coordinate in [f64::NAN, f64::INFINITY] {
        assert_eq!(
            with([point(0.0, 0.0), point(coordinate, 1.0)]),
            FreeFacesError::InvalidBarrier {
                barrier: 2,
                reason: BarrierDefect::NonFinitePoint
            }
        );
    }
    // As long as the tolerance is too short; longer is fine.
    for ends in [
        [point(0.0, 0.0), point(0.0, 0.0)],
        [point(0.0, 0.0), point(1e-3, 0.0)],
    ] {
        assert_eq!(
            with(ends),
            FreeFacesError::InvalidBarrier {
                barrier: 2,
                reason: BarrierDefect::TooShort
            }
        );
    }
    let barriers = vec![barrier((0.0, 0.0), (1.5e-3, 0.0))];
    assert!(free_faces(&input(Vec::new(), barriers), GAP, tolerance).is_ok());
}

#[test]
fn nothing_gives_nothing() {
    let faces = faces_of(Vec::new(), Vec::new(), GAP);
    assert!(faces.faces().is_empty() && faces.bridges().is_empty());
    assert_eq!(
        faces.locate(point(0.0, 0.0)).unwrap(),
        FaceLocation::Unbounded
    );
}

#[test]
fn errors_say_what_is_wrong() {
    let text = FreeFacesError::InvalidObstacle {
        polygon: 3,
        reason: OverlayError::ZeroArea,
    }
    .to_string();
    assert!(text.contains('3') && text.contains("no area"), "{text}");
    assert!(FreeFacesError::InvalidGap
        .to_string()
        .contains("closing width"));
    let source: &dyn std::error::Error = &FreeFacesError::Cancelled;
    assert!(source.source().is_none());
}

// --- budget and cancellation -----------------------------------------------

/// `n` staggered parallel barriers 1e-3 apart: every vertex is within the
/// gap of every barrier, but only the neighbours are unshadowed, so the
/// search keeps a foot for each of about `n` squared candidates.
fn comb(n: usize) -> FreeSpaceInput {
    input(
        Vec::new(),
        (0..n)
            .map(|i| {
                let y = i as f64 * 1e-3;
                if i % 2 == 0 {
                    barrier((0.0, y), (1.0, y))
                } else {
                    barrier((-0.5, y), (1.5, y))
                }
            })
            .collect(),
    )
}

#[test]
fn the_memory_budget_bounds_the_search_for_bridges() {
    let n = 60;
    // Without a budget the same input is searched to the end.
    let all = free_faces_within(&comb(n), 10.0, &ExecutionOptions::new(t())).unwrap();
    let kept = all.bridges().len();
    assert!(kept > 10);
    // The search keeps a foot for each of the pairs it tests, far more than
    // the bridges it keeps; a budget for the inputs, those bridges and their
    // feet, and a little more, refuses it.
    let options = ExecutionOptions::new(t()).with_memory_budget(256 * (n + 3 * kept));
    assert_eq!(
        free_faces_within(&comb(n), 10.0, &options).unwrap_err(),
        FreeFacesError::BudgetExceeded { resource: "memory" }
    );
    // Bridges between ends keep no foot, and count by themselves: `m` rows
    // of two barriers that nearly meet.
    let m = 20;
    let rows = input(
        Vec::new(),
        (0..m)
            .flat_map(|i| {
                let y = i as f64;
                [
                    barrier((0.0, y), (1.0, y)),
                    barrier((1.0 + 1e-3, y), (2.0, y)),
                ]
            })
            .collect(),
    );
    let few = ExecutionOptions::new(t()).with_memory_budget(2 * m * 256 + 3 * 256);
    assert_eq!(
        free_faces_within(&rows, 2e-3, &few).unwrap_err(),
        FreeFacesError::BudgetExceeded { resource: "memory" }
    );
    let all = free_faces_within(&rows, 2e-3, &ExecutionOptions::new(t())).unwrap();
    assert_eq!(all.bridges().len(), m);
    // A budget for everything the search keeps lets it finish.
    let ample = ExecutionOptions::new(t()).with_memory_budget(256 * (n + 4 * n * n));
    assert!(free_faces_within(&comb(n), 10.0, &ample).is_ok());
}

#[test]
fn a_tiny_memory_budget_is_refused_and_a_cancelled_token_stops() {
    let room = input(closed_room(), Vec::new());
    let small = ExecutionOptions::new(t()).with_memory_budget(1);
    assert_eq!(
        free_faces_within(&room, GAP, &small).unwrap_err(),
        FreeFacesError::BudgetExceeded { resource: "memory" }
    );
    let token = CancellationToken::new();
    token.cancel();
    let cancelled = ExecutionOptions::new(t()).with_cancellation(token);
    assert_eq!(
        free_faces_within(&room, GAP, &cancelled).unwrap_err(),
        FreeFacesError::Cancelled
    );
    assert!(free_faces_within(&room, GAP, &ExecutionOptions::new(t())).is_ok());
}

// --- adversarial -----------------------------------------------------------

#[test]
fn four_rooms_meeting_at_a_cross_with_gaps_at_the_centre_stay_four_faces() {
    let g = 2e-5;
    let walls = vec![walled(&[
        (-5.0, -5.0),
        (5.0, -5.0),
        (5.0, 5.0),
        (-5.0, 5.0),
    ])];
    let arms = vec![
        barrier((0.0, 5.0), (0.0, g)),
        barrier((5.0, 0.0), (g, 0.0)),
        barrier((0.0, -5.0), (0.0, -g)),
        barrier((-5.0, 0.0), (-g, 0.0)),
    ];
    let faces = faces_of(walls, arms, GAP);
    assert_eq!(faces.faces().len(), 4);
    // The bridges between the four tips cross at the centre and cut the
    // diamond into four triangles, each a closed gap.
    assert_eq!(faces.bridges().len(), 6);
    assert_eq!(faces.closed_gaps().len(), 4);
    for face in faces.faces() {
        assert!((area(&face.outer.points) - 25.0).abs() < 40.0 * GAP);
    }
}

#[test]
fn a_jog_shorter_than_the_gap_is_not_bridged_with_an_obtuse_corner() {
    let faces = faces_of(
        vec![walled(&[
            (0.0, 0.0),
            (2.0, 0.0),
            (2.0 + 3e-5, 3e-5),
            (4.0, 3e-5),
            (4.0, 3.0),
            (0.0, 3.0),
        ])],
        Vec::new(),
        GAP,
    );
    assert_eq!(faces.faces().len(), 1);
    assert_eq!(faces.faces()[0].outer.points.len(), 6);
    assert!(faces.bridges().is_empty());
    assert!(faces.closed_gaps().is_empty());
}

#[test]
fn a_jog_shorter_than_the_gap_along_a_shallow_line_is_not_cut_apart() {
    let faces = faces_of(
        vec![walled(&[
            (0.0, 0.0),
            (2.0, 0.0),
            (2.0 + 5.2e-5, 3e-5),
            (4.0, 3e-5),
            (4.0, 3.0),
            (0.0, 3.0),
        ])],
        Vec::new(),
        GAP,
    );
    assert_eq!(faces.faces().len(), 1);
    assert!(faces.bridges().is_empty());
    assert_eq!(faces.faces()[0].outer.points.len(), 6);
}

#[test]
fn a_pillar_corner_within_the_gap_of_a_wall_gets_one_dangling_bridge_and_stays_a_hole() {
    let pillar = polygon(&[(2.0, 5e-5), (2.5, 0.5), (2.0, 1.0), (1.5, 0.5)]);
    let faces = faces_of(
        vec![
            walled(&[(0.0, 0.0), (4.0, 0.0), (4.0, 3.0), (0.0, 3.0)]),
            pillar,
        ],
        Vec::new(),
        GAP,
    );
    assert_eq!(faces.bridges().len(), 1);
    assert_eq!(faces.faces().len(), 1);
    let room = &faces.faces()[0];
    assert_eq!(room.holes.len(), 1);
    assert_eq!(room.outer.points.len(), 4);
    let bridge = &faces.bridges()[0];
    assert_eq!(bridge.ends, [point(2.0, 0.0), point(2.0, 5e-5)]);
    assert_eq!(bridge.length, 5e-5);
    assert!(room
        .outer
        .sources
        .iter()
        .chain(&room.holes[0].sources)
        .flatten()
        .all(|source| !matches!(source, PieceSource::Bridge { .. })));
    // The bridge bounds nothing: its points are inside the room.
    assert_eq!(
        faces.locate(point(2.0, 2.5e-5)).unwrap(),
        FaceLocation::Inside(0)
    );
}

#[test]
fn a_pillar_face_parallel_to_a_wall_within_the_gap_is_a_closed_gap_and_the_pillar_joins_the_outer_ring(
) {
    let pillar = rect(1.0, 5e-5, 2.0, 1.0);
    let faces = faces_of(
        vec![
            walled(&[(0.0, 0.0), (4.0, 0.0), (4.0, 3.0), (0.0, 3.0)]),
            pillar,
        ],
        Vec::new(),
        GAP,
    );
    assert_eq!(faces.faces().len(), 1);
    let room = &faces.faces()[0];
    assert!(room.holes.is_empty());
    assert_eq!(faces.closed_gaps().len(), 1);
    assert_eq!(faces.bridges().len(), 2);
    let pillar_edges = room
        .outer
        .sources
        .iter()
        .filter(|sources| matches!(sources[..], [PieceSource::Obstacle { polygon: 1, .. }]))
        .count();
    assert_eq!(pillar_edges, 3);
    assert_eq!(count_edges(room, is_bridge), 2);
    // The strip under the pillar is a gap, not a face.
    assert!(matches!(
        faces.locate(point(1.5, 2.5e-5)).unwrap(),
        FaceLocation::InClosedGap(0)
    ));
}

#[test]
fn two_rooms_touching_at_a_corner_are_two_faces_sharing_a_bit_identical_point() {
    let wall = Polygon {
        outer: ring(&[(-1.0, -1.0), (9.0, -1.0), (9.0, 7.0), (-1.0, 7.0)]),
        holes: vec![
            ring(&[(0.0, 0.0), (4.0, 0.0), (4.0, 3.0), (0.0, 3.0)]),
            ring(&[(4.0, 3.0), (8.0, 3.0), (8.0, 6.0), (4.0, 6.0)]),
        ],
    };
    let faces = faces_of(vec![wall], Vec::new(), GAP);
    assert_eq!(faces.faces().len(), 2);
    for face in faces.faces() {
        assert!(face.outer.points.contains(&point(4.0, 3.0)));
    }
}

#[test]
fn an_island_touching_a_wall_at_a_vertex_is_one_ring_that_visits_it_twice() {
    let island = polygon(&[(2.0, 0.0), (3.0, 1.0), (2.0, 2.0), (1.0, 1.0)]);
    let faces = faces_of(
        vec![
            walled(&[(0.0, 0.0), (4.0, 0.0), (4.0, 3.0), (0.0, 3.0)]),
            island,
        ],
        Vec::new(),
        GAP,
    );
    assert_eq!(faces.faces().len(), 1);
    let room = &faces.faces()[0];
    assert!(room.holes.is_empty());
    let visits = room
        .outer
        .points
        .iter()
        .filter(|&&p| p == point(2.0, 0.0))
        .count();
    assert_eq!(visits, 2);
}

#[test]
fn a_courtyard_is_a_face_of_its_own_and_a_hole_of_the_room() {
    let courtyard = Polygon {
        outer: ring(&[(1.0, 1.0), (3.0, 1.0), (3.0, 2.0), (1.0, 2.0)]),
        holes: vec![ring(&[
            (1.25, 1.25),
            (2.75, 1.25),
            (2.75, 1.75),
            (1.25, 1.75),
        ])],
    };
    let faces = faces_of(
        vec![
            walled(&[(0.0, 0.0), (4.0, 0.0), (4.0, 3.0), (0.0, 3.0)]),
            courtyard,
        ],
        Vec::new(),
        GAP,
    );
    assert_eq!(faces.faces().len(), 2);
    let room = faces
        .faces()
        .iter()
        .position(|f| f.holes.len() == 1)
        .unwrap();
    let yard = 1 - room;
    assert_eq!(
        faces.locate(point(2.0, 1.5)).unwrap(),
        FaceLocation::Inside(yard)
    );
    assert_eq!(
        faces.locate(point(0.5, 0.5)).unwrap(),
        FaceLocation::Inside(room)
    );
    assert_eq!(
        faces.locate(point(1.1, 1.5)).unwrap(),
        FaceLocation::InObstacle(1)
    );
}

#[test]
fn two_flush_walls_on_one_line_keep_the_vertex_where_the_wall_changes() {
    let faces = faces_of(
        vec![
            rect(-0.25, -0.25, 2.0, 0.0),
            rect(2.0, -0.25, 4.25, 0.0),
            rect(-0.25, 3.0, 4.25, 3.25),
            rect(-0.25, 0.0, 0.0, 3.0),
            rect(4.0, 0.0, 4.25, 3.0),
        ],
        Vec::new(),
        GAP,
    );
    let room = &faces.faces()[0];
    assert_eq!(
        room.outer.points,
        vec![
            point(0.0, 0.0),
            point(2.0, 0.0),
            point(4.0, 0.0),
            point(4.0, 3.0),
            point(0.0, 3.0)
        ]
    );
    let owner = |edge: usize| match room.outer.sources[edge][..] {
        [PieceSource::Obstacle { polygon, .. }] => polygon,
        _ => panic!("edge {edge} lies on one obstacle edge"),
    };
    assert_eq!([owner(0), owner(1)], [0, 1]);
}

#[test]
fn a_point_on_the_wall_edge_of_a_closed_gap_is_in_the_gap() {
    let faces = faces_of(mitred_walls(D), Vec::new(), GAP);
    // The two mitre faces that bound a gap, the south wall's east end and
    // the east wall's south end: one of them is on the side of the gap that
    // an inside test by winding leaves out.
    for (a, b) in [
        (point(4.0 - D, 0.0), point(4.25 - D, -0.25)),
        (point(4.0, D), point(4.25, -0.25 + D)),
    ] {
        let middle = point((a.x + b.x) / 2.0, (a.y + b.y) / 2.0);
        assert!(matches!(
            faces.locate(middle).unwrap(),
            FaceLocation::InClosedGap(_)
        ));
    }
}

#[test]
fn a_ring_that_visits_its_smallest_point_twice_starts_with_the_lowest_rotation() {
    // The island shares the room's corner (0, 0), the smallest point of the
    // ring, which the ring visits once for each of the two wedges of free
    // space there.
    let island = polygon(&[(0.0, 0.0), (1.0, 0.5), (0.5, 1.0)]);
    let faces = faces_of(
        vec![
            walled(&[(0.0, 0.0), (4.0, 0.0), (4.0, 3.0), (0.0, 3.0)]),
            island,
        ],
        Vec::new(),
        GAP,
    );
    assert_eq!(faces.faces().len(), 1);
    let points = &faces.faces()[0].outer.points;
    let starts: Vec<usize> = (0..points.len())
        .filter(|&i| points[i] == point(0.0, 0.0))
        .collect();
    assert_eq!(starts.len(), 2);
    assert_eq!(starts[0], 0);
    let rotated = |from: usize| -> Vec<(f64, f64)> {
        (0..points.len())
            .map(|i| points[(from + i) % points.len()])
            .map(|p| (p.x, p.y))
            .collect()
    };
    let (own, other) = (rotated(0), rotated(starts[1]));
    assert!(own.partial_cmp(&other) == Some(std::cmp::Ordering::Less));
}

/// Two walls whose facing edges are `separation` apart, near the origin.
fn two_walls(separation: f64) -> Vec<Polygon> {
    vec![rect(-2.0, 0.0, 0.0, 1.0), rect(separation, 0.0, 2.0, 1.0)]
}

#[test]
fn a_distance_exactly_the_gap_is_not_bridged_and_one_ulp_below_is() {
    let gap = 1.0 / 1024.0;
    let at = faces_of(two_walls(gap), Vec::new(), gap);
    assert!(at.bridges().is_empty());
    let below = faces_of(two_walls(gap.next_down()), Vec::new(), gap);
    assert_eq!(below.bridges().len(), 2);
    // The strip between them is open at both ends: two bridges close it.
    assert_eq!(below.closed_gaps().len(), 1);
}

#[test]
fn gap_zero_closes_nothing() {
    let mitred = faces_of(mitred_walls(D), Vec::new(), 0.0);
    assert!(mitred.faces().is_empty() && mitred.bridges().is_empty());
    assert!(mitred.closed_gaps().is_empty());
    let room = faces_of(closed_room(), Vec::new(), 0.0);
    assert_eq!(room.faces().len(), 1);
    assert!(room.bridges().is_empty() && room.closed_gaps().is_empty());
}

#[test]
fn a_barrier_between_a_vertex_and_a_wall_takes_the_bridge() {
    let pillar = polygon(&[(0.0, 6e-5), (0.5, 0.5), (0.0, 1.0), (-0.5, 0.5)]);
    let faces = faces_of(
        vec![
            walled(&[(-4.0, 0.0), (4.0, 0.0), (4.0, 3.0), (-4.0, 3.0)]),
            pillar,
        ],
        vec![barrier((-1.0, 3e-5), (1.0, 3e-5))],
        GAP,
    );
    let to_pillar: Vec<_> = faces
        .bridges()
        .iter()
        .filter(|bridge| bridge.ends.contains(&point(0.0, 6e-5)))
        .collect();
    assert_eq!(to_pillar.len(), 1);
    assert!(to_pillar[0].ends.contains(&point(0.0, 3e-5)));
    assert!(to_pillar[0]
        .joins
        .iter()
        .any(|join| join.contains(&PieceSource::Barrier { barrier: 0 })));
    // The barrier's own ends reach the wall.
    assert_eq!(faces.bridges().len(), 3);
}

#[test]
fn a_wall_thinner_than_the_gap_lists_no_bridge_inside_it() {
    let faces = faces_of(
        vec![rect(0.0, 0.0, 4.0, 5e-5), rect(2.0, 5e-5, 2.5, 3.0)],
        Vec::new(),
        GAP,
    );
    assert!(faces.bridges().is_empty());
    assert!(faces.faces().is_empty());
}

#[test]
fn a_corridor_narrower_than_the_gap_between_two_rooms_is_a_closed_gap() {
    let first = Polygon {
        outer: ring(&[(-1.0, -1.0), (5.0, -1.0), (5.0, 4.0), (-1.0, 4.0)]),
        holes: vec![ring(&[(0.0, 0.0), (4.0, 0.0), (4.0, 3.0), (0.0, 3.0)])],
    };
    let second = Polygon {
        outer: ring(&[
            (5.0 + 5e-5, -1.0),
            (11.0, -1.0),
            (11.0, 4.0),
            (5.0 + 5e-5, 4.0),
        ]),
        holes: vec![ring(&[(6.0, 0.0), (10.0, 0.0), (10.0, 3.0), (6.0, 3.0)])],
    };
    let faces = faces_of(vec![first, second], Vec::new(), GAP);
    assert_eq!(faces.faces().len(), 2);
    assert_eq!(faces.closed_gaps().len(), 1);
    assert_eq!(faces.bridges().len(), 2);
}

/// An equilateral hole of inradius 1 round the origin and a concentric
/// island of inradius `1 - w`: a wall `w` wide all round, at distance `2 w`
/// corner to corner.
fn triangle_ring(inradius: f64) -> Vec<(f64, f64)> {
    let r = 2.0 * inradius;
    let c = 3f64.sqrt() / 2.0 * r;
    vec![(0.0, r), (-c, -r / 2.0), (c, -r / 2.0)]
}

fn wrapped_triangle(w: f64) -> FreeFaces {
    faces_of(
        vec![
            walled(&triangle_ring(1.0)),
            polygon(&triangle_ring(1.0 - w)),
        ],
        Vec::new(),
        GAP,
    )
}

#[test]
fn a_wall_narrower_than_the_gap_all_round_is_gaps_and_wider_corner_cells_are_faces() {
    // At a corner of angle 60 degrees wrapped at clearance w, bridging the
    // island's corner to the two outer sides cuts off a kite of width
    // 2 w cos^2(theta / 2): a gap while that is under the closing width, a
    // face once it exceeds it (ADR 0109).
    let narrow = wrapped_triangle(0.6e-4);
    assert_eq!(narrow.faces().len(), 0);
    assert_eq!(narrow.closed_gaps().len(), 6);
    let wide = wrapped_triangle(0.8e-4);
    assert_eq!(wide.faces().len(), 3);
    assert_eq!(wide.closed_gaps().len(), 3);
}

#[test]
fn a_gap_that_rounds_away_is_counted() {
    // A sliver 1e-8 high: two bridges close it, and its ring is gone after
    // snap rounding.
    let faces = faces_of(
        vec![walled(&[
            (1.0, 1.0),
            (1.5, 1.0),
            (5.0, 1.0),
            (5.0, 1.0 + 1e-8),
        ])],
        Vec::new(),
        GAP,
    );
    assert!(faces.faces().is_empty());
    assert_eq!(faces.closed_gaps().len(), 2);
    assert_eq!(faces.bridges().len(), 2);
    assert_eq!(faces.evidence().collapsed_gaps, 1);
}

#[test]
fn a_bridge_from_a_projection_on_a_vertex_joins_everything_at_that_vertex() {
    // (0, 1) projects on the end (0, 0) of the second barrier exactly; two
    // barriers end there.
    let faces = faces_of(
        Vec::new(),
        vec![
            barrier((0.0, 1.0), (-3.0, -0.5)),
            barrier((0.0, 0.0), (5.0, 0.0)),
            barrier((0.0, 0.0), (1.0, 3.0)),
        ],
        2.0,
    );
    let bridge = faces
        .bridges()
        .iter()
        .find(|bridge| bridge.ends == [point(0.0, 0.0), point(0.0, 1.0)])
        .expect("the bridge from (0, 1) to (0, 0)");
    assert_eq!(
        bridge.joins[0],
        vec![
            PieceSource::Barrier { barrier: 1 },
            PieceSource::Barrier { barrier: 2 }
        ]
    );
    assert_eq!(bridge.joins[1], vec![PieceSource::Barrier { barrier: 0 }]);
}

#[test]
fn a_sharp_corner_is_not_cut_but_a_collinear_vertex_within_the_gap_of_the_other_edge_cuts_it() {
    let sharp = |ring: &[(f64, f64)]| faces_of(vec![walled(ring)], Vec::new(), GAP);
    let plain = sharp(&[(0.0, 0.0), (4.0, 0.0), (4.0, 3e-4)]);
    assert_eq!(plain.faces().len(), 1);
    assert_eq!(plain.faces()[0].outer.points.len(), 3);
    assert!(plain.bridges().is_empty() && plain.closed_gaps().is_empty());
    let marked = sharp(&[(0.0, 0.0), (0.5, 0.0), (4.0, 0.0), (4.0, 3e-4)]);
    assert_eq!(marked.faces().len(), 1);
    assert_eq!(marked.faces()[0].outer.points.len(), 4);
    assert_eq!(marked.closed_gaps().len(), 1);
    assert_eq!(marked.bridges().len(), 1);
}

#[test]
fn overlapping_obstacles_and_barriers_along_edges_list_every_source() {
    // Walls that overlap a lot.
    let mut walls = closed_room();
    walls.push(rect(-0.25, -0.25, 1.0, 1.0));
    let faces = faces_of(walls, Vec::new(), GAP);
    assert_eq!(faces.faces().len(), 1);
    // The corner block cuts the room's corner: six corners.
    assert_eq!(faces.faces()[0].outer.points.len(), 6);
    // A barrier along the south wall's inner edge, twice.
    let barriers = vec![
        barrier((0.0, 0.0), (4.0, 0.0)),
        barrier((0.0, 0.0), (4.0, 0.0)),
    ];
    let faces = faces_of(closed_room(), barriers, GAP);
    let room = &faces.faces()[0];
    let south = room
        .outer
        .points
        .iter()
        .position(|&p| p == point(0.0, 0.0))
        .unwrap();
    let sources = &room.outer.sources[south];
    assert_eq!(sources.len(), 3);
    assert!(matches!(sources[0], PieceSource::Obstacle { .. }));
    assert_eq!(sources[1], PieceSource::Barrier { barrier: 0 });
    assert_eq!(sources[2], PieceSource::Barrier { barrier: 1 });
}

#[test]
fn a_barrier_through_walls_splits_the_room_and_only_its_free_part_is_an_edge() {
    let faces = faces_of(closed_room(), vec![barrier((2.0, -5.0), (2.0, 8.0))], GAP);
    assert_eq!(faces.faces().len(), 2);
    for face in faces.faces() {
        let barrier_edges = face
            .outer
            .sources
            .iter()
            .filter(|sources| matches!(sources[..], [PieceSource::Barrier { barrier: 0 }]))
            .count();
        assert_eq!(barrier_edges, 1);
        assert_eq!(face.outer.points.len(), 4);
    }
    // A barrier that crosses a wall and ends inside the room bounds nothing.
    let faces = faces_of(closed_room(), vec![barrier((2.0, -5.0), (2.0, 1.0))], GAP);
    assert_eq!(faces.faces().len(), 1);
    assert_eq!(faces.faces()[0].outer.points.len(), 4);
}

#[test]
fn the_mitre_fixture_far_from_the_origin_gives_the_same_result_translated() {
    let shift = 1_048_576.0;
    let moved = |polygon: &Polygon| Polygon {
        outer: Ring {
            points: polygon
                .outer
                .points
                .iter()
                .map(|p| point(p.x + shift, p.y + shift))
                .collect(),
        },
        holes: Vec::new(),
    };
    let near = faces_of(mitred_walls(D), Vec::new(), GAP);
    let far = faces_of(mitred_walls(D).iter().map(moved).collect(), Vec::new(), GAP);
    let shifted = |points: &[Point2]| -> Vec<Point2> {
        points
            .iter()
            .map(|p| point(p.x + shift, p.y + shift))
            .collect()
    };
    assert_eq!(far.faces().len(), 1);
    assert_eq!(
        far.faces()[0].outer.points,
        shifted(&near.faces()[0].outer.points)
    );
    assert_eq!(far.faces()[0].outer.sources, near.faces()[0].outer.sources);
    assert_eq!(far.bridges().len(), near.bridges().len());
    for (a, b) in far.bridges().iter().zip(near.bridges()) {
        assert_eq!(a.ends.to_vec(), shifted(&b.ends));
    }
    assert_eq!(far.closed_gaps().len(), 4);
}

// --- permuted input --------------------------------------------------------

/// What a result says, with polygon and barrier indices mapped back to the
/// original order: `polygons[k]` is the original index of input `k`.
type Summary = Vec<(Vec<Point2>, Vec<Vec<PieceSource>>)>;

fn summary(faces: &FreeFaces, polygons: &[usize], barriers: &[usize]) -> Summary {
    let map = |list: &[PieceSource]| -> Vec<PieceSource> {
        let mut mapped: Vec<PieceSource> = list
            .iter()
            .map(|&source| match source {
                PieceSource::Obstacle {
                    polygon,
                    ring,
                    edge,
                } => PieceSource::Obstacle {
                    polygon: polygons[polygon],
                    ring,
                    edge,
                },
                PieceSource::Barrier { barrier } => PieceSource::Barrier {
                    barrier: barriers[barrier],
                },
                other => other,
            })
            .collect();
        mapped.sort_unstable();
        mapped
    };
    let mut out = Vec::new();
    let mut ring = |ring: &FreeRing| {
        out.push((
            ring.points.clone(),
            ring.sources.iter().map(|s| map(s)).collect(),
        ));
    };
    for face in faces.faces().iter().chain(faces.closed_gaps()) {
        ring(&face.outer);
        face.holes.iter().for_each(&mut ring);
    }
    for bridge in faces.bridges() {
        out.push((
            bridge.ends.to_vec(),
            bridge.joins.iter().map(|join| map(join)).collect(),
        ));
    }
    out
}

fn identity(count: usize) -> Vec<usize> {
    (0..count).collect()
}

/// A seeded xorshift generator.
struct Rng(u64);

impl Rng {
    fn next(&mut self) -> u64 {
        self.0 ^= self.0 << 13;
        self.0 ^= self.0 >> 7;
        self.0 ^= self.0 << 17;
        self.0
    }

    fn unit(&mut self) -> f64 {
        (self.next() >> 11) as f64 / (1u64 << 53) as f64
    }

    fn between(&mut self, low: f64, high: f64) -> f64 {
        low + (high - low) * self.unit()
    }

    fn below(&mut self, bound: usize) -> usize {
        (self.next() % bound as u64) as usize
    }

    fn shuffle(&mut self, count: usize) -> Vec<usize> {
        let mut order = identity(count);
        for i in (1..count).rev() {
            order.swap(i, self.below(i + 1));
        }
        order
    }
}

/// The input in another order, and the original index of each entry.
fn permuted(source: &FreeSpaceInput, rng: &mut Rng) -> (FreeSpaceInput, Vec<usize>, Vec<usize>) {
    let polygons = rng.shuffle(source.obstacles.len());
    let barriers = rng.shuffle(source.barriers.len());
    let shuffled = FreeSpaceInput {
        obstacles: polygons
            .iter()
            .map(|&k| source.obstacles[k].clone())
            .collect(),
        barriers: barriers.iter().map(|&k| source.barriers[k]).collect(),
    };
    (shuffled, polygons, barriers)
}

fn assert_order_free(source: &FreeSpaceInput, gap: f64, seed: u64) {
    let direct = free_faces(source, gap, t()).unwrap();
    let mut rng = Rng(seed | 1);
    let (shuffled, polygons, barriers) = permuted(source, &mut rng);
    let other = free_faces(&shuffled, gap, t()).unwrap();
    assert_eq!(
        summary(
            &direct,
            &identity(source.obstacles.len()),
            &identity(source.barriers.len())
        ),
        summary(&other, &polygons, &barriers),
        "seed {seed}"
    );
}

#[test]
fn a_permuted_input_gives_the_same_rings_and_the_same_sources() {
    let mut obstacles = mitred_walls(D);
    obstacles.extend(mitred_room(6.0, 0.0, (4.0, 3.0), 0.25, [D; 8]));
    obstacles.push(rect(1.0, 1.0, 2.0, 2.0));
    let barriers = vec![
        barrier((2.0, 0.0), (2.0, 1.0)),
        barrier((3.0, 1.5), (3.5, 1.5)),
    ];
    let source = input(obstacles, barriers);
    for seed in 1..=20 {
        assert_order_free(&source, GAP, seed);
    }
}

// --- noisy floor plans -----------------------------------------------------

fn through_f32(value: f64) -> f64 {
    f64::from(value as f32)
}

fn cast_polygon(polygon: &Polygon) -> Polygon {
    let cast = |ring: &Ring| Ring {
        points: ring
            .points
            .iter()
            .map(|p| point(through_f32(p.x), through_f32(p.y)))
            .collect(),
    };
    Polygon {
        outer: cast(&polygon.outer),
        holes: polygon.holes.iter().map(cast).collect(),
    }
}

/// A floor plan: what blocks, and the rectangles of its rooms.
struct Plan {
    source: FreeSpaceInput,
    rooms: Vec<[f64; 4]>,
}

const HALF: f64 = 0.125;

fn column(i: usize) -> f64 {
    4.5 * i as f64
}

fn row(j: usize) -> f64 {
    3.5 * j as f64
}

/// Horizontal wall strips full width and vertical wall pieces between them,
/// each end of a vertical piece short of its strip by `shorten` (negative:
/// into it), in a grid of `nx` by `ny` rooms.
fn strip_plan(nx: usize, ny: usize, mut shorten: impl FnMut() -> f64) -> Plan {
    let mut obstacles = Vec::new();
    for j in 0..=ny {
        obstacles.push(rect(
            column(0) - HALF,
            row(j) - HALF,
            column(nx) + HALF,
            row(j) + HALF,
        ));
    }
    for i in 0..=nx {
        for j in 0..ny {
            let (low, high) = (shorten(), shorten());
            obstacles.push(rect(
                column(i) - HALF,
                row(j) + HALF + low,
                column(i) + HALF,
                row(j + 1) - HALF - high,
            ));
        }
    }
    let mut rooms = Vec::new();
    for i in 0..nx {
        for j in 0..ny {
            rooms.push([
                column(i) + HALF,
                row(j) + HALF,
                column(i + 1) - HALF,
                row(j + 1) - HALF,
            ]);
        }
    }
    Plan {
        source: FreeSpaceInput {
            obstacles: obstacles.iter().map(cast_polygon).collect(),
            barriers: Vec::new(),
        },
        rooms,
    }
}

/// A room of its own for each cell, walls mitred, every end pulled back by
/// its own noise.
fn mitre_plan(rng: &mut Rng) -> Plan {
    let mut obstacles = Vec::new();
    let mut rooms = Vec::new();
    for i in 0..3 {
        for j in 0..3 {
            let (x0, y0) = (5.0 * i as f64, 4.0 * j as f64);
            let mut pulls = [0.0; 8];
            pulls.iter_mut().for_each(|p| *p = rng.between(-1e-5, 1e-5));
            obstacles.extend(mitred_room(x0, y0, (4.0, 3.0), 0.25, pulls));
            rooms.push([x0, y0, x0 + 4.0, y0 + 3.0]);
        }
    }
    Plan {
        source: FreeSpaceInput {
            obstacles: obstacles.iter().map(cast_polygon).collect(),
            barriers: Vec::new(),
        },
        rooms,
    }
}

/// Wall strips round the floor and barriers for the partitions, each barrier
/// end short of the wall or past its near face by up to 1e-5.
fn barrier_plan(rng: &mut Rng) -> Plan {
    let (width, height) = (column(3), row(3));
    let mut obstacles = vec![
        rect(-HALF, -HALF, width + HALF, HALF),
        rect(-HALF, height - HALF, width + HALF, height + HALF),
        rect(-HALF, -HALF, HALF, height + HALF),
        rect(width - HALF, -HALF, width + HALF, height + HALF),
    ];
    obstacles = obstacles.iter().map(cast_polygon).collect();
    let end = |x: f64| through_f32(x);
    let mut barriers = Vec::new();
    for i in 1..3 {
        let x = column(i);
        barriers.push([
            point(end(x), end(HALF + rng.between(-1e-5, 1e-5))),
            point(end(x), end(height - HALF - rng.between(-1e-5, 1e-5))),
        ]);
    }
    for j in 1..3 {
        let y = row(j);
        barriers.push([
            point(end(HALF + rng.between(-1e-5, 1e-5)), end(y)),
            point(end(width - HALF - rng.between(-1e-5, 1e-5)), end(y)),
        ]);
    }
    let edges_x = [HALF, column(1), column(2), width - HALF];
    let edges_y = [HALF, row(1), row(2), height - HALF];
    let mut rooms = Vec::new();
    for i in 0..3 {
        for j in 0..3 {
            rooms.push([edges_x[i], edges_y[j], edges_x[i + 1], edges_y[j + 1]]);
        }
    }
    Plan {
        source: FreeSpaceInput {
            obstacles,
            barriers,
        },
        rooms,
    }
}

/// What a floor plan must give, for any noise: one face per room of about
/// its area, each bridge shorter than the gap, every room's centre inside
/// its own face, and the same rings for a shuffled input.
fn check_plan(plan: &Plan, seed: u64) -> usize {
    let faces = free_faces(&plan.source, GAP, t()).unwrap();
    assert_eq!(faces.faces().len(), plan.rooms.len(), "seed {seed}");
    assert!(
        faces
            .bridges()
            .iter()
            .all(|bridge| bridge.length <= GAP * (1.0 + 1e-12)),
        "seed {seed}"
    );
    let mut seen = Vec::new();
    for &[x0, y0, x1, y1] in &plan.rooms {
        let centre = point((x0 + x1) / 2.0, (y0 + y1) / 2.0);
        let FaceLocation::Inside(index) = faces.locate(centre).unwrap() else {
            panic!("seed {seed}: a room centre is not inside a face");
        };
        let found = area(&faces.faces()[index].outer.points);
        let nominal = (x1 - x0) * (y1 - y0);
        let perimeter = 2.0 * ((x1 - x0) + (y1 - y0));
        assert!(
            (found - nominal).abs() <= perimeter * GAP,
            "seed {seed}: area {found} for {nominal}"
        );
        assert!(faces.faces()[index].holes.is_empty(), "seed {seed}");
        seen.push(index);
    }
    seen.sort_unstable();
    seen.dedup();
    assert_eq!(seen.len(), plan.rooms.len(), "seed {seed}");
    assert_order_free(&plan.source, GAP, seed);
    faces.bridges().len()
}

#[test]
fn noisy_mitre_joints_give_one_face_per_room() {
    let mut bridges = 0;
    for seed in 1..=200u64 {
        let mut rng = Rng(seed.wrapping_mul(0x9E37_79B9_7F4A_7C15) | 1);
        bridges += check_plan(&mitre_plan(&mut rng), seed);
    }
    assert!(bridges > 0, "the corners are open");
}

#[test]
fn noisy_butt_joints_with_a_gap_give_one_face_per_room() {
    let mut bridges = 0;
    for seed in 1..=200u64 {
        let mut rng = Rng(seed.wrapping_mul(0x9E37_79B9_7F4A_7C15) | 1);
        bridges += check_plan(&strip_plan(3, 3, || rng.between(1e-6, 1e-5)), seed);
    }
    assert!(bridges >= 200 * 24, "every joint is open: {bridges}");
}

#[test]
fn noisy_butt_joints_that_overlap_give_one_face_per_room() {
    for seed in 1..=200u64 {
        let mut rng = Rng(seed.wrapping_mul(0x9E37_79B9_7F4A_7C15) | 1);
        check_plan(&strip_plan(3, 3, || -rng.between(1e-6, 1e-5)), seed);
    }
}

#[test]
fn noisy_t_junctions_give_one_face_per_room() {
    for seed in 1..=200u64 {
        let mut rng = Rng(seed.wrapping_mul(0x9E37_79B9_7F4A_7C15) | 1);
        check_plan(&strip_plan(3, 3, || rng.between(-1e-6, 1e-6)), seed);
    }
}

#[test]
fn barrier_ends_short_of_or_past_a_wall_give_one_face_per_room() {
    let mut bridges = 0;
    for seed in 1..=200u64 {
        let mut rng = Rng(seed.wrapping_mul(0x9E37_79B9_7F4A_7C15) | 1);
        bridges += check_plan(&barrier_plan(&mut rng), seed);
    }
    assert!(bridges > 0, "some barrier ends are short");
}

#[test]
fn the_candidate_pairs_of_a_floor_grow_with_its_size_and_not_its_square() {
    let pairs = |n: usize| {
        let plan = strip_plan(n, n, || D);
        let faces = free_faces(&plan.source, GAP, t()).unwrap();
        assert_eq!(faces.faces().len(), n * n);
        faces.evidence().candidate_pairs
    };
    let (small, large) = (pairs(10), pairs(20));
    assert!(
        large < 5 * small,
        "{small} pairs for 10x10, {large} for 20x20"
    );
}
