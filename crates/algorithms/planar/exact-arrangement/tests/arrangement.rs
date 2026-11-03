//! Arrangement and ring tracing through the public API only.

use axiolid_contracts::{CancellationToken, ExecutionOptions, GeomError};
use axiolid_core::{Point2, Tolerance};
use axiolid_exact_arrangement::{
    arrange, trace, Arrangement, ExactPoints, HalfEdge, PairBudget, PieceId, PointId, Segment,
    VertexId,
};

type Line = ((f64, f64), (f64, f64));

fn options() -> ExecutionOptions {
    ExecutionOptions::new(Tolerance::METRE)
}

fn budget(max: usize) -> PairBudget {
    PairBudget::new(max, "test pairs")
}

/// Segments between positions; each position becomes a point of its own.
fn build(
    corners: &[Line],
    options: &ExecutionOptions,
    max: usize,
) -> (
    ExactPoints,
    Vec<Segment>,
    Result<(Arrangement, usize), GeomError>,
) {
    let mut points = ExactPoints::new();
    let segments: Vec<Segment> = corners
        .iter()
        .map(|&((ax, ay), (bx, by))| Segment {
            ends: [
                points.push_flat(Point2::new(ax, ay)).unwrap(),
                points.push_flat(Point2::new(bx, by)).unwrap(),
            ],
        })
        .collect();
    let mut pairs = budget(max);
    let result = arrange(&mut points, &segments, &mut pairs, 0, options).map(|a| (a, pairs.used));
    (points, segments, result)
}

fn arranged(corners: &[Line]) -> (ExactPoints, Arrangement) {
    let (points, _, result) = build(corners, &options(), 1000);
    (points, result.expect("arrangement").0)
}

#[test]
fn two_crossing_segments_make_four_pieces_around_a_crossing() {
    let (points, arrangement) = arranged(&[((0.0, 0.0), (2.0, 2.0)), ((2.0, 0.0), (0.0, 2.0))]);
    assert_eq!(arrangement.pieces.len(), 4);
    assert_eq!(arrangement.crossings, 1);
    let crossing = arrangement
        .vertices
        .iter()
        .find(|vertex| vertex.inputs.is_empty())
        .expect("crossing vertex");
    assert_eq!(points.round(crossing.point).unwrap(), Point2::new(1.0, 1.0));
    assert_eq!(crossing.crossing.unwrap().segments, [0, 1]);
    // Every segment runs over two pieces, in order.
    assert!(arrangement.paths.iter().all(|path| path.len() == 2));
}

#[test]
fn a_t_junction_splits_the_segment_it_ends_on() {
    let (_, arrangement) = arranged(&[((0.0, 0.0), (4.0, 0.0)), ((2.0, 0.0), (2.0, 3.0))]);
    assert_eq!(arrangement.pieces.len(), 3);
    assert_eq!(arrangement.crossings, 0);
    assert_eq!(arrangement.paths[0].len(), 2);
    assert_eq!(arrangement.paths[1].len(), 1);
}

#[test]
fn a_collinear_partial_overlap_gives_three_pieces_and_the_shared_one_has_both_segments() {
    let (_, arrangement) = arranged(&[((0.0, 0.0), (3.0, 0.0)), ((1.0, 0.0), (5.0, 0.0))]);
    assert_eq!(arrangement.pieces.len(), 3);
    let shared: Vec<_> = arrangement
        .pieces
        .iter()
        .filter(|piece| piece.segments == [0, 1])
        .collect();
    assert_eq!(shared.len(), 1);
}

#[test]
fn equal_points_with_different_ids_share_a_vertex() {
    // Three segments through (1, 1), one of them ending there.
    let (_, arrangement) = arranged(&[
        ((0.0, 0.0), (2.0, 2.0)),
        ((2.0, 0.0), (0.0, 2.0)),
        ((1.0, 1.0), (1.0, 5.0)),
        ((1.0, 1.0), (3.0, 1.0)),
    ]);
    let centre: Vec<_> = arrangement
        .vertices
        .iter()
        .filter(|vertex| vertex.inputs.len() == 2)
        .collect();
    assert_eq!(centre.len(), 1);
    assert_eq!(arrangement.crossings, 0);
    // Four spokes leave the shared vertex... and the third and fourth ends.
    assert_eq!(
        arrangement.paths.iter().map(Vec::len).collect::<Vec<_>>(),
        [2, 2, 1, 1]
    );
}

#[test]
fn vertices_are_numbered_by_input_point_then_by_crossing() {
    let (_, arrangement) = arranged(&[((0.0, 0.0), (2.0, 2.0)), ((2.0, 0.0), (0.0, 2.0))]);
    let inputs: Vec<bool> = arrangement
        .vertices
        .iter()
        .map(|vertex| !vertex.inputs.is_empty())
        .collect();
    assert_eq!(inputs, [true, true, true, true, false]);
    let firsts: Vec<PointId> = arrangement.vertices[..4]
        .iter()
        .map(|v| v.inputs[0])
        .collect();
    assert!(firsts.windows(2).all(|pair| pair[0] < pair[1]));
    assert!(arrangement
        .pieces
        .windows(2)
        .all(|pair| pair[0].ends < pair[1].ends));
}

#[test]
fn the_fan_of_a_vertex_runs_counter_clockwise() {
    // A plus sign; the fan at its centre starts at +x... half-plane order.
    let (points, arrangement) = arranged(&[((-1.0, 0.0), (1.0, 0.0)), ((0.0, -1.0), (0.0, 1.0))]);
    let centre = arrangement
        .vertices
        .iter()
        .position(|vertex| vertex.inputs.is_empty())
        .unwrap();
    let fan = arrangement.fan(&points, VertexId(centre)).unwrap();
    let directions: Vec<Point2> = fan
        .iter()
        .map(|&(piece, forward)| {
            let to = arrangement.pieces[piece.0].ends[usize::from(forward)];
            points.round(arrangement.vertices[to.0].point).unwrap()
        })
        .collect();
    assert_eq!(
        directions,
        [
            Point2::new(1.0, 0.0),
            Point2::new(0.0, 1.0),
            Point2::new(-1.0, 0.0),
            Point2::new(0.0, -1.0)
        ]
    );
}

#[test]
fn every_examined_pair_counts_and_the_budget_refuses_beyond_it() {
    let square = [
        ((0.0, 0.0), (1.0, 0.0)),
        ((1.0, 0.0), (1.0, 1.0)),
        ((1.0, 1.0), (0.0, 1.0)),
        ((0.0, 1.0), (0.0, 0.0)),
    ];
    let (_, _, result) = build(&square, &options(), 1000);
    let used = result.expect("arrangement").1;
    assert!(used > 0);
    let (_, _, cramped) = build(&square, &options(), used - 1);
    assert!(matches!(
        cramped,
        Err(GeomError::BudgetExceeded {
            resource: "test pairs"
        })
    ));
    assert!(build(&square, &options(), used).2.is_ok());
}

#[test]
fn the_crossing_charge_refuses_beyond_the_memory_budget() {
    let cross = [((0.0, 0.0), (2.0, 2.0)), ((2.0, 0.0), (0.0, 2.0))];
    let fits = options().with_memory_budget(1 << 30);
    assert!(build(&cross, &fits, 1000).2.is_ok());
    let tight = options().with_memory_budget(100);
    assert!(matches!(
        build(&cross, &tight, 1000).2,
        Err(GeomError::BudgetExceeded { resource: "memory" })
    ));
    // Segments that do not cross are not charged.
    let apart = [((0.0, 0.0), (1.0, 0.0)), ((0.0, 5.0), (1.0, 5.0))];
    assert!(build(&apart, &tight, 1000).2.is_ok());
}

#[test]
fn a_cancelled_arrangement_stops() {
    let token = CancellationToken::new();
    token.cancel();
    let cancelled = options().with_cancellation(token);
    let cross = [((0.0, 0.0), (2.0, 2.0)), ((2.0, 0.0), (0.0, 2.0))];
    assert!(matches!(
        build(&cross, &cancelled, 1000).2,
        Err(GeomError::Cancelled)
    ));
}

/// The half-edges of the pieces whose both ends lie on the given closed
/// polylines, directed along them.
fn ring_half_edges(
    points: &ExactPoints,
    arrangement: &Arrangement,
    rings: &[&[(f64, f64)]],
) -> Vec<HalfEdge> {
    let position = |vertex: VertexId| points.round(arrangement.vertices[vertex.0].point).unwrap();
    let mut out = Vec::new();
    for ring in rings {
        for index in 0..ring.len() {
            let (from, to) = (ring[index], ring[(index + 1) % ring.len()]);
            let found = arrangement
                .pieces
                .iter()
                .enumerate()
                .find_map(|(piece, p)| {
                    let (a, b) = (position(p.ends[0]), position(p.ends[1]));
                    let (from, to) = (Point2::new(from.0, from.1), Point2::new(to.0, to.1));
                    if (a, b) == (from, to) {
                        Some(HalfEdge {
                            piece: PieceId(piece),
                            forward: true,
                        })
                    } else if (a, b) == (to, from) {
                        Some(HalfEdge {
                            piece: PieceId(piece),
                            forward: false,
                        })
                    } else {
                        None
                    }
                });
            out.push(found.expect("piece of the ring"));
        }
    }
    out
}

fn closed(ring: &[(f64, f64)]) -> Vec<Line> {
    (0..ring.len())
        .map(|i| (ring[i], ring[(i + 1) % ring.len()]))
        .collect()
}

#[test]
fn a_square_is_one_outer_ring_and_collinear_vertices_are_dropped_unless_kept() {
    let ring = [(0.0, 0.0), (1.0, 0.0), (2.0, 0.0), (2.0, 2.0), (0.0, 2.0)];
    let (points, arrangement) = arranged(&closed(&ring));
    let half_edges = ring_half_edges(&points, &arrangement, &[&ring]);
    let (rings, regions) = trace(
        &points,
        &arrangement,
        &half_edges,
        &|_| false,
        &mut budget(100),
        &options(),
    )
    .unwrap();
    assert_eq!(rings.len(), 1);
    assert!(!rings[0].hole);
    assert_eq!(rings[0].half_edges.len(), 5);
    assert_eq!(rings[0].vertices.len(), 4);
    assert_eq!((regions[0].outer, regions[0].holes.len()), (0, 0));
    // Keeping every vertex keeps the one on the straight stretch.
    let (kept, _) = trace(
        &points,
        &arrangement,
        &half_edges,
        &|_| true,
        &mut budget(100),
        &options(),
    )
    .unwrap();
    assert_eq!(kept[0].vertices.len(), 5);
}

#[test]
fn a_hole_belongs_to_the_outer_ring_around_it_and_counts_pairs() {
    let outer = [(0.0, 0.0), (6.0, 0.0), (6.0, 6.0), (0.0, 6.0)];
    // Clockwise: the hole's own left is the material.
    let hole = [(2.0, 2.0), (2.0, 4.0), (4.0, 4.0), (4.0, 2.0)];
    let mut segments = closed(&outer);
    segments.extend(closed(&hole));
    let (points, arrangement) = arranged(&segments);
    let half_edges = ring_half_edges(&points, &arrangement, &[&outer, &hole]);
    let mut pairs = budget(100);
    let (rings, regions) = trace(
        &points,
        &arrangement,
        &half_edges,
        &|_| false,
        &mut pairs,
        &options(),
    )
    .unwrap();
    assert_eq!(rings.len(), 2);
    assert_eq!(rings.iter().filter(|ring| ring.hole).count(), 1);
    assert_eq!(regions.len(), 1);
    assert_eq!(regions[0].holes.len(), 1);
    assert!(rings[regions[0].holes[0]].hole);
    assert!(pairs.used > 0);
    let mut cramped = budget(0);
    assert!(matches!(
        trace(
            &points,
            &arrangement,
            &half_edges,
            &|_| false,
            &mut cramped,
            &options()
        ),
        Err(GeomError::BudgetExceeded { .. })
    ));
}

#[test]
fn two_squares_touching_at_a_corner_are_two_rings_not_one() {
    let first = [(0.0, 0.0), (1.0, 0.0), (1.0, 1.0), (0.0, 1.0)];
    let second = [(1.0, 1.0), (2.0, 1.0), (2.0, 2.0), (1.0, 2.0)];
    let mut segments = closed(&first);
    segments.extend(closed(&second));
    let (points, arrangement) = arranged(&segments);
    let half_edges = ring_half_edges(&points, &arrangement, &[&first, &second]);
    let (rings, regions) = trace(
        &points,
        &arrangement,
        &half_edges,
        &|_| true,
        &mut budget(100),
        &options(),
    )
    .unwrap();
    assert_eq!((rings.len(), regions.len()), (2, 2));
    assert!(rings.iter().all(|ring| ring.vertices.len() == 4));
}

#[test]
fn a_mirrored_frame_reads_a_counter_clockwise_ring_as_a_hole() {
    let ring = [(0.0, 0.0), (1.0, 0.0), (1.0, 1.0), (0.0, 1.0)];
    let mut points = ExactPoints::with_frame_sign(axiolid_contracts::Sign::Negative);
    let segments: Vec<Segment> = closed(&ring)
        .iter()
        .map(|&((ax, ay), (bx, by))| Segment {
            ends: [
                points.push_flat(Point2::new(ax, ay)).unwrap(),
                points.push_flat(Point2::new(bx, by)).unwrap(),
            ],
        })
        .collect();
    let arrangement = arrange(&mut points, &segments, &mut budget(100), 0, &options()).unwrap();
    let half_edges = ring_half_edges(&points, &arrangement, &[&ring]);
    // Counter-clockwise in (x, y) is clockwise in the mirrored frame: a hole
    // with no outer ring around it.
    assert!(matches!(
        trace(
            &points,
            &arrangement,
            &half_edges,
            &|_| false,
            &mut budget(100),
            &options()
        ),
        Err(GeomError::Degenerate(_))
    ));
}

/// The positions of the pieces' ends, each piece as a sorted pair of x.
fn pieces_along_x(points: &ExactPoints, arrangement: &Arrangement) -> Vec<(f64, f64)> {
    let mut found: Vec<(f64, f64)> = arrangement
        .pieces
        .iter()
        .map(|piece| {
            let [a, b] = piece
                .ends
                .map(|v| points.round(arrangement.vertices[v.0].point).unwrap().x);
            (a.min(b), a.max(b))
        })
        .collect();
    found.sort_by(|a, b| a.partial_cmp(b).unwrap());
    found
}

#[test]
fn segments_sharing_a_point_on_one_line_are_cut_where_the_shorter_ends() {
    for reversed in [false, true] {
        let mut points = ExactPoints::new();
        let a = points.push_flat(Point2::new(0.0, 0.0)).unwrap();
        let b = points.push_flat(Point2::new(1.0, 0.0)).unwrap();
        let c = points.push_flat(Point2::new(2.0, 0.0)).unwrap();
        let mut segments = vec![Segment { ends: [a, b] }, Segment { ends: [a, c] }];
        if reversed {
            segments.reverse();
        }
        let arrangement = arrange(&mut points, &segments, &mut budget(100), 0, &options()).unwrap();
        assert_eq!(
            pieces_along_x(&points, &arrangement),
            vec![(0.0, 1.0), (1.0, 2.0)],
            "reversed {reversed}"
        );
        // The vertex at B has one piece each way, and the fan lists both.
        let at_b = arrangement
            .vertices
            .iter()
            .position(|v| points.round(v.point).unwrap().x == 1.0)
            .unwrap();
        assert_eq!(arrangement.fan(&points, VertexId(at_b)).unwrap().len(), 2);
    }
}

#[test]
fn segments_sharing_their_far_end_are_cut_too() {
    let mut points = ExactPoints::new();
    let a = points.push_flat(Point2::new(0.0, 0.0)).unwrap();
    let b = points.push_flat(Point2::new(1.0, 0.0)).unwrap();
    let c = points.push_flat(Point2::new(2.0, 0.0)).unwrap();
    let segments = [Segment { ends: [b, c] }, Segment { ends: [a, c] }];
    let arrangement = arrange(&mut points, &segments, &mut budget(100), 0, &options()).unwrap();
    assert_eq!(
        pieces_along_x(&points, &arrangement),
        vec![(0.0, 1.0), (1.0, 2.0)]
    );
}

#[test]
fn the_same_segment_twice_is_one_piece() {
    let mut points = ExactPoints::new();
    let a = points.push_flat(Point2::new(0.0, 0.0)).unwrap();
    let b = points.push_flat(Point2::new(1.0, 0.0)).unwrap();
    let segments = [Segment { ends: [a, b] }, Segment { ends: [b, a] }];
    let arrangement = arrange(&mut points, &segments, &mut budget(100), 0, &options()).unwrap();
    assert_eq!(arrangement.pieces.len(), 1);
}

#[test]
fn splits_count_against_the_memory_budget() {
    // One long segment with many ends of other segments inside it.
    let mut points = ExactPoints::new();
    let a = points.push_flat(Point2::new(0.0, 0.0)).unwrap();
    let z = points.push_flat(Point2::new(100.0, 0.0)).unwrap();
    let mut segments = vec![Segment { ends: [a, z] }];
    for k in 1..50 {
        let on = points.push_flat(Point2::new(f64::from(k), 0.0)).unwrap();
        let off = points.push_flat(Point2::new(f64::from(k), 1.0)).unwrap();
        segments.push(Segment { ends: [on, off] });
    }
    let free = arrange(&mut points, &segments, &mut budget(100_000), 0, &options());
    assert!(free.is_ok());
    let tight = options().with_memory_budget(40 * axiolid_exact_arrangement::SPLIT_SCRATCH_BYTES);
    let refused = arrange(&mut points, &segments, &mut budget(100_000), 0, &tight);
    assert!(matches!(
        refused,
        Err(GeomError::BudgetExceeded { resource: "memory" })
    ));
}
