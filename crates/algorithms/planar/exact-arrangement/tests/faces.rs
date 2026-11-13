//! The faces of the complement of an arrangement, through the public API.

use axiolid_contracts::{CancellationToken, ExecutionOptions, GeomError, Sign};
use axiolid_core::{Point2, Tolerance};
use axiolid_exact_arrangement::{
    arrange, faces, Arrangement, ExactPoints, FaceId, Faces, PairBudget, PieceId, Segment,
};

type Pos = (f64, f64);
type Line = (Pos, Pos);

const OUT: FaceId = Faces::UNBOUNDED;

fn options() -> ExecutionOptions {
    ExecutionOptions::new(Tolerance::METRE)
}

struct Scene {
    points: ExactPoints,
    arrangement: Arrangement,
    faces: Faces,
}

fn scene_in(mut points: ExactPoints, lines: &[Line]) -> Scene {
    let segments: Vec<Segment> = lines
        .iter()
        .map(|&((ax, ay), (bx, by))| Segment {
            ends: [
                points.push_flat(Point2::new(ax, ay)).unwrap(),
                points.push_flat(Point2::new(bx, by)).unwrap(),
            ],
        })
        .collect();
    let mut pairs = PairBudget::new(1_000_000, "test pairs");
    let arrangement = arrange(&mut points, &segments, &mut pairs, 0, &options()).expect("arrange");
    let faces = faces(&points, &arrangement, &mut pairs, &options()).expect("faces");
    Scene {
        points,
        arrangement,
        faces,
    }
}

fn scene(lines: &[Line]) -> Scene {
    scene_in(ExactPoints::new(), lines)
}

/// The edges of a closed polyline.
fn ring(corners: &[Pos]) -> Vec<Line> {
    (0..corners.len())
        .map(|i| (corners[i], corners[(i + 1) % corners.len()]))
        .collect()
}

fn square(x: f64, y: f64, side: f64) -> Vec<Line> {
    ring(&[(x, y), (x + side, y), (x + side, y + side), (x, y + side)])
}

impl Scene {
    fn at(&self, vertex: usize) -> Pos {
        let point = self
            .points
            .round(self.arrangement.vertices[vertex].point)
            .unwrap();
        (point.x, point.y)
    }

    /// The faces left and right of the piece that runs from `from` to `to`.
    fn sides(&self, from: Pos, to: Pos) -> (FaceId, FaceId) {
        for (index, piece) in self.arrangement.pieces.iter().enumerate() {
            let ends = [self.at(piece.ends[0].0), self.at(piece.ends[1].0)];
            let id = PieceId(index);
            if ends == [from, to] {
                return (self.faces.left(id), self.faces.right(id));
            }
            if ends == [to, from] {
                return (self.faces.right(id), self.faces.left(id));
            }
        }
        panic!("no piece from {from:?} to {to:?}");
    }

    fn left(&self, from: Pos, to: Pos) -> FaceId {
        self.sides(from, to).0
    }

    fn right(&self, from: Pos, to: Pos) -> FaceId {
        self.sides(from, to).1
    }
}

#[test]
fn one_segment_has_the_unbounded_face_on_both_sides() {
    let s = scene(&[((0.0, 0.0), (2.0, 1.0))]);
    assert_eq!(s.faces.count(), 1);
    assert_eq!(s.sides((0.0, 0.0), (2.0, 1.0)), (OUT, OUT));
}

#[test]
fn a_triangle_has_its_inside_on_the_left_when_counter_clockwise() {
    let s = scene(&ring(&[(0.0, 0.0), (4.0, 0.0), (0.0, 4.0)]));
    assert_eq!(s.faces.count(), 2);
    assert_eq!(s.sides((0.0, 0.0), (4.0, 0.0)), (FaceId(1), OUT));
    assert_eq!(s.sides((4.0, 0.0), (0.0, 4.0)), (FaceId(1), OUT));
    assert_eq!(s.sides((0.0, 4.0), (0.0, 0.0)), (FaceId(1), OUT));
    assert_eq!(s.sides((4.0, 0.0), (0.0, 0.0)), (OUT, FaceId(1)));
}

#[test]
fn a_dangling_piece_has_one_face_on_both_sides() {
    let mut lines = ring(&[(0.0, 0.0), (8.0, 0.0), (0.0, 8.0)]);
    lines.push(((1.0, 1.0), (3.0, 1.0)));
    lines.push(((0.0, 0.0), (-2.0, -1.0)));
    // The inward spike hangs from the bottom edge's neighbour at (1, 0).
    lines.push(((1.0, 0.0), (1.0, 1.0)));
    let s = scene(&lines);
    assert_eq!(s.faces.count(), 2);
    assert_eq!(s.sides((1.0, 0.0), (1.0, 1.0)), (FaceId(1), FaceId(1)));
    assert_eq!(s.sides((1.0, 1.0), (3.0, 1.0)), (FaceId(1), FaceId(1)));
    assert_eq!(s.sides((0.0, 0.0), (-2.0, -1.0)), (OUT, OUT));
}

#[test]
fn disjoint_triangles_are_two_faces_and_the_unbounded_one() {
    let mut lines = ring(&[(0.0, 0.0), (2.0, 0.0), (0.0, 2.0)]);
    lines.extend(ring(&[(5.0, 0.0), (7.0, 0.0), (5.0, 2.0)]));
    let s = scene(&lines);
    assert_eq!(s.faces.count(), 3);
    assert_eq!(s.left((0.0, 0.0), (2.0, 0.0)), FaceId(1));
    assert_eq!(s.left((5.0, 0.0), (7.0, 0.0)), FaceId(2));
    assert_eq!(s.right((0.0, 0.0), (2.0, 0.0)), OUT);
    assert_eq!(s.right((5.0, 0.0), (7.0, 0.0)), OUT);
}

#[test]
fn a_triangle_inside_a_triangle_lies_in_the_face_of_the_outer_one() {
    let mut lines = ring(&[(0.0, 0.0), (20.0, 0.0), (0.0, 20.0)]);
    lines.extend(ring(&[(2.0, 2.0), (6.0, 2.0), (2.0, 6.0)]));
    let s = scene(&lines);
    assert_eq!(s.faces.count(), 3);
    let outer_inside = s.left((0.0, 0.0), (20.0, 0.0));
    let inner_inside = s.left((2.0, 2.0), (6.0, 2.0));
    assert_ne!(outer_inside, OUT);
    assert_ne!(inner_inside, outer_inside);
    assert_eq!(s.right((2.0, 2.0), (6.0, 2.0)), outer_inside);
    assert_eq!(s.right((0.0, 0.0), (20.0, 0.0)), OUT);
}

#[test]
fn the_innermost_of_three_nesting_levels_holds_a_component() {
    let mut lines = square(0.0, 0.0, 30.0);
    lines.extend(square(2.0, 2.0, 20.0));
    lines.extend(square(4.0, 4.0, 10.0));
    // Listed before the rings it lies in, so the order of the input does not
    // decide the nesting.
    let mut lines_reversed = lines.clone();
    lines_reversed.reverse();
    for lines in [lines, lines_reversed] {
        let s = scene(&lines);
        assert_eq!(s.faces.count(), 4);
        let first = s.left((0.0, 0.0), (30.0, 0.0));
        let second = s.left((2.0, 2.0), (22.0, 2.0));
        let third = s.left((4.0, 4.0), (14.0, 4.0));
        assert_eq!(s.right((0.0, 0.0), (30.0, 0.0)), OUT);
        assert_eq!(s.right((2.0, 2.0), (22.0, 2.0)), first);
        assert_eq!(s.right((4.0, 4.0), (14.0, 4.0)), second);
        let all = [OUT, first, second, third];
        for (i, face) in all.iter().enumerate() {
            assert!(all[i + 1..].iter().all(|other| other != face));
        }
    }
}

#[test]
fn an_isolated_segment_in_a_square_has_the_square_inside_on_both_sides() {
    let mut lines = square(0.0, 0.0, 10.0);
    lines.push(((3.0, 5.0), (7.0, 5.0)));
    let s = scene(&lines);
    assert_eq!(s.faces.count(), 2);
    assert_eq!(s.sides((3.0, 5.0), (7.0, 5.0)), (FaceId(1), FaceId(1)));
}

#[test]
fn a_component_inside_a_face_with_a_hole_picks_the_face_not_the_hole() {
    // A square with an inner square: the segment lies in the ring between
    // them, not in the inner square.
    let mut lines = square(0.0, 0.0, 10.0);
    lines.extend(square(4.0, 4.0, 2.0));
    lines.push(((1.0, 1.0), (2.0, 1.0)));
    let s = scene(&lines);
    assert_eq!(s.faces.count(), 3);
    let between = s.left((0.0, 0.0), (10.0, 0.0));
    assert_eq!(s.sides((1.0, 1.0), (2.0, 1.0)), (between, between));
    assert_ne!(s.left((4.0, 4.0), (6.0, 4.0)), between);
}

#[test]
fn two_squares_sharing_a_vertex_are_two_faces() {
    let mut lines = square(0.0, 0.0, 2.0);
    lines.extend(square(2.0, 2.0, 2.0));
    let s = scene(&lines);
    assert_eq!(s.faces.count(), 3);
    assert_ne!(
        s.left((0.0, 0.0), (2.0, 0.0)),
        s.left((2.0, 2.0), (4.0, 2.0))
    );
}

#[test]
fn two_squares_sharing_an_edge_have_different_faces_on_it() {
    let mut lines = square(0.0, 0.0, 2.0);
    lines.extend(square(2.0, 0.0, 2.0));
    let s = scene(&lines);
    assert_eq!(s.faces.count(), 3);
    let (left, right) = s.sides((2.0, 0.0), (2.0, 2.0));
    assert_ne!(left, right);
    assert_ne!(left, OUT);
    assert_ne!(right, OUT);
}

#[test]
fn a_four_by_four_grid_has_sixteen_cells() {
    let mut lines = Vec::new();
    for i in 0..5 {
        let i = f64::from(i);
        lines.push(((i, 0.0), (i, 4.0)));
        lines.push(((0.0, i), (4.0, i)));
    }
    let s = scene(&lines);
    assert_eq!(s.faces.count(), 17);
    let mut cells = Vec::new();
    for x in 0..4 {
        for y in 0..4 {
            let (x, y) = (f64::from(x), f64::from(y));
            cells.push(s.left((x, y), (x + 1.0, y)));
        }
    }
    cells.sort();
    cells.dedup();
    assert_eq!(cells.len(), 16);
    assert!(!cells.contains(&OUT));
}

#[test]
fn a_component_whose_extreme_vertex_has_one_piece_still_has_an_outer_contour() {
    let mut lines = vec![((0.0, 1.0), (1.0, 1.0))];
    lines.extend(ring(&[(1.0, 1.0), (3.0, 0.0), (3.0, 2.0)]));
    let s = scene(&lines);
    assert_eq!(s.faces.count(), 2);
    assert_eq!(s.sides((0.0, 1.0), (1.0, 1.0)), (OUT, OUT));
    assert_eq!(s.left((1.0, 1.0), (3.0, 0.0)), FaceId(1));
}

#[test]
fn a_mirrored_frame_reads_the_same_faces_with_the_sides_swapped() {
    let mut lines = square(0.0, 0.0, 30.0);
    lines.extend(square(2.0, 2.0, 10.0));
    lines.push(((14.0, 4.0), (20.0, 4.0)));
    let plain = scene(&lines);
    let mirrored = scene_in(ExactPoints::with_frame_sign(Sign::Negative), &lines);
    assert_eq!(mirrored.faces.count(), plain.faces.count());
    // Counter-clockwise in (x, y) is clockwise in the mirrored frame.
    assert_eq!(mirrored.left((0.0, 0.0), (30.0, 0.0)), OUT);
    let between = mirrored.right((0.0, 0.0), (30.0, 0.0));
    assert_ne!(between, OUT);
    assert_eq!(mirrored.left((2.0, 2.0), (12.0, 2.0)), between);
    assert_eq!(mirrored.sides((14.0, 4.0), (20.0, 4.0)), (between, between));
}

#[test]
fn a_segment_ending_at_a_computed_crossing_shares_its_vertex() {
    let mut points = ExactPoints::new();
    let mut push = |x: f64, y: f64| points.push_flat(Point2::new(x, y)).unwrap();
    let ends = [
        push(0.0, 0.0),
        push(4.0, 4.0),
        push(0.0, 4.0),
        push(4.0, 0.0),
    ];
    let top = push(2.0, 5.0);
    let crossing = points.segment_crossing(ends).unwrap();
    let segments = [
        Segment {
            ends: [ends[0], ends[1]],
        },
        Segment {
            ends: [ends[2], ends[3]],
        },
        Segment {
            ends: [crossing, top],
        },
    ];
    let mut pairs = PairBudget::new(1000, "test pairs");
    let arrangement = arrange(&mut points, &segments, &mut pairs, 0, &options()).expect("arrange");
    assert_eq!(arrangement.vertices.len(), 6);
    assert_eq!(arrangement.pieces.len(), 5);
    let found = faces(&points, &arrangement, &mut pairs, &options()).expect("faces");
    assert_eq!(found.count(), 1);
}

/// A seeded xorshift generator.
struct Xorshift(u64);

impl Xorshift {
    fn next(&mut self) -> u64 {
        self.0 ^= self.0 << 13;
        self.0 ^= self.0 >> 7;
        self.0 ^= self.0 << 17;
        self.0
    }

    fn below(&mut self, bound: u64) -> u64 {
        self.next() % bound
    }
}

#[test]
fn euler_counts_hold_for_random_soups() {
    let mut random = Xorshift(0x9E37_79B9_7F4A_7C15);
    for _ in 0..500 {
        let count = 1 + random.below(9);
        let mut lines = Vec::new();
        while lines.len() < count as usize {
            let mut at = || (random.below(7) as f64, random.below(7) as f64);
            let (a, b) = (at(), at());
            if a != b {
                lines.push((a, b));
            }
        }
        let s = scene(&lines);
        let vertices = s.arrangement.vertices.len();
        let pieces = s.arrangement.pieces.len();
        let mut parent: Vec<usize> = (0..vertices).collect();
        fn find(parent: &mut [usize], mut i: usize) -> usize {
            while parent[i] != i {
                parent[i] = parent[parent[i]];
                i = parent[i];
            }
            i
        }
        for piece in &s.arrangement.pieces {
            let (a, b) = (
                find(&mut parent, piece.ends[0].0),
                find(&mut parent, piece.ends[1].0),
            );
            parent[a] = b;
        }
        let components = (0..vertices).filter(|&v| find(&mut parent, v) == v).count();
        assert_eq!(
            vertices + s.faces.count(),
            1 + components + pieces,
            "V - E + F = 1 + C for {lines:?}"
        );
        for piece in 0..pieces {
            assert!(s.faces.left(PieceId(piece)).0 < s.faces.count());
            assert!(s.faces.right(PieceId(piece)).0 < s.faces.count());
        }
    }
}

fn nested_squares() -> Vec<Line> {
    let mut lines = square(0.0, 0.0, 10.0);
    lines.extend(square(2.0, 2.0, 2.0));
    lines
}

#[test]
fn every_nesting_decision_counts_one_pair_and_the_budget_refuses_beyond_it() {
    let mut points = ExactPoints::new();
    let segments: Vec<Segment> = nested_squares()
        .iter()
        .map(|&((ax, ay), (bx, by))| Segment {
            ends: [
                points.push_flat(Point2::new(ax, ay)).unwrap(),
                points.push_flat(Point2::new(bx, by)).unwrap(),
            ],
        })
        .collect();
    let mut pairs = PairBudget::new(1000, "test pairs");
    let arrangement = arrange(&mut points, &segments, &mut pairs, 0, &options()).unwrap();
    let before = pairs.used;
    faces(&points, &arrangement, &mut pairs, &options()).unwrap();
    assert_eq!(pairs.used - before, 1);
    let mut none = PairBudget::new(0, "test pairs");
    assert!(matches!(
        faces(&points, &arrangement, &mut none, &options()),
        Err(GeomError::BudgetExceeded { .. })
    ));
}

#[test]
fn a_cancelled_request_stops() {
    let mut points = ExactPoints::new();
    let segments: Vec<Segment> = nested_squares()
        .iter()
        .map(|&((ax, ay), (bx, by))| Segment {
            ends: [
                points.push_flat(Point2::new(ax, ay)).unwrap(),
                points.push_flat(Point2::new(bx, by)).unwrap(),
            ],
        })
        .collect();
    let mut pairs = PairBudget::new(1000, "test pairs");
    let arrangement = arrange(&mut points, &segments, &mut pairs, 0, &options()).unwrap();
    let token = CancellationToken::new();
    token.cancel();
    let cancelled = options().with_cancellation(token);
    assert!(matches!(
        faces(&points, &arrangement, &mut pairs, &cancelled),
        Err(GeomError::Cancelled)
    ));
}
