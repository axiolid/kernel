//! The Bentley–Ottmann sweep against closed-form cases and an independent
//! brute-force exact oracle (#146).
//!
//! The oracle tests every pair of segments in dyadic arithmetic, collects
//! every endpoint and proper crossing, and applies the reporting rule of
//! `segment_sweep.rs` directly, so it shares no code with the sweep.

use std::cmp::Ordering;

use axiolid_core::Point2;
use axiolid_exact::{Arith, Dyadic};
use axiolid_guarantees::Sign;
use axiolid_overlay::{
    segment_intersections, ExactPoint2, Incidence, SegmentIntersections, SegmentLocation,
    SegmentSweepError,
};

// ------------------------------------------------------------- helpers

fn p(x: f64, y: f64) -> Point2 {
    Point2::new(x, y)
}

fn s(ax: f64, ay: f64, bx: f64, by: f64) -> [Point2; 2] {
    [p(ax, ay), p(bx, by)]
}

fn run(segments: &[[Point2; 2]]) -> SegmentIntersections {
    segment_intersections(segments).expect("finite input")
}

fn inc(segment: usize, location: SegmentLocation) -> Incidence {
    Incidence { segment, location }
}

use SegmentLocation::{Degenerate, End, Interior, Start};

// -------------------------------------------------------------- oracle

/// `(x, y, w)` with `w > 0`.
#[derive(Debug, Clone)]
struct Q(Dyadic, Dyadic, Dyadic);

fn d(v: f64) -> Dyadic {
    Dyadic::from_f64(v)
}

fn sg(v: &Dyadic) -> i8 {
    match v.sign().unwrap() {
        Sign::Positive => 1,
        Sign::Negative => -1,
        _ => 0,
    }
}

fn q(point: Point2) -> Q {
    Q(d(point.x), d(point.y), d(1.0))
}

fn qcmp(a: &Q, b: &Q) -> Ordering {
    let x = sg(&a.0.mul(&b.2).sub(&b.0.mul(&a.2)));
    let y = sg(&a.1.mul(&b.2).sub(&b.1.mul(&a.2)));
    x.cmp(&0).then(y.cmp(&0))
}

fn orient(a: Point2, b: Point2, c: &Q) -> i8 {
    let (ax, ay) = (d(a.x), d(a.y));
    let v = d(b.x)
        .sub(&ax)
        .mul(&c.1.sub(&ay.mul(&c.2)))
        .sub(&d(b.y).sub(&ay).mul(&c.0.sub(&ax.mul(&c.2))));
    sg(&v)
}

fn cross(s: &[Point2; 2], t: &[Point2; 2]) -> Dyadic {
    let ux = d(s[1].x).sub(&d(s[0].x));
    let uy = d(s[1].y).sub(&d(s[0].y));
    let vx = d(t[1].x).sub(&d(t[0].x));
    let vy = d(t[1].y).sub(&d(t[0].y));
    ux.mul(&vy).sub(&uy.mul(&vx))
}

fn lex(a: Point2, b: Point2) -> Ordering {
    a.x.partial_cmp(&b.x)
        .unwrap()
        .then(a.y.partial_cmp(&b.y).unwrap())
}

fn degenerate(s: &[Point2; 2]) -> bool {
    s[0] == s[1]
}

fn contains(s: &[Point2; 2], c: &Q) -> bool {
    if degenerate(s) {
        return qcmp(c, &q(s[0])) == Ordering::Equal;
    }
    let (lo, hi) = if lex(s[0], s[1]) == Ordering::Less {
        (s[0], s[1])
    } else {
        (s[1], s[0])
    };
    orient(s[0], s[1], c) == 0
        && qcmp(&q(lo), c) != Ordering::Greater
        && qcmp(c, &q(hi)) != Ordering::Greater
}

/// The crossing when it is interior to both segments.
fn proper(s: &[Point2; 2], t: &[Point2; 2]) -> Option<Q> {
    if degenerate(s) || degenerate(t) {
        return None;
    }
    let den = cross(s, t);
    let ds = sg(&den);
    if ds == 0 {
        return None;
    }
    let acx = d(t[0].x).sub(&d(s[0].x));
    let acy = d(t[0].y).sub(&d(s[0].y));
    let ux = d(s[1].x).sub(&d(s[0].x));
    let uy = d(s[1].y).sub(&d(s[0].y));
    let vx = d(t[1].x).sub(&d(t[0].x));
    let vy = d(t[1].y).sub(&d(t[0].y));
    let num_s = acx.mul(&vy).sub(&acy.mul(&vx));
    let num_t = acx.mul(&uy).sub(&acy.mul(&ux));
    // 0 < num / den < 1 for both parameters.
    let inside = |num: &Dyadic| {
        let (n, m) = if ds > 0 {
            (num.clone(), den.sub(num))
        } else {
            (num.neg(), num.sub(&den))
        };
        sg(&n) > 0 && sg(&m) > 0
    };
    if !inside(&num_s) || !inside(&num_t) {
        return None;
    }
    let x = d(s[0].x).mul(&den).add(&num_s.mul(&ux));
    let y = d(s[0].y).mul(&den).add(&num_s.mul(&uy));
    Some(if ds > 0 {
        Q(x, y, den)
    } else {
        Q(x.neg(), y.neg(), den.neg())
    })
}

struct Expected {
    points: Vec<(Q, Vec<Incidence>)>,
    overlaps: Vec<(Point2, Point2, Vec<usize>)>,
}

fn oracle(segments: &[[Point2; 2]]) -> Expected {
    let n = segments.len();
    let mut candidates: Vec<Q> = Vec::new();
    for seg in segments {
        candidates.push(q(seg[0]));
        candidates.push(q(seg[1]));
    }
    for i in 0..n {
        for j in i + 1..n {
            if let Some(c) = proper(&segments[i], &segments[j]) {
                candidates.push(c);
            }
        }
    }
    candidates.sort_by(qcmp);
    candidates.dedup_by(|a, b| qcmp(a, b) == Ordering::Equal);

    let mut points = Vec::new();
    for c in candidates {
        let mut incidences = Vec::new();
        for (index, seg) in segments.iter().enumerate() {
            if !contains(seg, &c) {
                continue;
            }
            let location = if degenerate(seg) {
                Degenerate
            } else if qcmp(&c, &q(seg[0])) == Ordering::Equal {
                Start
            } else if qcmp(&c, &q(seg[1])) == Ordering::Equal {
                End
            } else {
                Interior
            };
            incidences.push(inc(index, location));
        }
        if incidences.len() < 2 {
            continue;
        }
        let all_interior = incidences.iter().all(|i| i.location == Interior);
        let first = &segments[incidences[0].segment];
        let collinear = incidences
            .iter()
            .all(|i| sg(&cross(first, &segments[i.segment])) == 0);
        if all_interior && collinear {
            continue;
        }
        points.push((c, incidences));
    }

    // Overlaps: group segments by carrier line, then cut each line at every
    // endpoint on it.
    let mut group: Vec<usize> = (0..n).collect();
    fn find(group: &mut [usize], i: usize) -> usize {
        let mut r = i;
        while group[r] != r {
            r = group[r];
        }
        group[i] = r;
        r
    }
    for i in 0..n {
        for j in i + 1..n {
            let (a, b) = (&segments[i], &segments[j]);
            if degenerate(a) || degenerate(b) {
                continue;
            }
            if sg(&cross(a, b)) == 0 && orient(a[0], a[1], &q(b[0])) == 0 {
                let (ri, rj) = (find(&mut group, i), find(&mut group, j));
                group[ri] = rj;
            }
        }
    }
    let mut overlaps = Vec::new();
    for root in 0..n {
        let members: Vec<usize> = (0..n)
            .filter(|&i| !degenerate(&segments[i]) && find(&mut group, i) == root)
            .collect();
        if members.len() < 2 {
            continue;
        }
        let mut ends: Vec<Point2> = members.iter().flat_map(|&i| segments[i]).collect();
        ends.sort_by(|a, b| lex(*a, *b));
        ends.dedup_by(|a, b| lex(*a, *b) == Ordering::Equal);
        let mut pieces: Vec<(Point2, Point2, Vec<usize>)> = Vec::new();
        for w in ends.windows(2) {
            let (u, v) = (w[0], w[1]);
            let covering: Vec<usize> = members
                .iter()
                .copied()
                .filter(|&i| {
                    let seg = segments[i];
                    let (lo, hi) = if lex(seg[0], seg[1]) == Ordering::Less {
                        (seg[0], seg[1])
                    } else {
                        (seg[1], seg[0])
                    };
                    lex(lo, u) != Ordering::Greater && lex(v, hi) != Ordering::Greater
                })
                .collect();
            if covering.len() < 2 {
                continue;
            }
            if let Some(last) = pieces.last_mut() {
                if last.1 == u && last.2 == covering {
                    last.1 = v;
                    continue;
                }
            }
            pieces.push((u, v, covering));
        }
        overlaps.extend(pieces);
    }
    Expected { points, overlaps }
}

fn overlap_key(o: &(Point2, Point2, Vec<usize>)) -> (u64, u64, u64, u64, Vec<usize>) {
    // Sort key only; exact equality is checked on the values themselves.
    (
        o.0.x.to_bits(),
        o.0.y.to_bits(),
        o.1.x.to_bits(),
        o.1.y.to_bits(),
        o.2.clone(),
    )
}

/// The rounded point is the nearest double to the exact one.
fn check_rounding(exact: &ExactPoint2, rounded: Point2) {
    let (x, y, w) = exact.homogeneous();
    for (n, r) in [(x, rounded.x), (y, rounded.y)] {
        let half = d(0.5);
        let low = d(r.next_down()).add(&d(r)).mul(&half).mul(&w);
        let high = d(r).add(&d(r.next_up())).mul(&half).mul(&w);
        assert!(
            sg(&n.sub(&low)) >= 0 && sg(&high.sub(&n)) >= 0,
            "coordinate {r} is not the nearest double"
        );
    }
}

fn compare(segments: &[[Point2; 2]], label: &str) -> SegmentIntersections {
    let got = run(segments);
    let want = oracle(segments);
    assert_eq!(
        got.points.len(),
        want.points.len(),
        "{label}: point count; got {:?}",
        got.points
            .iter()
            .map(|p| (p.rounded, p.incidences.clone()))
            .collect::<Vec<_>>()
    );
    for (index, (have, (point, incidences))) in got.points.iter().zip(&want.points).enumerate() {
        let (x, y, w) = have.point.homogeneous();
        assert_eq!(
            qcmp(&Q(x, y, w), point),
            Ordering::Equal,
            "{label}: point {index} at {:?}",
            have.rounded
        );
        assert_eq!(
            &have.incidences, incidences,
            "{label}: incidences at {:?}",
            have.rounded
        );
        assert_eq!(have.rounded, have.point.rounded());
        check_rounding(&have.point, have.rounded);
        if let Some(input) = have.point.as_input() {
            assert_eq!(
                input, have.rounded,
                "{label}: input points come back exactly"
            );
        }
    }
    let mut got_overlaps: Vec<_> = got
        .overlaps
        .iter()
        .map(|o| (o.from, o.to, o.segments.clone()))
        .collect();
    let mut want_overlaps = want.overlaps;
    got_overlaps.sort_by_key(overlap_key);
    want_overlaps.sort_by_key(overlap_key);
    assert_eq!(got_overlaps, want_overlaps, "{label}: overlaps");
    got
}

// ------------------------------------------------------------ generators

struct Rng(u64);

impl Rng {
    fn next(&mut self) -> u64 {
        self.0 ^= self.0 << 13;
        self.0 ^= self.0 >> 7;
        self.0 ^= self.0 << 17;
        self.0
    }

    fn below(&mut self, n: u64) -> u64 {
        self.next() % n
    }

    fn unit(&mut self) -> f64 {
        (self.next() >> 11) as f64 / (1u64 << 53) as f64
    }
}

fn grid_segments(rng: &mut Rng, n: usize, size: u64, zero_length: bool) -> Vec<[Point2; 2]> {
    (0..n)
        .map(|_| {
            let a = p(rng.below(size) as f64, rng.below(size) as f64);
            let b = if zero_length && rng.below(8) == 0 {
                a
            } else {
                p(rng.below(size) as f64, rng.below(size) as f64)
            };
            [a, b]
        })
        .collect()
}

fn scaled(segments: &[[Point2; 2]], factor: f64) -> Vec<[Point2; 2]> {
    segments
        .iter()
        .map(|[a, b]| [p(a.x * factor, a.y * factor), p(b.x * factor, b.y * factor)])
        .collect()
}

// ------------------------------------------------------ closed-form cases

#[test]
fn empty_and_single_inputs_report_nothing() {
    let got = run(&[]);
    assert!(got.points.is_empty() && got.overlaps.is_empty());
    let got = run(&[s(0.0, 0.0, 1.0, 1.0)]);
    assert!(got.points.is_empty() && got.overlaps.is_empty());
    assert_eq!(got.evidence.segments, 1);
    assert_eq!(got.evidence.events, 2);
}

#[test]
fn an_x_crosses_at_its_centre() {
    let got = run(&[s(0.0, 0.0, 2.0, 2.0), s(0.0, 2.0, 2.0, 0.0)]);
    assert_eq!(got.points.len(), 1);
    assert_eq!(got.points[0].rounded, p(1.0, 1.0));
    assert_eq!(
        got.points[0].incidences,
        vec![inc(0, Interior), inc(1, Interior)]
    );
    assert_eq!(got.evidence.crossings, 1);
}

#[test]
fn a_crossing_at_a_third_is_exact_and_correctly_rounded() {
    // y = 2x meets x + y = 1 at (1/3, 2/3), which no double holds.
    let got = run(&[s(0.0, 0.0, 1.0, 2.0), s(0.0, 1.0, 1.0, 0.0)]);
    assert_eq!(got.points.len(), 1);
    let point = &got.points[0];
    assert!(point.point.as_input().is_none());
    let (x, y, w) = point.point.homogeneous();
    assert_eq!(sg(&x.mul(&d(3.0)).sub(&w)), 0, "x is exactly 1/3");
    assert_eq!(
        sg(&y.mul(&d(3.0)).sub(&w.mul(&d(2.0)))),
        0,
        "y is exactly 2/3"
    );
    // IEEE division is correctly rounded, so it is the reference.
    assert_eq!(point.rounded, p(1.0 / 3.0, 2.0 / 3.0));
}

#[test]
fn a_vertical_segment_crosses_a_horizontal_one() {
    let got = run(&[s(1.0, -1.0, 1.0, 1.0), s(0.0, 0.0, 2.0, 0.0)]);
    assert_eq!(got.points.len(), 1);
    assert_eq!(got.points[0].rounded, p(1.0, 0.0));
    assert_eq!(
        got.points[0].incidences,
        vec![inc(0, Interior), inc(1, Interior)]
    );
}

#[test]
fn vertical_segments_sharing_ends_and_overlapping() {
    let segments = [
        s(0.0, 0.0, 0.0, 2.0),
        s(0.0, 2.0, 0.0, 3.0),
        s(0.0, 1.0, 0.0, 4.0),
        s(-1.0, 1.5, 1.0, 1.5),
    ];
    let got = compare(&segments, "vertical");
    let pieces: Vec<_> = got
        .overlaps
        .iter()
        .map(|o| (o.from, o.to, o.segments.clone()))
        .collect();
    assert_eq!(
        pieces,
        vec![
            (p(0.0, 1.0), p(0.0, 2.0), vec![0, 2]),
            (p(0.0, 2.0), p(0.0, 3.0), vec![1, 2]),
        ]
    );
}

#[test]
fn lines_crossing_beyond_a_segment_schedule_nothing() {
    // The vertical's endpoints straddle the horizontal's line, but the
    // lines meet at (3, 0), past the horizontal's end.
    let got = run(&[s(0.0, 0.0, 1.0, 0.0), s(3.0, -1.0, 3.0, 1.0)]);
    assert!(got.points.is_empty());
    assert_eq!(got.evidence.crossings, 0);
    assert_eq!(got.evidence.events, 4);
}

#[test]
fn a_t_junction_is_an_end_on_an_interior() {
    let got = run(&[s(0.0, 0.0, 2.0, 0.0), s(1.0, 1.0, 1.0, 0.0)]);
    assert_eq!(got.points.len(), 1);
    assert_eq!(
        got.points[0].incidences,
        vec![inc(0, Interior), inc(1, End)]
    );
    assert_eq!(
        got.evidence.crossings, 0,
        "an endpoint is not a crossing event"
    );
}

#[test]
fn shared_endpoints_are_reported_with_their_ends() {
    // A V, and a chain continuing it.
    let got = run(&[
        s(0.0, 1.0, 1.0, 0.0),
        s(2.0, 1.0, 1.0, 0.0),
        s(2.0, 1.0, 3.0, 1.0),
    ]);
    let found: Vec<_> = got
        .points
        .iter()
        .map(|p| (p.rounded, p.incidences.clone()))
        .collect();
    assert_eq!(
        found,
        vec![
            (p(1.0, 0.0), vec![inc(0, End), inc(1, End)]),
            (p(2.0, 1.0), vec![inc(1, Start), inc(2, Start)]),
        ]
    );
}

#[test]
fn collinear_overlap_is_a_piece_with_reported_ends() {
    let got = run(&[s(0.0, 0.0, 4.0, 2.0), s(6.0, 3.0, 2.0, 1.0)]);
    let pieces: Vec<_> = got
        .overlaps
        .iter()
        .map(|o| (o.from, o.to, o.segments.clone()))
        .collect();
    assert_eq!(pieces, vec![(p(2.0, 1.0), p(4.0, 2.0), vec![0, 1])]);
    let found: Vec<_> = got
        .points
        .iter()
        .map(|p| (p.rounded, p.incidences.clone()))
        .collect();
    assert_eq!(
        found,
        vec![
            (p(2.0, 1.0), vec![inc(0, Interior), inc(1, End)]),
            (p(4.0, 2.0), vec![inc(0, End), inc(1, Interior)]),
        ]
    );
}

#[test]
fn collinear_segments_touching_end_to_end_do_not_overlap() {
    let got = run(&[s(0.0, 0.0, 1.0, 1.0), s(1.0, 1.0, 2.0, 2.0)]);
    assert!(got.overlaps.is_empty());
    assert_eq!(got.points.len(), 1);
    assert_eq!(got.points[0].incidences, vec![inc(0, End), inc(1, Start)]);
}

#[test]
fn a_crossing_inside_an_overlap_keeps_one_piece() {
    // Two segments share [1, 3] on y = 0, and a third crosses the shared
    // piece at x = 1/3 + 2 (non-dyadic): the piece stays whole and the
    // crossing is reported once with all three segments.
    let segments = [
        s(0.0, 0.0, 3.0, 0.0),
        s(1.0, 0.0, 4.0, 0.0),
        s(2.0, -1.0, 3.0, 2.0),
    ];
    let got = compare(&segments, "crossing in overlap");
    assert_eq!(got.overlaps.len(), 1);
    assert_eq!(got.overlaps[0].from, p(1.0, 0.0));
    assert_eq!(got.overlaps[0].to, p(3.0, 0.0));
    let crossing = got
        .points
        .iter()
        .find(|p| p.point.as_input().is_none())
        .expect("the crossing");
    assert_eq!(
        crossing.incidences,
        vec![inc(0, Interior), inc(1, Interior), inc(2, Interior)]
    );
}

#[test]
fn nested_overlaps_split_where_the_covering_set_changes() {
    let segments = [
        s(0.0, 0.0, 10.0, 0.0),
        s(2.0, 0.0, 8.0, 0.0),
        s(4.0, 0.0, 6.0, 0.0),
        s(10.0, 0.0, 0.0, 0.0),
    ];
    let got = compare(&segments, "nested");
    let pieces: Vec<_> = got
        .overlaps
        .iter()
        .map(|o| (o.from.x, o.to.x, o.segments.clone()))
        .collect();
    assert_eq!(
        pieces,
        vec![
            (0.0, 2.0, vec![0, 3]),
            (2.0, 4.0, vec![0, 1, 3]),
            (4.0, 6.0, vec![0, 1, 2, 3]),
            (6.0, 8.0, vec![0, 1, 3]),
            (8.0, 10.0, vec![0, 3]),
        ]
    );
    // The duplicate's ends are reported, reversed.
    assert_eq!(got.points[0].incidences, vec![inc(0, Start), inc(3, End)]);
}

#[test]
fn many_segments_through_one_point() {
    let spokes: Vec<[Point2; 2]> = (0..12)
        .map(|k| {
            let t = f64::from(k);
            s(-t - 1.0, 12.0 - 2.0 * t, t + 1.0, 2.0 * t - 12.0)
        })
        .collect();
    let got = compare(&spokes, "star");
    let centre: Vec<_> = got.points.iter().filter(|p| p.rounded == p0()).collect();
    assert_eq!(centre.len(), 1);
    assert_eq!(centre[0].incidences.len(), 12);
}

fn p0() -> Point2 {
    p(0.0, 0.0)
}

#[test]
fn zero_length_segments_are_points() {
    let segments = [
        s(0.0, 0.0, 2.0, 2.0),
        s(1.0, 1.0, 1.0, 1.0), // on segment 0
        s(5.0, 5.0, 5.0, 5.0), // alone
        s(7.0, 0.0, 7.0, 0.0), // with its duplicate
        s(7.0, 0.0, 7.0, 0.0),
        s(2.0, 2.0, 2.0, 2.0), // on an end
    ];
    let got = compare(&segments, "zero length");
    let found: Vec<_> = got
        .points
        .iter()
        .map(|p| (p.rounded, p.incidences.clone()))
        .collect();
    assert_eq!(
        found,
        vec![
            (p(1.0, 1.0), vec![inc(0, Interior), inc(1, Degenerate)]),
            (p(2.0, 2.0), vec![inc(0, End), inc(5, Degenerate)]),
            (p(7.0, 0.0), vec![inc(3, Degenerate), inc(4, Degenerate)]),
        ]
    );
    assert!(got.overlaps.is_empty());
}

#[test]
fn non_finite_input_is_refused_by_index() {
    let err = segment_intersections(&[s(0.0, 0.0, 1.0, 1.0), s(0.0, f64::NAN, 1.0, 0.0)]);
    assert_eq!(err, Err(SegmentSweepError::NonFinite { segment: 1 }));
    let err = segment_intersections(&[s(f64::INFINITY, 0.0, 1.0, 1.0)]);
    assert_eq!(err, Err(SegmentSweepError::NonFinite { segment: 0 }));
}

#[test]
fn negative_zero_is_the_same_point() {
    let got = run(&[s(-0.0, 0.0, 1.0, 1.0), s(0.0, -0.0, 1.0, -1.0)]);
    assert_eq!(got.points.len(), 1);
    assert_eq!(got.points[0].incidences, vec![inc(0, Start), inc(1, Start)]);
}

#[test]
fn exact_points_order_lexicographically() {
    let third = ExactPoint2::from_homogeneous(d(1.0), d(5.0), d(3.0)).unwrap();
    let same = ExactPoint2::from_homogeneous(d(-2.0), d(-10.0), d(-6.0)).unwrap();
    assert_eq!(third, same);
    let below = ExactPoint2::from_homogeneous(d(1.0), d(4.0), d(3.0)).unwrap();
    assert!(below < third);
    assert!(ExactPoint2::from_homogeneous(d(1.0), d(1.0), d(0.0)).is_none());
    let (_, _, w) = same.homogeneous();
    assert_eq!(sg(&w), 1, "the weight is kept positive");
}

// ---------------------------------------------------- oracle comparisons

#[test]
fn random_grid_segments_match_the_oracle() {
    // A small grid forces shared ends, collinear overlaps, verticals,
    // T-junctions and several segments through one point.
    let mut rng = Rng(0x5eed_1234_abcd_ef01);
    for round in 0..300 {
        let n = 2 + rng.below(30) as usize;
        let size = 2 + rng.below(7);
        let segments = grid_segments(&mut rng, n, size, true);
        compare(&segments, &format!("grid round {round}"));
    }
}

#[test]
fn random_float_segments_match_the_oracle() {
    let mut rng = Rng(0x0123_4567_89ab_cdef);
    for round in 0..40 {
        let n = 10 + rng.below(60) as usize;
        let segments: Vec<[Point2; 2]> = (0..n)
            .map(|_| {
                [
                    p(rng.unit() * 2.0 - 1.0, rng.unit() * 2.0 - 1.0),
                    p(rng.unit() * 2.0 - 1.0, rng.unit() * 2.0 - 1.0),
                ]
            })
            .collect();
        let got = compare(&segments, &format!("float round {round}"));
        assert!(got.evidence.crossings > 0);
    }
}

#[test]
fn extreme_scales_match_the_oracle() {
    // Scaling by a power of two is exact, so the answer must be the same
    // scene; at these magnitudes the interval filter overflows or
    // underflows and every decision is exact.
    let mut rng = Rng(0xfeed_face_cafe_beef);
    for round in 0..30 {
        let base = grid_segments(&mut rng, 12, 6, false);
        for factor in [
            2f64.powi(-1070),
            2f64.powi(-600),
            2f64.powi(600),
            2f64.powi(1015),
        ] {
            let segments = scaled(&base, factor);
            let got = compare(&segments, &format!("scale {factor:e} round {round}"));
            let reference = run(&base);
            assert_eq!(got.points.len(), reference.points.len());
            assert_eq!(got.overlaps.len(), reference.overlaps.len());
        }
    }
}

#[test]
fn grids_of_lines_match_the_oracle() {
    let mut segments = Vec::new();
    for i in 0..8 {
        let t = f64::from(i);
        segments.push(s(0.0, t, 7.0, t)); // horizontal
        segments.push(s(t, 0.0, t, 7.0)); // vertical
        segments.push(s(0.0, t, 7.0 - t, 7.0)); // diagonal
        segments.push(s(t, 0.0, 0.0, t)); // anti-diagonal
    }
    // Overlapping copies of some grid lines.
    segments.push(s(1.0, 3.0, 5.0, 3.0));
    segments.push(s(3.0, 2.0, 3.0, 9.0));
    let got = compare(&segments, "line grid");
    assert!(got.points.len() > 64);
    assert!(!got.overlaps.is_empty());
}

#[test]
fn stars_through_a_non_dyadic_point_match_the_oracle() {
    // Q -> 1 - 2Q passes through (1/3, 1/3) in its interior for every
    // integer Q, so the whole star meets at a point no double holds.
    let mut rng = Rng(0xabad_cafe_0000_0001);
    for round in 0..20 {
        let mut segments = Vec::new();
        for _ in 0..(3 + rng.below(12)) {
            let (i, j) = (rng.below(9) as f64 - 4.0, rng.below(9) as f64 - 4.0);
            if i == 0.0 && j == 0.0 {
                continue;
            }
            segments.push(s(i, j, 1.0 - 2.0 * i, 1.0 - 2.0 * j));
        }
        let got = compare(&segments, &format!("third star round {round}"));
        let centre = got.points.iter().find(|p| {
            let (x, _, w) = p.point.homogeneous();
            sg(&x.mul(&d(3.0)).sub(&w)) == 0 && p.point.as_input().is_none()
        });
        if segments.len() >= 2 {
            assert!(centre.is_some(), "round {round}: the star centre");
        }
    }
}

#[test]
fn stars_through_an_input_point_match_the_oracle() {
    let mut rng = Rng(0x1357_9bdf_2468_ace0);
    for round in 0..30 {
        let mut segments = Vec::new();
        for _ in 0..(2 + rng.below(16)) {
            let (dx, dy) = (rng.below(7) as f64 - 3.0, rng.below(7) as f64 - 3.0);
            let (a, b) = (rng.below(3) as f64, rng.below(3) as f64);
            // Some spokes end at the centre, some pass through it.
            segments.push(s(
                2.0 - a * dx,
                2.0 - a * dy,
                2.0 + (b + 1.0) * dx,
                2.0 + (b + 1.0) * dy,
            ));
        }
        compare(&segments, &format!("input star round {round}"));
    }
}

#[test]
fn collinear_chains_with_transversals_match_the_oracle() {
    let mut rng = Rng(0x2468_1357_0f0f_f0f0);
    for round in 0..40 {
        let mut segments = Vec::new();
        // Pieces of the line y = 2x + 1, chained, overlapping and nested.
        for _ in 0..(3 + rng.below(10)) {
            let a = rng.below(12) as f64 - 6.0;
            let b = rng.below(12) as f64 - 6.0;
            if rng.below(2) == 0 {
                segments.push(s(a, 2.0 * a + 1.0, b, 2.0 * b + 1.0));
            } else {
                segments.push(s(b, 2.0 * b + 1.0, a, 2.0 * a + 1.0));
            }
        }
        // And of the vertical x = 1.
        for _ in 0..rng.below(4) {
            let a = rng.below(8) as f64;
            segments.push(s(1.0, a, 1.0, a + 1.0 + rng.below(4) as f64));
        }
        // Transversals, some through chain vertices.
        for _ in 0..(1 + rng.below(5)) {
            let x = rng.below(12) as f64 - 6.0;
            let y = rng.below(20) as f64 - 10.0;
            segments.push(s(x - 3.0, y + 1.0, x + 2.0, y - 7.0));
        }
        compare(&segments, &format!("chain round {round}"));
    }
}

#[test]
fn many_short_segments_match_the_oracle() {
    // A denser, larger scene, closer to the sweep's intended use.
    let mut rng = Rng(0x7777_8888_9999_aaaa);
    let segments: Vec<[Point2; 2]> = (0..400)
        .map(|_| {
            let (x, y) = (rng.below(64) as f64, rng.below(64) as f64);
            let (dx, dy) = (rng.below(9) as f64 - 4.0, rng.below(9) as f64 - 4.0);
            s(x, y, x + dx, y + dy)
        })
        .collect();
    let got = compare(&segments, "many short");
    assert!(got.points.len() > 100);
}
