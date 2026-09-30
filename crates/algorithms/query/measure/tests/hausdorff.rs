//! Hausdorff distance between polylines (#147, ledger row H4).
//!
//! Oracles are closed forms, plus an independent one for pseudo-random
//! curves: a dense resampling's discrete max-of-min, which brackets the
//! continuous one-sided distance within the sampling step (the same
//! strategy as `tests/frechet.rs`).

use axiolid_core::{Point2, Point3};
use axiolid_measure::{
    one_sided_polyline_hausdorff_distance, one_sided_polyline_hausdorff_distance_2d,
    polyline_hausdorff_distance_2d, FrechetError,
};

fn p2(points: &[(f64, f64)]) -> Vec<Point2> {
    points.iter().map(|&(x, y)| Point2::new(x, y)).collect()
}

#[test]
fn identical_curves_are_zero_apart() {
    let a = p2(&[(0.0, 0.0), (1.0, 2.0), (3.0, -1.0), (4.0, 4.0)]);
    assert_eq!(
        one_sided_polyline_hausdorff_distance_2d(&a, &a).unwrap(),
        0.0
    );
    assert_eq!(polyline_hausdorff_distance_2d(&a, &a).unwrap(), 0.0);
}

#[test]
fn parallel_segments_of_equal_length_are_their_offset_apart() {
    let a = p2(&[(0.0, 0.0), (10.0, 0.0)]);
    let b = p2(&[(0.0, 3.0), (10.0, 3.0)]);
    assert_eq!(
        one_sided_polyline_hausdorff_distance_2d(&a, &b).unwrap(),
        3.0
    );
    assert_eq!(
        one_sided_polyline_hausdorff_distance_2d(&b, &a).unwrap(),
        3.0
    );
    assert_eq!(polyline_hausdorff_distance_2d(&a, &b).unwrap(), 3.0);
}

#[test]
fn a_shadowed_segment_is_asymmetric() {
    // The short segment sits directly beneath the middle of the long one,
    // offset by d. Every point of the short segment has a point of the
    // long one directly across (distance d), but the long segment's
    // uncovered ends are farther than d from the short one.
    let d = 2.0;
    let short = p2(&[(4.0, d), (6.0, d)]);
    let long = p2(&[(0.0, 0.0), (10.0, 0.0)]);
    let short_to_long = one_sided_polyline_hausdorff_distance_2d(&short, &long).unwrap();
    let long_to_short = one_sided_polyline_hausdorff_distance_2d(&long, &short).unwrap();
    assert_eq!(short_to_long, d);
    // The long curve's far end (0, 0) is distance sqrt(4^2 + d^2) from the
    // short segment's near vertex (4, d), which is its closest point.
    let expected_long_to_short = (4.0_f64 * 4.0 + d * d).sqrt();
    assert!(
        (long_to_short - expected_long_to_short).abs() < 1e-9,
        "{long_to_short} vs {expected_long_to_short}"
    );
    assert!(long_to_short > short_to_long);
    assert_eq!(
        polyline_hausdorff_distance_2d(&short, &long).unwrap(),
        long_to_short
    );
}

#[test]
fn a_point_to_segment_is_ordinary_point_to_segment_distance() {
    let point = p2(&[(0.0, 0.0)]);
    let segment = p2(&[(3.0, 4.0), (3.0, -4.0)]);
    // The nearest point on the segment to the origin is (3, 0): distance 3.
    assert_eq!(
        one_sided_polyline_hausdorff_distance_2d(&point, &segment).unwrap(),
        3.0
    );
}

#[test]
fn a_segment_to_a_point_is_the_farthest_vertex() {
    // Distance-to-a-single-point is convex along a segment, so its
    // supremum sits at an endpoint, not the segment's midpoint (which is
    // closer to the point than one of its ends).
    let point = p2(&[(0.0, 0.0)]);
    let segment = p2(&[(1.0, 0.0), (5.0, 0.0)]);
    let expected = 5.0; // farthest vertex (5, 0), not the midpoint (3, 0)
    assert_eq!(
        one_sided_polyline_hausdorff_distance_2d(&segment, &point).unwrap(),
        expected
    );
}

#[test]
fn a_zig_zag_differs_from_its_frechet_distance() {
    // A zig-zag against the straight line beneath it: every vertex of the
    // zig-zag is either on the line (y = 0, matching a line vertex or
    // interior point) or at height h, so the one-sided distance from
    // zig-zag to line is exactly h (its own Frechet distance too, in this
    // case -- see `a_zig_zag_is_its_amplitude_from_the_line_beneath` in
    // `tests/frechet.rs`).
    let h = 0.75;
    let zig = p2(&[(0.0, 0.0), (1.0, h), (2.0, 0.0), (3.0, h), (4.0, 0.0)]);
    let line = p2(&[(0.0, 0.0), (4.0, 0.0)]);
    assert_eq!(
        one_sided_polyline_hausdorff_distance_2d(&zig, &line).unwrap(),
        h
    );
    // The reverse direction is where Hausdorff and Frechet genuinely part
    // ways: the line's point directly under a zig-zag peak, say (1, 0), is
    // NOT on the zig-zag, and its nearest zig-zag point is the foot on the
    // rising segment (0,0)-(1,h), not a vertex. With h = 3/4, that foot is
    // at parameter 16/25 along the segment, giving point (16/25, 12/25) and
    // distance 3/5 from (1, 0) -- a genuine segment-tie/perpendicular-foot
    // candidate, not a vertex-vertex distance, and it is the line's
    // farthest point from the zig-zag (worked by hand and confirmed with a
    // dense brute-force scan offline).
    let expected = 0.6;
    let got = one_sided_polyline_hausdorff_distance_2d(&line, &zig).unwrap();
    assert!((got - expected).abs() < 1e-9, "{got} vs {expected}");
    assert_ne!(got, h, "Hausdorff and Frechet must differ here");
}

#[test]
fn a_segment_vs_segment_tie_candidate_does_not_spoil_the_answer() {
    // Two unit-length B segments meeting at the origin, perpendicular: one
    // along +x from (0,0) to (1,0), the other along +y from (0,0) to
    // (0,1). An A segment from (-1, 2) to (2, -1) crosses the line x = y
    // at its own midpoint (0.5, 0.5), which is exactly where the two B
    // segments' infinite lines are equidistant from A's own line -- a
    // genuine line-vs-line tie candidate that the algorithm must evaluate
    // without it corrupting the true answer. The true supremum here sits
    // at A's own endpoints instead: (-1, 2) and (2, -1) are each
    // sqrt(2) from B (nearest points (0,1) and (1,0) respectively), worked
    // by hand and confirmed with a dense brute-force scan offline.
    let a = p2(&[(-1.0, 2.0), (2.0, -1.0)]);
    let b = p2(&[(0.0, 1.0), (0.0, 0.0), (1.0, 0.0)]);
    let d = one_sided_polyline_hausdorff_distance_2d(&a, &b).unwrap();
    let expected = 2.0_f64.sqrt();
    assert!((d - expected).abs() < 1e-9, "{d} vs {expected}");
}

#[test]
fn a_vertex_vs_line_tie_sets_the_answer() {
    // B: a vertex V = (5, -3), reached by a detour (V - Q - P1) that stays
    // clear of the region under test, then the segment P1-P2 = (0,0)-(10,0)
    // ("the line"). Q = (0, -10) is chosen so the detour segment V-Q's
    // *unclamped* line does not coincide with V's own distance function
    // near A (its foot stays outside [0, 1] throughout A, confirmed with a
    // dense scan offline), so this case cannot be found through a
    // line-vs-line tie in disguise: it genuinely needs the vertex feature.
    // A is the vertical segment x = 5 from y = -2.5 to y = 1.0. Along A,
    // distance to the line is |y| (its foot at x = 5 is always active);
    // distance to V is |y + 3|. These tie at y = -1.5, value 1.5 -- an
    // interior local maximum of the pointwise minimum, since the distance
    // to V is decreasing there (moving away from -3) while the distance to
    // the line is increasing, and each is smaller than 1.5 immediately on
    // either side. Neither of A's own endpoints (values 0.5 and 1.0)
    // reaches 1.5. V is not an endpoint of the line segment, so this is a
    // genuine vertex-vs-line tie, not a same-segment endpoint clamp.
    let b = p2(&[(5.0, -3.0), (0.0, -10.0), (0.0, 0.0), (10.0, 0.0)]);
    let a = p2(&[(5.0, -2.5), (5.0, 1.0)]);
    let d = one_sided_polyline_hausdorff_distance_2d(&a, &b).unwrap();
    assert!((d - 1.5).abs() < 1e-9, "{d}");
}

#[test]
fn a_line_vs_line_tie_sets_the_answer() {
    // B: the segment (0,0)-(10,0) ("s1"), a detour to (20,-6) that stays
    // clear of the region under test, then the segment (2,-6)-(8,-2)
    // ("s2"). A is the vertical segment x = 5 from y = -6.5 to y = 0.5.
    // Along A, s1's foot at x = 5 is always active (t = 0.5), giving
    // distance |y|; s2's foot is also active throughout the tie region.
    // Solving `y^2 == distance_to_s2_line(y)^2` gives the tie at
    // `y = 9 - 3 sqrt(13)`, value `3 sqrt(13) - 9 ~ 1.8167` -- an interior
    // local maximum (worked with `sympy` and confirmed with a dense scan
    // offline), higher than either of A's own endpoints (each 0.5). Both
    // feet are strictly interior there, so this is a genuine line-vs-line
    // tie, not reducible to any vertex feature.
    let b = p2(&[
        (0.0, 0.0),
        (10.0, 0.0),
        (20.0, -6.0),
        (2.0, -6.0),
        (8.0, -2.0),
    ]);
    let a = p2(&[(5.0, -6.5), (5.0, 0.5)]);
    let d = one_sided_polyline_hausdorff_distance_2d(&a, &b).unwrap();
    let expected = 3.0 * 13.0_f64.sqrt() - 9.0;
    assert!((d - expected).abs() < 1e-9, "{d} vs {expected}");
}

#[test]
fn the_plane_and_space_agree() {
    let a = p2(&[(0.0, 0.0), (2.0, 3.0), (5.0, -1.0)]);
    let b = p2(&[(1.0, 1.0), (3.0, 4.0)]);
    let lift = |c: &[Point2]| -> Vec<Point3> { c.iter().map(|p| p.extend(0.0)).collect() };
    assert_eq!(
        one_sided_polyline_hausdorff_distance_2d(&a, &b).unwrap(),
        one_sided_polyline_hausdorff_distance(&lift(&a), &lift(&b)).unwrap()
    );
    assert_eq!(
        polyline_hausdorff_distance_2d(&a, &b).unwrap(),
        axiolid_measure::polyline_hausdorff_distance(&lift(&a), &lift(&b)).unwrap()
    );
}

#[test]
fn inputs_are_refused_not_guessed() {
    let a = p2(&[(0.0, 0.0), (1.0, 0.0)]);
    let empty: Vec<Point2> = Vec::new();
    let nan = p2(&[(0.0, f64::NAN), (1.0, 0.0)]);
    let inf = p2(&[(0.0, 0.0), (f64::INFINITY, 0.0)]);
    for (x, y) in [(&a, &empty), (&empty, &a)] {
        assert_eq!(
            one_sided_polyline_hausdorff_distance_2d(x, y),
            Err(FrechetError::EmptyPolyline)
        );
        assert_eq!(
            polyline_hausdorff_distance_2d(x, y),
            Err(FrechetError::EmptyPolyline)
        );
    }
    for bad in [&nan, &inf] {
        assert_eq!(
            one_sided_polyline_hausdorff_distance_2d(&a, bad),
            Err(FrechetError::NonFiniteInput)
        );
        assert_eq!(
            polyline_hausdorff_distance_2d(bad, &a),
            Err(FrechetError::NonFiniteInput)
        );
    }
}

/// A deterministic pseudo-random walk of `n` points.
fn walk(seed: u64, n: usize) -> Vec<Point2> {
    let mut state = seed;
    let mut next = || {
        state = state
            .wrapping_mul(6364136223846793005)
            .wrapping_add(1442695040888963407);
        ((state >> 11) as f64 / (1u64 << 53) as f64) * 2.0 - 1.0
    };
    let mut p = Point2::ZERO;
    (0..n)
        .map(|_| {
            p += Point2::new(1.0 + next(), 2.0 * next());
            p
        })
        .collect()
}

/// Points along `curve` no more than `step` apart, keeping its vertices.
fn resample(curve: &[Point2], step: f64) -> Vec<Point2> {
    let mut out = vec![curve[0]];
    for s in curve.windows(2) {
        let pieces = ((s[1] - s[0]).length() / step).ceil().max(1.0) as usize;
        for k in 1..=pieces {
            out.push(s[0] + (s[1] - s[0]) * (k as f64 / pieces as f64));
        }
    }
    out
}

/// Dense-resampling oracle: max over samples of `a` of min over samples of
/// `b` of the Euclidean distance.
fn dense_one_sided(a: &[Point2], b: &[Point2], step: f64) -> f64 {
    let sa = resample(a, step);
    let sb = resample(b, step);
    sa.iter()
        .map(|&p| {
            sb.iter()
                .map(|&q| (p - q).length())
                .fold(f64::INFINITY, f64::min)
        })
        .fold(f64::NEG_INFINITY, f64::max)
}

#[test]
fn random_curves_agree_with_dense_resampling() {
    let step = 0.01;
    for seed in 0..12 {
        let a = walk(seed, 3 + (seed as usize % 5));
        let b = walk(seed + 100, 2 + (seed as usize % 7));
        let exact = one_sided_polyline_hausdorff_distance_2d(&a, &b).unwrap();
        let dense = dense_one_sided(&a, &b, step);
        // The exact supremum is the true one; the dense oracle only ever
        // finds a lower-or-equal value (it samples a subset of points), but
        // it can also overshoot slightly because its own nearest sample of
        // b is farther than the true nearest point of b's segments.
        assert!(exact >= dense - step, "{seed}: {exact} < {dense} - {step}");
        assert!(exact <= dense + step, "{seed}: {exact} > {dense} + {step}");
    }
}
