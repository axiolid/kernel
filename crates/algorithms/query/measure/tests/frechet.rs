//! Fréchet distance between polylines (#147, ledger row H4).
//!
//! Oracles are closed forms, plus an independent one for pseudo-random
//! curves: the discrete distance of dense resamplings, which is within the
//! sample spacing of the continuous distance.

use axiolid_core::{Point2, Point3};
use axiolid_measure::{
    discrete_frechet_distance, discrete_frechet_distance_2d, frechet_at_most, frechet_at_most_2d,
    frechet_distance, frechet_distance_2d, FrechetError,
};

fn p2(points: &[(f64, f64)]) -> Vec<Point2> {
    points.iter().map(|&(x, y)| Point2::new(x, y)).collect()
}

fn both(a: &[Point2], b: &[Point2]) -> (f64, f64) {
    (
        frechet_distance_2d(a, b).unwrap(),
        discrete_frechet_distance_2d(a, b).unwrap(),
    )
}

#[test]
fn identical_curves_are_zero_apart() {
    let a = p2(&[(0.0, 0.0), (1.0, 2.0), (3.0, -1.0), (4.0, 4.0)]);
    assert_eq!(both(&a, &a), (0.0, 0.0));
}

#[test]
fn a_translated_curve_is_the_translation_away() {
    let a = p2(&[(0.0, 0.0), (1.0, 2.0), (3.0, -1.0), (4.0, 4.0)]);
    let b: Vec<Point2> = a.iter().map(|p| *p + Point2::new(3.0, 4.0)).collect();
    assert_eq!(both(&a, &b), (5.0, 5.0));
}

#[test]
fn a_translated_3d_curve_is_the_translation_away() {
    let a = vec![
        Point3::new(0.0, 0.0, 0.0),
        Point3::new(1.0, 0.0, 2.0),
        Point3::new(1.0, 3.0, 2.0),
        Point3::new(-2.0, 3.0, 5.0),
    ];
    let b: Vec<Point3> = a.iter().map(|p| *p + Point3::new(1.0, 2.0, 2.0)).collect();
    assert_eq!(frechet_distance(&a, &b).unwrap(), 3.0);
    assert_eq!(discrete_frechet_distance(&a, &b).unwrap(), 3.0);
}

#[test]
fn parallel_segments_are_as_far_as_their_ends() {
    let a = p2(&[(0.0, 0.0), (10.0, 0.0)]);
    let b = p2(&[(0.0, 1.0), (10.0, 1.0)]);
    assert_eq!(both(&a, &b), (1.0, 1.0));
    // Shorter and shifted: the ends set it, sqrt(2^2 + 1).
    let c = p2(&[(2.0, 1.0), (8.0, 1.0)]);
    let root5 = 5.0_f64.sqrt();
    assert_eq!(both(&a, &c), (root5, root5));
}

#[test]
fn walking_along_segments_beats_hopping_between_vertices() {
    // b has a vertex half way; a does not. The walker on a stops under it
    // continuously (distance 1), but must stand on an end discretely.
    let a = p2(&[(0.0, 0.0), (2.0, 0.0)]);
    let b = p2(&[(0.0, 1.0), (1.0, 1.0), (2.0, 1.0)]);
    let (continuous, discrete) = both(&a, &b);
    assert_eq!(continuous, 1.0);
    assert_eq!(discrete, 2.0_f64.sqrt());
}

#[test]
fn direction_matters() {
    let a = p2(&[(0.0, 0.0), (10.0, 0.0)]);
    let reversed = p2(&[(10.0, 0.0), (0.0, 0.0)]);
    assert_eq!(both(&a, &a), (0.0, 0.0));
    assert_eq!(both(&a, &reversed), (10.0, 10.0));
}

#[test]
fn a_backtrack_costs_half_its_length() {
    // b runs to 6, back to 4, then on to 10. The walker on a waits at 5
    // while b backtracks: leash 1, where the two vertices 6 and 4 are
    // equidistant from a point of a's segment.
    let a = p2(&[(0.0, 0.0), (10.0, 0.0)]);
    let b = p2(&[(0.0, 0.0), (6.0, 0.0), (4.0, 0.0), (10.0, 0.0)]);
    let (continuous, discrete) = both(&a, &b);
    assert_eq!(continuous, 1.0);
    // Discretely a has only its two ends: 6 and 4 must sit on one of them.
    assert_eq!(discrete, 6.0);
    assert!(frechet_at_most_2d(&a, &b, 1.0).unwrap());
    assert!(!frechet_at_most_2d(&a, &b, 0.999_999).unwrap());
    // Either curve may be the one that backtracks.
    assert_eq!(both(&b, &a), (1.0, 6.0));
    assert!(!frechet_at_most_2d(&b, &a, 0.999_999).unwrap());
}

#[test]
fn a_detour_at_the_start_is_walked_while_the_other_waits() {
    // b climbs 3 and comes back before setting off along a. a's walker can
    // only wait at the start, so the leash is the climb.
    let a = p2(&[(0.0, 0.0), (10.0, 0.0)]);
    let b = p2(&[(0.0, 0.0), (0.0, 3.0), (0.0, 0.0), (10.0, 0.0)]);
    assert_eq!(both(&a, &b), (3.0, 3.0));
    assert_eq!(both(&b, &a), (3.0, 3.0));
    assert!(!frechet_at_most_2d(&a, &b, 2.999).unwrap());
    assert!(!frechet_at_most_2d(&b, &a, 2.999).unwrap());
}

#[test]
fn a_zig_zag_is_its_amplitude_from_the_line_beneath() {
    let h = 0.75;
    let zig = p2(&[(0.0, 0.0), (1.0, h), (2.0, 0.0), (3.0, h), (4.0, 0.0)]);
    let line = p2(&[(0.0, 0.0), (4.0, 0.0)]);
    let dense = p2(&[(0.0, 0.0), (1.0, 0.0), (2.0, 0.0), (3.0, 0.0), (4.0, 0.0)]);
    assert_eq!(frechet_distance_2d(&zig, &line).unwrap(), h);
    assert_eq!(both(&zig, &dense), (h, h));
    // With only the line's ends, the zig-zag's middle vertex (2, 0) must
    // stand on one of them, 2 away either way.
    assert_eq!(discrete_frechet_distance_2d(&zig, &line).unwrap(), 2.0);
}

#[test]
fn a_point_is_as_far_as_the_farthest_vertex() {
    let point = p2(&[(0.0, 0.0)]);
    // The farthest vertex is inside the curve, not at an end.
    let curve = p2(&[(3.0, 4.0), (6.0, -8.0), (0.0, 2.0)]);
    assert_eq!(both(&point, &curve), (10.0, 10.0));
    assert_eq!(both(&curve, &point), (10.0, 10.0));
    assert_eq!(both(&point, &p2(&[(0.0, 0.5)])), (0.5, 0.5));
}

#[test]
fn inputs_are_refused_not_guessed() {
    let a = p2(&[(0.0, 0.0), (1.0, 0.0)]);
    let empty: Vec<Point2> = Vec::new();
    let nan = p2(&[(0.0, f64::NAN), (1.0, 0.0)]);
    let inf = p2(&[(0.0, 0.0), (f64::INFINITY, 0.0)]);
    for (x, y) in [(&a, &empty), (&empty, &a)] {
        assert_eq!(frechet_distance_2d(x, y), Err(FrechetError::EmptyPolyline));
        assert_eq!(
            discrete_frechet_distance_2d(x, y),
            Err(FrechetError::EmptyPolyline)
        );
        assert_eq!(
            frechet_at_most_2d(x, y, 1.0),
            Err(FrechetError::EmptyPolyline)
        );
    }
    for bad in [&nan, &inf] {
        assert_eq!(
            frechet_distance_2d(&a, bad),
            Err(FrechetError::NonFiniteInput)
        );
        assert_eq!(
            discrete_frechet_distance_2d(bad, &a),
            Err(FrechetError::NonFiniteInput)
        );
    }
    for eps in [-1.0, f64::NAN, f64::INFINITY] {
        assert_eq!(
            frechet_at_most_2d(&a, &a, eps),
            Err(FrechetError::InvalidBound)
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

#[test]
fn random_curves_agree_with_dense_discrete_walks() {
    let step = 0.01;
    for seed in 0..12 {
        let a = walk(seed, 3 + (seed as usize % 5));
        let b = walk(seed + 100, 2 + (seed as usize % 7));
        let (continuous, discrete) = both(&a, &b);
        let dense = discrete_frechet_distance_2d(&resample(&a, step), &resample(&b, step)).unwrap();
        // Sampling only restricts the walk, and a walk between samples is
        // within half a step of one along the segments on each curve.
        assert!(
            continuous <= dense + 1e-12,
            "{seed}: {continuous} > {dense}"
        );
        assert!(
            dense <= continuous + step,
            "{seed}: {dense} vs {continuous}"
        );
        assert!(continuous <= discrete + 1e-12, "{seed}");
        // Symmetric, and the decision agrees at the answer.
        assert_eq!(frechet_distance_2d(&b, &a).unwrap(), continuous, "{seed}");
        assert!(frechet_at_most_2d(&a, &b, continuous).unwrap(), "{seed}");
        assert!(
            !frechet_at_most_2d(&a, &b, continuous * (1.0 - 1e-9)).unwrap(),
            "{seed}"
        );
    }
}

#[test]
fn the_plane_and_space_agree() {
    let a = walk(7, 6);
    let b = walk(8, 5);
    let lift = |c: &[Point2]| -> Vec<Point3> { c.iter().map(|p| p.extend(0.0)).collect() };
    assert_eq!(
        frechet_distance_2d(&a, &b).unwrap(),
        frechet_distance(&lift(&a), &lift(&b)).unwrap()
    );
    assert_eq!(
        discrete_frechet_distance_2d(&a, &b).unwrap(),
        discrete_frechet_distance(&lift(&a), &lift(&b)).unwrap()
    );
    assert_eq!(
        frechet_at_most_2d(&a, &b, 1.5).unwrap(),
        frechet_at_most(&lift(&a), &lift(&b), 1.5).unwrap()
    );
}

#[test]
fn long_curves_decide_at_their_answer() {
    // The same walk twice, one copy with its interior jittered by another
    // walk's steps: shared ends, so the answer is not the ends' distance.
    let a = walk(21, 120);
    let jitter = walk(22, 120);
    let b: Vec<Point2> = (0..120)
        .map(|i| {
            if i == 0 || i == 119 {
                a[i]
            } else {
                a[i] + (jitter[i] - jitter[i - 1]) * 0.3
            }
        })
        .collect();
    let continuous = frechet_distance_2d(&a, &b).unwrap();
    // Not settled by the endpoints alone, which return before any search.
    let ends = (a[0] - b[0]).length().max((a[119] - b[119]).length());
    assert!(continuous > ends, "{continuous} vs ends {ends}");
    assert!(frechet_at_most_2d(&a, &b, continuous).unwrap());
    assert!(!frechet_at_most_2d(&a, &b, continuous * (1.0 - 1e-9)).unwrap());
    assert!(continuous <= discrete_frechet_distance_2d(&a, &b).unwrap());
}
