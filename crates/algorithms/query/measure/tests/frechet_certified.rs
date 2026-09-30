//! Certified Fréchet decision (#147, ledger row H4).
//!
//! Oracles: comfortable cases agree with `frechet_at_most`; a case built
//! from a known-exact Fréchet distance, queried right at (or a tiny
//! fraction of an ULP from) that distance, must be `Undecided`; and a
//! property test checks that whenever the certified decision is not
//! `Undecided`, it agrees with `frechet_at_most` at the same `eps`.

use axiolid_core::Point2;
use axiolid_measure::{
    frechet_at_most_2d, frechet_decide_certified_2d, FrechetDecision, FrechetError,
};

fn p2(points: &[(f64, f64)]) -> Vec<Point2> {
    points.iter().map(|&(x, y)| Point2::new(x, y)).collect()
}

#[test]
fn comfortable_cases_agree_with_frechet_at_most() {
    let a = p2(&[(0.0, 0.0), (1.0, 2.0), (3.0, -1.0), (4.0, 4.0)]);
    let b: Vec<Point2> = a.iter().map(|p| *p + Point2::new(3.0, 4.0)).collect();
    // The exact Fréchet distance here is 5.0 (see `tests/frechet.rs`).
    assert_eq!(
        frechet_decide_certified_2d(&a, &b, 100.0).unwrap(),
        FrechetDecision::AtMost
    );
    assert!(frechet_at_most_2d(&a, &b, 100.0).unwrap());
    assert_eq!(
        frechet_decide_certified_2d(&a, &b, 0.001).unwrap(),
        FrechetDecision::MoreThan
    );
    assert!(!frechet_at_most_2d(&a, &b, 0.001).unwrap());
}

#[test]
fn an_eps_within_the_margin_of_the_true_distance_is_undecided() {
    let a = p2(&[(0.0, 0.0), (1.0, 2.0), (3.0, -1.0), (4.0, 4.0)]);
    let b: Vec<Point2> = a.iter().map(|p| *p + Point2::new(3.0, 4.0)).collect();
    let true_distance = 5.0;

    // Far from the true distance: never Undecided.
    assert_ne!(
        frechet_decide_certified_2d(&a, &b, true_distance * 2.0).unwrap(),
        FrechetDecision::Undecided
    );
    assert_ne!(
        frechet_decide_certified_2d(&a, &b, true_distance * 0.5).unwrap(),
        FrechetDecision::Undecided
    );

    // A delta far smaller than any plausible margin, built from the true
    // distance's own ULP spacing.
    let delta = 1e-13;
    let just_below = frechet_decide_certified_2d(&a, &b, true_distance - delta).unwrap();
    let just_above = frechet_decide_certified_2d(&a, &b, true_distance + delta).unwrap();
    let at = frechet_decide_certified_2d(&a, &b, true_distance).unwrap();
    assert_eq!(just_below, FrechetDecision::Undecided, "{just_below:?}");
    assert_eq!(just_above, FrechetDecision::Undecided, "{just_above:?}");
    assert_eq!(at, FrechetDecision::Undecided, "{at:?}");
}

#[test]
fn a_decision_when_not_undecided_agrees_with_frechet_at_most() {
    let a = p2(&[(0.0, 0.0), (1.0, 2.0), (3.0, -1.0), (4.0, 4.0)]);
    let b: Vec<Point2> = a.iter().map(|p| *p + Point2::new(3.0, 4.0)).collect();
    let true_distance = 5.0;
    let mut saw_at_most = false;
    let mut saw_more_than = false;
    for i in -2000..=2000 {
        let eps = true_distance + (i as f64) * 0.01;
        if eps < 0.0 {
            continue;
        }
        match frechet_decide_certified_2d(&a, &b, eps).unwrap() {
            FrechetDecision::AtMost => {
                saw_at_most = true;
                assert!(frechet_at_most_2d(&a, &b, eps).unwrap(), "eps={eps}");
            }
            FrechetDecision::MoreThan => {
                saw_more_than = true;
                assert!(!frechet_at_most_2d(&a, &b, eps).unwrap(), "eps={eps}");
            }
            FrechetDecision::Undecided => {}
        }
    }
    assert!(saw_at_most);
    assert!(saw_more_than);
}

#[test]
fn inputs_are_refused_not_guessed() {
    let a = p2(&[(0.0, 0.0), (1.0, 0.0)]);
    let empty: Vec<Point2> = Vec::new();
    let nan = p2(&[(0.0, f64::NAN), (1.0, 0.0)]);
    let inf = p2(&[(0.0, 0.0), (f64::INFINITY, 0.0)]);
    for (x, y) in [(&a, &empty), (&empty, &a)] {
        assert_eq!(
            frechet_decide_certified_2d(x, y, 1.0),
            Err(FrechetError::EmptyPolyline)
        );
    }
    for bad in [&nan, &inf] {
        assert_eq!(
            frechet_decide_certified_2d(&a, bad, 1.0),
            Err(FrechetError::NonFiniteInput)
        );
    }
    for eps in [-1.0, f64::NAN, f64::INFINITY] {
        assert_eq!(
            frechet_decide_certified_2d(&a, &a, eps),
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

#[test]
fn random_curves_agree_with_frechet_at_most_when_decided() {
    use axiolid_measure::frechet_distance_2d;
    for seed in 0..8 {
        let a = walk(seed, 3 + (seed as usize % 5));
        let b = walk(seed + 100, 2 + (seed as usize % 7));
        let distance = frechet_distance_2d(&a, &b).unwrap();
        for i in -20..=20 {
            let eps = (distance + (i as f64) * (distance.max(1.0) * 0.1)).max(0.0);
            match frechet_decide_certified_2d(&a, &b, eps).unwrap() {
                FrechetDecision::AtMost => {
                    assert!(frechet_at_most_2d(&a, &b, eps).unwrap(), "{seed} eps={eps}");
                }
                FrechetDecision::MoreThan => {
                    assert!(
                        !frechet_at_most_2d(&a, &b, eps).unwrap(),
                        "{seed} eps={eps}"
                    );
                }
                FrechetDecision::Undecided => {}
            }
        }
    }
}
