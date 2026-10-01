//! Real polynomial roots: simple, clustered and multiple roots against
//! their closed forms, completeness, rigorous evaluation bounds, refusals.

use axiolid_numeric::{NumericError, Polynomial, PolynomialRoot, RootKind};

fn roots_of(true_roots: &[f64], lower: f64, upper: f64) -> Vec<PolynomialRoot> {
    Polynomial::from_roots(true_roots)
        .expect("finite")
        .real_roots(lower, upper)
        .expect("non-zero polynomial")
}

fn contains(root: &PolynomialRoot, x: f64) -> bool {
    root.lower <= x && x <= root.upper
}

/// Every true root in the interval lies in exactly one reported entry, and
/// every simple entry contains exactly one true root.
fn assert_complete(true_roots: &[f64], lower: f64, upper: f64, found: &[PolynomialRoot]) {
    for &r in true_roots.iter().filter(|r| (lower..=upper).contains(*r)) {
        let holders = found.iter().filter(|f| contains(f, r)).count();
        assert!(holders >= 1, "root {r} not enclosed by {found:?}");
    }
    for f in found {
        assert!(
            f.lower <= f.x && f.x <= f.upper,
            "estimate outside its enclosure: {f:?}"
        );
        let inside = true_roots.iter().filter(|r| contains(f, **r)).count();
        match f.kind {
            RootKind::Simple => assert_eq!(inside, 1, "simple entry {f:?} holds {inside} roots"),
            RootKind::Unresolved { max_count } => {
                assert!(inside <= max_count, "{f:?} holds {inside} roots");
            }
        }
    }
}

#[test]
fn distinct_roots_are_simple_and_tightly_enclosed() {
    let truth = [1.0, 2.0, 3.0];
    let found = roots_of(&truth, 0.0, 4.0);
    assert_eq!(found.len(), 3);
    for (f, r) in found.iter().zip(truth) {
        assert_eq!(f.kind, RootKind::Simple);
        assert!(contains(f, r));
        assert!(f.error_bound() <= 5e-13, "{f:?}");
        assert!((f.x - r).abs() <= f.error_bound());
    }
    assert_complete(&truth, 0.0, 4.0, &found);
}

#[test]
fn roots_outside_the_interval_are_not_reported() {
    let found = roots_of(&[-5.0, 1.5, 7.0], 0.0, 4.0);
    assert_eq!(found.len(), 1);
    assert!(contains(&found[0], 1.5));
    assert!(roots_of(&[-5.0, 7.0], 0.0, 4.0).is_empty());
}

#[test]
fn no_real_roots() {
    let p = Polynomial::new(&[1.0, 0.0, 1.0]).expect("finite");
    assert!(p.real_roots(-10.0, 10.0).expect("non-zero").is_empty());
    let p = Polynomial::new(&[3.0]).expect("finite");
    assert!(p.real_roots(-10.0, 10.0).expect("non-zero").is_empty());
}

#[test]
fn wilkinson_degree_ten_roots_are_all_simple() {
    let truth: Vec<f64> = (1..=10).map(f64::from).collect();
    let found = roots_of(&truth, 0.5, 10.5);
    assert_eq!(found.len(), 10);
    assert!(found.iter().all(|f| f.kind == RootKind::Simple));
    for (f, r) in found.iter().zip(&truth) {
        assert!(contains(f, *r), "{r} not in {f:?}");
        assert!(f.error_bound() <= 1e-6, "{f:?}");
    }
    assert_complete(&truth, 0.5, 10.5, &found);
}

#[test]
fn chebyshev_roots_match_the_closed_form() {
    // T_8(x) = 128x^8 - 256x^6 + 160x^4 - 32x^2 + 1.
    let p =
        Polynomial::new(&[1.0, 0.0, -32.0, 0.0, 160.0, 0.0, -256.0, 0.0, 128.0]).expect("finite");
    let found = p.real_roots(-1.0, 1.0).expect("non-zero");
    let mut truth: Vec<f64> = (1..=8)
        .map(|k| ((2 * k - 1) as f64 * std::f64::consts::PI / 16.0).cos())
        .collect();
    truth.sort_by(f64::total_cmp);
    assert_eq!(found.len(), 8);
    for (f, r) in found.iter().zip(&truth) {
        assert_eq!(f.kind, RootKind::Simple);
        assert!((f.x - r).abs() <= 1e-14, "{} vs {r}", f.x);
    }
    assert_complete(&truth, -1.0, 1.0, &found);
}

#[test]
fn close_but_separable_cluster_is_resolved() {
    let truth = [1.0, 1.001, 1.002];
    let found = roots_of(&truth, 0.0, 2.0);
    assert_eq!(found.len(), 3, "{found:?}");
    assert!(found.iter().all(|f| f.kind == RootKind::Simple));
    assert_complete(&truth, 0.0, 2.0, &found);
}

#[test]
fn inseparable_cluster_is_one_unresolved_region() {
    // Roots 1e-9 apart: |p| between them is 2.5e-19, below the rounding
    // bound of any f64 evaluation there.
    let truth = [1.0, 1.0 + 1e-9];
    let found = roots_of(&truth, 0.0, 2.0);
    assert_eq!(found.len(), 1, "{found:?}");
    assert_eq!(found[0].kind, RootKind::Unresolved { max_count: 2 });
    assert_complete(&truth, 0.0, 2.0, &found);
}

#[test]
fn multiple_roots_are_unresolved_with_their_multiplicity_bound() {
    let found = roots_of(&[1.0, 1.0, 1.0], 0.0, 2.0);
    assert_eq!(found.len(), 1, "{found:?}");
    assert_eq!(found[0].kind, RootKind::Unresolved { max_count: 3 });
    assert!(contains(&found[0], 1.0));
    // The region is small: |p| = |x - 1|^3 exceeds rounding beyond ~1e-5.
    assert!(found[0].error_bound() < 1e-4, "{:?}", found[0]);

    let truth = [1.0, 1.0, 3.0];
    let found = roots_of(&truth, 0.0, 4.0);
    assert_eq!(found.len(), 2, "{found:?}");
    assert_eq!(found[0].kind, RootKind::Unresolved { max_count: 2 });
    assert_eq!(found[1].kind, RootKind::Simple);
    assert!(contains(&found[1], 3.0));
    assert_complete(&truth, 0.0, 4.0, &found);
}

#[test]
fn a_double_root_beside_simple_roots() {
    let truth = [-2.0, 0.5, 0.5, 4.0];
    let found = roots_of(&truth, -3.0, 5.0);
    assert_eq!(found.len(), 3, "{found:?}");
    assert_eq!(found[0].kind, RootKind::Simple);
    assert_eq!(found[1].kind, RootKind::Unresolved { max_count: 2 });
    assert_eq!(found[2].kind, RootKind::Simple);
    assert_complete(&truth, -3.0, 5.0, &found);
}

#[test]
fn a_root_at_the_interval_end_is_reported() {
    let found = roots_of(&[0.0, 1.0], 0.0, 1.0);
    assert_eq!(found.len(), 2, "{found:?}");
    assert!(contains(&found[0], 0.0) && contains(&found[1], 1.0));
    assert!(found
        .iter()
        .all(|f| f.kind == RootKind::Unresolved { max_count: 1 }));
}

/// Deterministic pseudo-random polynomials with known roots, some coincident.
#[test]
fn reported_roots_are_complete_on_random_root_sets() {
    let mut state = 0x2545_f491_4f6c_dd1du64;
    let mut next = || {
        state ^= state << 13;
        state ^= state >> 7;
        state ^= state << 17;
        (state >> 11) as f64 / (1u64 << 53) as f64
    };
    for case in 0..300 {
        let degree = 1 + case % 7;
        let mut truth: Vec<f64> = Vec::new();
        for _ in 0..degree {
            if !truth.is_empty() && next() < 0.2 {
                truth.push(truth[truth.len() - 1]);
            } else {
                truth.push((next() * 8.0 - 4.0 * 1.0).round() / 4.0 + next() * 0.01);
            }
        }
        let found = roots_of(&truth, -3.0, 3.0);
        assert_complete(&truth, -3.0, 3.0, &found);
    }
}

#[test]
fn evaluation_bound_covers_the_true_rounding_error() {
    // (x - 1)^7 expanded: catastrophic cancellation near x = 1.
    let p = Polynomial::from_roots(&[1.0; 7]).expect("finite");
    for i in 0..200 {
        let x = 0.99 + i as f64 * 1e-4;
        let (value, bound) = p.evaluate_with_bound(x);
        assert_eq!(value, p.evaluate(x));
        let exact = compensated_horner(p.coefficients(), x);
        assert!(
            (value - exact).abs() <= bound,
            "x = {x}: error {} exceeds bound {bound}",
            (value - exact).abs()
        );
    }
}

/// Horner's rule in double-double arithmetic: a reference far more accurate
/// than the plain evaluation.
fn compensated_horner(coefficients: &[f64], x: f64) -> f64 {
    let two_sum = |a: f64, b: f64| {
        let s = a + b;
        let bb = s - a;
        (s, (a - (s - bb)) + (b - bb))
    };
    let two_product = |a: f64, b: f64| {
        let p = a * b;
        (p, a.mul_add(b, -p))
    };
    let (mut hi, mut lo) = (0.0, 0.0);
    for &a in coefficients.iter().rev() {
        let (p, pe) = two_product(hi, x);
        let (s, se) = two_sum(p, a);
        lo = lo * x + pe + se;
        hi = s;
    }
    hi + lo
}

#[test]
fn derivative_and_degree() {
    let p = Polynomial::new(&[1.0, 2.0, 3.0, 0.0, 0.0]).expect("finite");
    assert_eq!(p.degree(), Some(2));
    assert_eq!(p.derivative().coefficients(), &[2.0, 6.0]);
    assert_eq!(p.evaluate(2.0), 17.0);
    assert_eq!(Polynomial::new(&[0.0]).expect("finite").degree(), None);
}

#[test]
fn refuses_bad_input_by_name() {
    assert_eq!(
        Polynomial::new(&[1.0, f64::NAN]),
        Err(NumericError::NonFiniteInput {
            name: "polynomial coefficient"
        })
    );
    let zero = Polynomial::new(&[0.0, 0.0]).expect("finite");
    assert_eq!(zero.real_roots(0.0, 1.0), Err(NumericError::ZeroPolynomial));
    let p = Polynomial::new(&[1.0, 1.0]).expect("finite");
    assert_eq!(
        p.real_roots(f64::NEG_INFINITY, 1.0),
        Err(NumericError::NonFiniteInput {
            name: "lower bound"
        })
    );
    assert!(matches!(
        p.real_roots(2.0, 1.0),
        Err(NumericError::InvalidArgument {
            name: "interval",
            ..
        })
    ));
}
