//! Bracketed root finding: closed forms, the bracket as a proven bound,
//! the iteration count on hard cases, and refusals.

/// Function, bracket ends, expected answer.
type Case = (fn(f64) -> f64, f64, f64, f64);
/// Function, bracket ends, iteration budget.
type Budgeted = (fn(f64) -> f64, f64, f64, usize);

use axiolid_numeric::{find_root, NumericError, Status};

#[test]
fn finds_closed_form_roots_inside_the_reported_bracket() {
    let cases: [Case; 5] = [
        (f64::cos, 0.0, 3.0, std::f64::consts::FRAC_PI_2),
        // Brent's own example.
        (
            |x| x * x * x - 2.0 * x - 5.0,
            2.0,
            3.0,
            2.094_551_481_542_326_5,
        ),
        (|x| x * x - 2.0, 0.0, 2.0, std::f64::consts::SQRT_2),
        (|x| x.exp() - 10.0, 0.0, 5.0, 10f64.ln()),
        (|x| (x - 1e-3).powi(3), -1.0, 1.0, 1e-3),
    ];
    for (f, a, b, exact) in cases {
        let root = find_root(f, a, b, 0.0).expect("bracketed");
        assert_eq!(root.status, Status::Converged);
        assert!(root.lower <= root.x && root.x <= root.upper);
        assert!(
            root.lower <= exact && exact <= root.upper
                || (root.x - exact).abs() <= 2.0 * f64::EPSILON * exact.abs(),
            "root {exact} outside [{}, {}]",
            root.lower,
            root.upper
        );
        // Full precision was asked for: the bracket is a few ulps wide.
        assert!(
            root.error_bound() <= 4.0 * f64::EPSILON * exact.abs() + 1e-300,
            "bracket {} too wide for {exact}",
            root.error_bound()
        );
        // f changes sign across the bracket.
        assert!(f(root.lower) * f(root.upper) <= 0.0);
    }
}

#[test]
fn honours_an_absolute_tolerance() {
    let root = find_root(|x| x * x - 2.0, 0.0, 2.0, 1e-3).expect("bracketed");
    assert!(root.error_bound() <= 1e-3 + 4.0 * f64::EPSILON * 2.0);
    assert!((root.x - std::f64::consts::SQRT_2).abs() <= root.error_bound());
    // Fewer evaluations than the full-precision search.
    let full = find_root(|x| x * x - 2.0, 0.0, 2.0, 0.0).expect("bracketed");
    assert!(root.evaluations < full.evaluations);
}

#[test]
fn an_exact_zero_at_an_end_is_returned_with_a_zero_width_bracket() {
    let root = find_root(|x| x - 1.0, 1.0, 3.0, 0.0).expect("zero at lower end");
    assert_eq!(
        (root.x, root.lower, root.upper, root.value),
        (1.0, 1.0, 1.0, 0.0)
    );
    let root = find_root(|x| x - 3.0, 1.0, 3.0, 0.0).expect("zero at upper end");
    assert_eq!(root.x, 3.0);
    assert_eq!(root.error_bound(), 0.0);
}

/// Functions on which interpolation crawls: flat near the root, or a root
/// far below the bracket's scale. Brent's acceptance test falls back to
/// bisection, so the counts stay near the bisection count (the last case
/// needs about 1100 halvings to walk the exponent range down to 1e-100).
#[test]
fn interpolation_safeguard_bounds_the_iteration_count() {
    let hard: [Budgeted; 5] = [
        (|x| x.powi(25), -1.0, 4.0, 200),
        (|x| (x - 0.3).powi(19), -1.0, 3.0, 200),
        (
            |x| (x - 0.3).signum() * (x - 0.3).abs().powf(20.0),
            0.0,
            1.0,
            200,
        ),
        (
            |x| (x - 0.3).signum() * (x - 0.3).abs().powf(0.001),
            0.0,
            1.0,
            100,
        ),
        (|x| x * x * x - 1e-300, -1.0, 1e10, 1200),
    ];
    for (f, a, b, budget) in hard {
        let root = find_root(f, a, b, 0.0).expect("bracketed");
        assert_eq!(root.status, Status::Converged);
        assert!(f(root.lower) * f(root.upper) <= 0.0);
        assert!(
            root.iterations <= budget,
            "{} iterations on [{a}, {b}]",
            root.iterations
        );
    }
    let root = find_root(|x| x * x * x - 1e-300, -1.0, 1e10, 0.0).expect("bracketed");
    assert!((root.x - 1e-100).abs() <= 1e-114);
}

#[test]
fn refuses_bad_input_by_name() {
    assert_eq!(
        find_root(|x| x, f64::NAN, 1.0, 0.0),
        Err(NumericError::NonFiniteInput {
            name: "lower bound"
        })
    );
    assert_eq!(
        find_root(|x| x, -1.0, f64::INFINITY, 0.0),
        Err(NumericError::NonFiniteInput {
            name: "upper bound"
        })
    );
    assert!(matches!(
        find_root(|x| x, -1.0, 1.0, -1.0),
        Err(NumericError::InvalidArgument {
            name: "x tolerance",
            ..
        })
    ));
    assert!(matches!(
        find_root(|x| x, 1.0, -1.0, 0.0),
        Err(NumericError::InvalidArgument {
            name: "bracket",
            ..
        })
    ));
    assert_eq!(
        find_root(|x| x * x + 1.0, -1.0, 1.0, 0.0),
        Err(NumericError::NoSignChange {
            f_lower: 2.0,
            f_upper: 2.0
        })
    );
    assert!(matches!(
        find_root(|x| if x > 0.4 { f64::NAN } else { x - 0.5 }, 0.0, 1.0, 0.0),
        Err(NumericError::NonFiniteEvaluation {
            name: "function",
            ..
        })
    ));
}
