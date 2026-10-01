//! Quadrature: closed-form integrals of smooth, oscillatory and
//! endpoint-singular functions, with the reported error estimate checked
//! against the true error; fixed Gauss-Legendre rules; refusals.

use axiolid_numeric::{integrate, GaussLegendre, IntegrationOptions, NumericError, Status};
use std::f64::consts::{E, PI};

fn check(f: fn(f64) -> f64, a: f64, b: f64, exact: f64, tolerance: f64) -> usize {
    let result = integrate(f, a, b, IntegrationOptions::new(0.0, tolerance)).expect("integrable");
    assert_eq!(result.status, Status::Converged, "{result:?}");
    let true_error = (result.value - exact).abs();
    assert!(
        true_error <= result.error_estimate,
        "true error {true_error} exceeds estimate {} for exact {exact}",
        result.error_estimate
    );
    assert!(
        result.error_estimate <= tolerance * exact.abs(),
        "{result:?}"
    );
    result.evaluations
}

#[test]
fn smooth_integrands_converge_on_the_first_panels() {
    let evaluations = check(f64::exp, 0.0, 1.0, E - 1.0, 1e-12);
    assert!(evaluations <= 45, "{evaluations} evaluations");
    check(f64::sin, 0.0, PI, 2.0, 1e-12);
    check(|x| 1.0 / (1.0 + x * x), -1.0, 1.0, PI / 2.0, 1e-12);
    check(
        |x| x.powi(5) - 3.0 * x * x,
        -2.0,
        3.0,
        (3f64.powi(6) - 64.0) / 6.0 - 35.0,
        1e-13,
    );
}

#[test]
fn oscillatory_and_peaked_integrands() {
    check(|x| x.sin().powi(2), 0.0, 20.0 * PI, 10.0 * PI, 1e-10);
    // Runge-like peak: integral of 1/(1 + 1e4 x^2) over [-1, 1] = 2 atan(100)/100.
    check(
        |x| 1.0 / (1.0 + 1e4 * x * x),
        -1.0,
        1.0,
        2.0 * 100f64.atan() / 100.0,
        1e-10,
    );
}

#[test]
fn endpoint_singularities_are_refined_toward() {
    check(|x| 1.0 / x.sqrt(), 0.0, 1.0, 2.0, 1e-9);
    check(f64::ln, 0.0, 1.0, -1.0, 1e-9);
    check(|x| x.powf(-0.75), 0.0, 1.0, 4.0, 1e-8);
    check(|x| (1.0 - x * x).sqrt(), -1.0, 1.0, PI / 2.0, 1e-10);
    check(|x| x.sqrt() * x.ln(), 0.0, 1.0, -4.0 / 9.0, 1e-10);
}

#[test]
fn estimate_bounds_the_error_even_when_the_budget_runs_out() {
    let mut options = IntegrationOptions::new(0.0, 1e-12);
    options.max_panels = 4;
    let result = integrate(|x| x.powf(-0.9), 0.0, 1.0, options).expect("integrable");
    assert_eq!(result.status, Status::BudgetExhausted);
    assert!(result.panels <= 4);
    assert!((result.value - 10.0).abs() <= result.error_estimate);
    assert!(matches!(
        result.converged(),
        Err(NumericError::NotConverged {
            name: "integrate",
            ..
        })
    ));
}

#[test]
fn absolute_tolerance_and_orientation() {
    let options = IntegrationOptions::new(1e-10, 0.0);
    let forward = integrate(f64::cos, 0.0, 2.0, options).expect("smooth");
    let backward = integrate(f64::cos, 2.0, 0.0, options).expect("smooth");
    assert_eq!(forward.value, -backward.value);
    assert!((forward.value - 2f64.sin()).abs() <= forward.error_estimate);
    assert!(forward.error_estimate <= 1e-10);
    let empty = integrate(f64::cos, 1.0, 1.0, options).expect("empty");
    assert_eq!((empty.value, empty.error_estimate), (0.0, 0.0));
    // Integral exactly zero: only the absolute tolerance can be met.
    let zero = integrate(f64::sin, -1.0, 1.0, options).expect("odd");
    assert_eq!(zero.status, Status::Converged);
    assert!(zero.value.abs() <= zero.error_estimate.max(1e-16));
}

#[test]
fn gauss_legendre_matches_the_eight_point_constants() {
    let rule = GaussLegendre::new(8).expect("valid");
    let positive = [
        0.183_434_642_495_649_8,
        0.525_532_409_916_329,
        0.796_666_477_413_626_7,
        0.960_289_856_497_536_3,
    ];
    let weights = [
        0.362_683_783_378_362,
        0.313_706_645_877_887_3,
        0.222_381_034_453_374_5,
        0.101_228_536_290_376_3,
    ];
    for k in 0..4 {
        assert!((rule.nodes()[4 + k] - positive[k]).abs() <= 1e-15);
        assert!((rule.nodes()[3 - k] + positive[k]).abs() <= 1e-15);
        assert!((rule.weights()[4 + k] - weights[k]).abs() <= 1e-15);
    }
}

#[test]
fn gauss_legendre_is_exact_to_degree_2n_minus_1() {
    for n in [1usize, 2, 3, 5, 8, 13, 20, 64] {
        let rule = GaussLegendre::new(n).expect("valid");
        let sum: f64 = rule.weights().iter().sum();
        assert!((sum - 2.0).abs() <= 1e-13, "n = {n}: weights sum {sum}");
        assert!(rule.nodes().windows(2).all(|w| w[0] < w[1]));
        for degree in 0..2 * n {
            // integral of x^d over [0, 2] = 2^(d+1) / (d+1).
            let got = rule
                .apply(|x| x.powi(degree as i32), 0.0, 2.0)
                .expect("finite");
            let exact = 2f64.powi(degree as i32 + 1) / (degree as f64 + 1.0);
            assert!(
                (got - exact).abs() <= 1e-13 * exact,
                "n = {n}, degree {degree}: {got} vs {exact}"
            );
        }
    }
    // One degree higher is not exact.
    let rule = GaussLegendre::new(3).expect("valid");
    let got = rule.apply(|x| x.powi(6), -1.0, 1.0).expect("finite");
    assert!((got - 2.0 / 7.0).abs() > 1e-3);
}

#[test]
fn refuses_bad_input_by_name() {
    let options = IntegrationOptions::new(1e-10, 1e-10);
    assert_eq!(
        integrate(f64::sin, f64::NAN, 1.0, options),
        Err(NumericError::NonFiniteInput {
            name: "lower bound"
        })
    );
    assert_eq!(
        integrate(f64::sin, 0.0, f64::INFINITY, options),
        Err(NumericError::NonFiniteInput {
            name: "upper bound"
        })
    );
    assert!(matches!(
        integrate(f64::sin, 0.0, 1.0, IntegrationOptions::new(-1.0, 0.0)),
        Err(NumericError::InvalidArgument {
            name: "absolute tolerance",
            ..
        })
    ));
    assert!(matches!(
        integrate(f64::sin, 0.0, 1.0, IntegrationOptions::new(0.0, 1e-17)),
        Err(NumericError::InvalidArgument {
            name: "tolerances",
            ..
        })
    ));
    // 1/x is evaluated at the panel centre 0.
    assert_eq!(
        integrate(|x| 1.0 / x, -1.0, 1.0, options),
        Err(NumericError::NonFiniteEvaluation {
            name: "integrand",
            at: Some(0.0)
        })
    );
    assert!(GaussLegendre::new(0).is_err());
    assert!(GaussLegendre::new(513).is_err());
    assert!(GaussLegendre::new(512).is_ok());
}
