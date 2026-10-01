//! Minimisation: bounded scalar minima against closed forms, and
//! Levenberg-Marquardt on the standard test problems (Rosenbrock, Powell's
//! singular function, exponential and circle fits, an equation system),
//! with the reported estimates checked against the true error.

/// Function, bracket ends, expected answer.
type Case = (fn(f64) -> f64, f64, f64, f64);

use axiolid_numeric::{
    levenberg_marquardt, minimize_scalar, LevenbergMarquardtOptions, Matrix, NumericError,
    Residuals, Status, Termination,
};
use std::f64::consts::{E, PI};

#[test]
fn scalar_minima_lie_in_the_reported_bracket() {
    let cases: [Case; 6] = [
        (|x| (x - 2.0) * (x - 2.0) + 1.0, 0.0, 5.0, 2.0),
        (f64::cos, 0.0, 2.0 * PI, PI),
        (|x| x * x.ln(), 0.05, 1.0, 1.0 / E),
        (|x| (x - 0.3).abs(), 0.0, 1.0, 0.3),
        (|x| x.powi(4) - 3.0 * x, -2.0, 3.0, 0.75f64.cbrt()),
        (|x| (x - 1e-3).powi(2), -1.0, 1.0, 1e-3),
    ];
    for (f, a, b, exact) in cases {
        let m = minimize_scalar(f, a, b, 1e-10).expect("valid");
        assert_eq!(m.status, Status::Converged);
        assert!(m.lower <= m.x && m.x <= m.upper);
        assert!(
            (m.x - exact).abs() <= m.error_bound(),
            "minimiser {exact}: got {} with bound {}",
            m.x,
            m.error_bound()
        );
        assert!((m.x - exact).abs() <= 1e-7, "{} vs {exact}", m.x);
        assert!(m.evaluations <= 80, "{} evaluations", m.evaluations);
    }
}

#[test]
fn scalar_minimiser_tolerance_controls_the_bracket() {
    // exp(x) - 2x: minimum at ln 2, not a parabola, so the tolerance matters.
    let f = |x: f64| x.exp() - 2.0 * x;
    let coarse = minimize_scalar(f, -3.0, 5.0, 1e-2).expect("valid");
    let fine = minimize_scalar(f, -3.0, 5.0, 1e-9).expect("valid");
    assert!(coarse.error_bound() <= 4.0 * (1e-2 / 3.0 + 2.0 * f64::EPSILON.sqrt()));
    assert!((coarse.x - 2f64.ln()).abs() <= coarse.error_bound());
    assert!((fine.x - 2f64.ln()).abs() <= fine.error_bound());
    assert!(coarse.evaluations < fine.evaluations);
    assert!(fine.error_bound() < coarse.error_bound());
}

#[test]
fn scalar_minimiser_on_a_monotone_function_approaches_the_end() {
    let m = minimize_scalar(|x| x, 1.0, 2.0, 1e-8).expect("valid");
    assert!(m.x - 1.0 <= 1e-7 && m.x > 1.0);
}

#[test]
fn scalar_minimiser_refuses_bad_input_by_name() {
    assert_eq!(
        minimize_scalar(|x| x, f64::NAN, 1.0, 0.0),
        Err(NumericError::NonFiniteInput {
            name: "lower bound"
        })
    );
    assert!(matches!(
        minimize_scalar(|x| x, 1.0, 1.0, 0.0),
        Err(NumericError::InvalidArgument {
            name: "interval",
            ..
        })
    ));
    assert!(matches!(
        minimize_scalar(|x| x, 0.0, 1.0, f64::NAN),
        Err(NumericError::NonFiniteInput {
            name: "x tolerance"
        })
    ));
    assert!(matches!(
        minimize_scalar(|x| x.ln(), -1.0, 1.0, 0.0),
        Err(NumericError::NonFiniteEvaluation {
            name: "function",
            ..
        })
    ));
}

/// Rosenbrock as residuals: r = (10 (y - x^2), 1 - x).
struct Rosenbrock;
impl Residuals for Rosenbrock {
    fn residual_count(&self) -> usize {
        2
    }
    fn residuals(&mut self, x: &[f64], out: &mut [f64]) {
        out[0] = 10.0 * (x[1] - x[0] * x[0]);
        out[1] = 1.0 - x[0];
    }
    fn jacobian(&mut self, x: &[f64], out: &mut Matrix) {
        out[(0, 0)] = -20.0 * x[0];
        out[(0, 1)] = 10.0;
        out[(1, 0)] = -1.0;
        out[(1, 1)] = 0.0;
    }
}

fn distance(a: &[f64], b: &[f64]) -> f64 {
    a.iter()
        .zip(b)
        .map(|(p, q)| (p - q).powi(2))
        .sum::<f64>()
        .sqrt()
}

#[test]
fn rosenbrock_from_the_standard_start() {
    let result = levenberg_marquardt(
        &mut Rosenbrock,
        &[-1.2, 1.0],
        LevenbergMarquardtOptions::new(1e-12),
    )
    .expect("valid");
    assert_eq!(result.status, Status::Converged, "{result:?}");
    let error = distance(&result.x, &[1.0, 1.0]);
    assert!(error <= 1e-10, "{result:?}");
    let estimate = result.error_estimate.expect("isolated minimum");
    assert!(
        error <= 2.0 * estimate + 1e-15,
        "error {error} vs estimate {estimate}"
    );
    assert!(result.cost <= 1e-20);
    assert!(result.iterations <= 50, "{} iterations", result.iterations);
}

/// The same problem with the default forward-difference Jacobian.
struct RosenbrockDifferenced;
impl Residuals for RosenbrockDifferenced {
    fn residual_count(&self) -> usize {
        2
    }
    fn residuals(&mut self, x: &[f64], out: &mut [f64]) {
        Rosenbrock.residuals(x, out);
    }
}

#[test]
fn rosenbrock_with_a_differenced_jacobian() {
    let result = levenberg_marquardt(
        &mut RosenbrockDifferenced,
        &[-1.2, 1.0],
        LevenbergMarquardtOptions::new(1e-10),
    )
    .expect("valid");
    assert_eq!(result.status, Status::Converged, "{result:?}");
    assert!(distance(&result.x, &[1.0, 1.0]) <= 1e-7, "{result:?}");
}

/// Powell's singular function: the Jacobian is singular at the minimiser
/// x = 0, so convergence is only linear and the minimiser is not isolated
/// to first order.
struct PowellSingular;
impl Residuals for PowellSingular {
    fn residual_count(&self) -> usize {
        4
    }
    fn residuals(&mut self, x: &[f64], out: &mut [f64]) {
        out[0] = x[0] + 10.0 * x[1];
        out[1] = 5f64.sqrt() * (x[2] - x[3]);
        out[2] = (x[1] - 2.0 * x[2]).powi(2);
        out[3] = 10f64.sqrt() * (x[0] - x[3]).powi(2);
    }
    fn jacobian(&mut self, x: &[f64], out: &mut Matrix) {
        let s5 = 5f64.sqrt();
        let s10 = 10f64.sqrt();
        let rows = [
            [1.0, 10.0, 0.0, 0.0],
            [0.0, 0.0, s5, -s5],
            [
                0.0,
                2.0 * (x[1] - 2.0 * x[2]),
                -4.0 * (x[1] - 2.0 * x[2]),
                0.0,
            ],
            [
                2.0 * s10 * (x[0] - x[3]),
                0.0,
                0.0,
                -2.0 * s10 * (x[0] - x[3]),
            ],
        ];
        for (i, row) in rows.iter().enumerate() {
            for (j, v) in row.iter().enumerate() {
                out[(i, j)] = *v;
            }
        }
    }
}

#[test]
fn powell_singular_converges_to_the_origin() {
    let result = levenberg_marquardt(
        &mut PowellSingular,
        &[3.0, -1.0, 0.0, 1.0],
        LevenbergMarquardtOptions::new(1e-14),
    )
    .expect("valid");
    assert_eq!(result.status, Status::Converged, "{result:?}");
    assert!(distance(&result.x, &[0.0; 4]) <= 1e-3, "{result:?}");
    assert!(result.cost <= 1e-18, "{result:?}");
    // The singular Jacobian shows in the condition estimate, the signal
    // that the step-length estimate is not reliable here.
    assert!(result.jacobian_condition > 1e4, "{result:?}");
}

/// Fit y = a exp(b t) to exact samples.
struct ExponentialFit {
    t: Vec<f64>,
    y: Vec<f64>,
}
impl Residuals for ExponentialFit {
    fn residual_count(&self) -> usize {
        self.t.len()
    }
    fn residuals(&mut self, x: &[f64], out: &mut [f64]) {
        for (i, (t, y)) in self.t.iter().zip(&self.y).enumerate() {
            out[i] = x[0] * (x[1] * t).exp() - y;
        }
    }
    fn jacobian(&mut self, x: &[f64], out: &mut Matrix) {
        for (i, t) in self.t.iter().enumerate() {
            let e = (x[1] * t).exp();
            out[(i, 0)] = e;
            out[(i, 1)] = x[0] * t * e;
        }
    }
}

#[test]
fn exponential_fit_recovers_exact_parameters() {
    let t: Vec<f64> = (0..20).map(|i| i as f64 * 0.1).collect();
    let y = t.iter().map(|t| 2.5 * (-1.3 * t).exp()).collect();
    let mut problem = ExponentialFit { t, y };
    let result = levenberg_marquardt(
        &mut problem,
        &[1.0, 0.0],
        LevenbergMarquardtOptions::new(1e-13),
    )
    .expect("valid");
    assert_eq!(result.status, Status::Converged, "{result:?}");
    let error = distance(&result.x, &[2.5, -1.3]);
    assert!(error <= 1e-10, "{result:?}");
    let estimate = result.error_estimate.expect("isolated");
    assert!(
        error <= 2.0 * estimate + 1e-14,
        "error {error} vs estimate {estimate}"
    );
}

/// Geometric circle fit: r_i = |p_i - c| - R, a nonzero-residual problem.
struct CircleFit {
    points: Vec<(f64, f64)>,
}
impl Residuals for CircleFit {
    fn residual_count(&self) -> usize {
        self.points.len()
    }
    fn residuals(&mut self, x: &[f64], out: &mut [f64]) {
        for (i, (px, py)) in self.points.iter().enumerate() {
            out[i] = (px - x[0]).hypot(py - x[1]) - x[2];
        }
    }
    fn jacobian(&mut self, x: &[f64], out: &mut Matrix) {
        for (i, (px, py)) in self.points.iter().enumerate() {
            let d = (px - x[0]).hypot(py - x[1]);
            out[(i, 0)] = (x[0] - px) / d;
            out[(i, 1)] = (x[1] - py) / d;
            out[(i, 2)] = -1.0;
        }
    }
}

#[test]
fn circle_fit_with_symmetric_noise_recovers_the_circle() {
    // Eight points alternately 0.1 inside and outside a circle of radius 2
    // about (1, -1): by symmetry the best fit is that circle, radius 2.
    let points = (0..8)
        .map(|k| {
            let angle = k as f64 * PI / 4.0;
            let radius = if k % 2 == 0 { 2.1 } else { 1.9 };
            (1.0 + radius * angle.cos(), -1.0 + radius * angle.sin())
        })
        .collect();
    let mut problem = CircleFit { points };
    let result = levenberg_marquardt(
        &mut problem,
        &[0.0, 0.0, 1.0],
        LevenbergMarquardtOptions::new(1e-12),
    )
    .expect("valid");
    assert_eq!(result.status, Status::Converged, "{result:?}");
    assert!(distance(&result.x, &[1.0, -1.0, 2.0]) <= 1e-9, "{result:?}");
    // Nonzero residual: eight residuals of 0.1.
    assert!(
        (result.cost - 0.5 * 8.0 * 0.01).abs() <= 1e-12,
        "{result:?}"
    );
    assert!(result.gradient_norm <= 1e-9, "{result:?}");
}

/// x^2 + y^2 = 4 and y = x^2: a square nonlinear system.
struct CircleParabola;
impl Residuals for CircleParabola {
    fn residual_count(&self) -> usize {
        2
    }
    fn residuals(&mut self, x: &[f64], out: &mut [f64]) {
        out[0] = x[0] * x[0] + x[1] * x[1] - 4.0;
        out[1] = x[1] - x[0] * x[0];
    }
    fn jacobian(&mut self, x: &[f64], out: &mut Matrix) {
        out[(0, 0)] = 2.0 * x[0];
        out[(0, 1)] = 2.0 * x[1];
        out[(1, 0)] = -2.0 * x[0];
        out[(1, 1)] = 1.0;
    }
}

#[test]
fn solves_a_nonlinear_system() {
    let mut options = LevenbergMarquardtOptions::new(0.0);
    options.cost_tolerance = 1e-28;
    let result = levenberg_marquardt(&mut CircleParabola, &[1.0, 1.0], options).expect("valid");
    assert_eq!(result.termination, Termination::Cost, "{result:?}");
    // y^2 + y - 4 = 0, y > 0; x = sqrt(y).
    let y = (-1.0 + 17f64.sqrt()) / 2.0;
    let exact = [y.sqrt(), y];
    assert!(distance(&result.x, &exact) <= 1e-13, "{result:?}");
}

#[test]
fn reports_an_exhausted_budget() {
    let mut options = LevenbergMarquardtOptions::new(1e-14);
    options.max_iterations = 3;
    let result = levenberg_marquardt(&mut Rosenbrock, &[-1.2, 1.0], options).expect("valid");
    assert_eq!(result.termination, Termination::IterationLimit);
    assert_eq!(result.status, Status::BudgetExhausted);
    assert_eq!(result.iterations, 3);
    assert!(matches!(
        result.converged(),
        Err(NumericError::NotConverged {
            name: "levenberg_marquardt",
            ..
        })
    ));
}

struct NanAtStart;
impl Residuals for NanAtStart {
    fn residual_count(&self) -> usize {
        1
    }
    fn residuals(&mut self, x: &[f64], out: &mut [f64]) {
        out[0] = x[0].sqrt();
    }
}

#[test]
fn levenberg_marquardt_refuses_bad_input_by_name() {
    let options = LevenbergMarquardtOptions::new(1e-10);
    assert_eq!(
        levenberg_marquardt(&mut Rosenbrock, &[f64::NAN, 1.0], options),
        Err(NumericError::NonFiniteInput {
            name: "initial point"
        })
    );
    assert!(matches!(
        levenberg_marquardt(&mut Rosenbrock, &[], options),
        Err(NumericError::InvalidArgument {
            name: "initial point",
            ..
        })
    ));
    assert_eq!(
        levenberg_marquardt(&mut NanAtStart, &[-1.0], options),
        Err(NumericError::NonFiniteEvaluation {
            name: "residuals",
            at: None
        })
    );
    let mut bad = options;
    bad.gradient_tolerance = -1.0;
    assert!(matches!(
        levenberg_marquardt(&mut Rosenbrock, &[0.0, 0.0], bad),
        Err(NumericError::InvalidArgument {
            name: "gradient tolerance",
            ..
        })
    ));
}

/// Records the cost at every accepted point (the Jacobian is evaluated
/// only there).
struct Recording {
    accepted_costs: Vec<f64>,
}
impl Residuals for Recording {
    fn residual_count(&self) -> usize {
        2
    }
    fn residuals(&mut self, x: &[f64], out: &mut [f64]) {
        Rosenbrock.residuals(x, out);
    }
    fn jacobian(&mut self, x: &[f64], out: &mut Matrix) {
        let mut r = [0.0; 2];
        Rosenbrock.residuals(x, &mut r);
        self.accepted_costs.push(0.5 * (r[0] * r[0] + r[1] * r[1]));
        Rosenbrock.jacobian(x, out);
    }
}

#[test]
fn accepted_steps_never_raise_the_cost() {
    for start in [[-1.2, 1.0], [3.0, -4.0], [-2.0, 6.0], [0.0, 10.0]] {
        let mut problem = Recording {
            accepted_costs: Vec::new(),
        };
        let result =
            levenberg_marquardt(&mut problem, &start, LevenbergMarquardtOptions::new(1e-12))
                .expect("valid");
        assert_eq!(result.status, Status::Converged);
        assert!(problem.accepted_costs.len() >= 3);
        assert!(
            problem.accepted_costs.windows(2).all(|w| w[1] < w[0]),
            "start {start:?}: {:?}",
            problem.accepted_costs
        );
    }
}
