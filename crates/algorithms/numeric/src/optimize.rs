//! Minimisation: a bounded scalar minimiser and damped Gauss-Newton
//! (Levenberg-Marquardt) for nonlinear least squares.
//!
//! - [`minimize_scalar`]: Brent's `localmin` (golden section with
//!   parabolic interpolation), Brent 1973 chapter 5, OCCT `math_BrentMinimum`.
//! - [`levenberg_marquardt`]: minimises `1/2 ||r(x)||^2` with Moré's
//!   column scaling (Moré, *The Levenberg-Marquardt algorithm:
//!   implementation and theory*, 1978) and Nielsen's damping update.
//!   Each step solves the damped system as an augmented least-squares
//!   problem by QR rather than forming normal equations, so the Jacobian's
//!   conditioning is not squared. A square system with a zero-residual
//!   solution makes this a solver for nonlinear equations (OCCT
//!   `math_FunctionSetRoot`).
//!
//! Both are local methods: they find a local minimiser near the start (or
//! inside the bracket), and say how they stopped.

use crate::error::{all_finite, finite, tolerance, NumericError, NumericResult, Status};
use crate::linalg::{Matrix, Qr};

/// A local minimum of a scalar function on an interval.
#[derive(Debug, Clone, Copy, PartialEq)]
#[must_use]
pub struct ScalarMinimum {
    /// The minimiser estimate.
    pub x: f64,
    /// `f(x)`.
    pub value: f64,
    /// Lower end of the final bracket.
    pub lower: f64,
    /// Upper end of the final bracket. If `f` is unimodal on the starting
    /// interval, the minimiser lies in `[lower, upper]`.
    pub upper: f64,
    /// Iterations performed.
    pub iterations: usize,
    /// Function evaluations performed.
    pub evaluations: usize,
    /// How the search stopped.
    pub status: Status,
}

impl ScalarMinimum {
    /// A bound on `|x - minimiser|` for a unimodal function: the final
    /// bracket width.
    pub fn error_bound(&self) -> f64 {
        self.upper - self.lower
    }

    /// The minimum if the tolerance was met.
    ///
    /// # Errors
    ///
    /// [`NumericError::NotConverged`] with the bracket width reached.
    pub fn converged(self) -> NumericResult<Self> {
        if self.status == Status::Converged {
            Ok(self)
        } else {
            Err(NumericError::NotConverged {
                name: "minimize_scalar",
                error_estimate: self.error_bound(),
            })
        }
    }
}

/// Minimise `f` on `[lower, upper]` by Brent's method.
///
/// Stops when the bracket around the estimate is within
/// `2 (sqrt(eps) |x| + x_tolerance / 3)` of it on both sides. Relative
/// accuracy is limited to `sqrt(eps)`: near a smooth minimum, `f` is flat to
/// second order, so `f64` cannot place the minimiser more closely. The
/// function is never evaluated at the interval ends; for a monotone `f` the
/// estimate converges to within the tolerance of the lower end.
///
/// # Errors
///
/// Refuses non-finite bounds or tolerance, `lower >= upper`, a negative
/// tolerance, and any non-finite value of `f`.
pub fn minimize_scalar<F>(
    mut f: F,
    lower: f64,
    upper: f64,
    x_tolerance: f64,
) -> NumericResult<ScalarMinimum>
where
    F: FnMut(f64) -> f64,
{
    finite(lower, "lower bound")?;
    finite(upper, "upper bound")?;
    tolerance(x_tolerance, "x tolerance")?;
    if lower >= upper {
        return Err(NumericError::InvalidArgument {
            name: "interval",
            reason: "must have lower < upper",
        });
    }
    const MAX_ITERATIONS: usize = 1000;
    let golden = 0.5 * (3.0 - 5.0f64.sqrt());
    let sqrt_eps = f64::EPSILON.sqrt();
    let mut evaluations = 0usize;
    let mut eval = |x: f64| -> NumericResult<f64> {
        evaluations += 1;
        let v = f(x);
        if v.is_finite() {
            Ok(v)
        } else {
            Err(NumericError::NonFiniteEvaluation {
                name: "function",
                at: Some(x),
            })
        }
    };

    let (mut a, mut b) = (lower, upper);
    let mut x = a + golden * (b - a);
    let (mut w, mut v) = (x, x);
    let mut fx = eval(x)?;
    let (mut fw, mut fv) = (fx, fx);
    let (mut d, mut e) = (0.0f64, 0.0f64);
    let mut iterations = 0;
    let status = loop {
        let xm = 0.5 * (a + b);
        let tol1 = sqrt_eps * x.abs() + x_tolerance / 3.0;
        let tol2 = 2.0 * tol1;
        if (x - xm).abs() <= tol2 - 0.5 * (b - a) {
            break Status::Converged;
        }
        if iterations >= MAX_ITERATIONS {
            break Status::BudgetExhausted;
        }
        iterations += 1;
        let mut use_golden = true;
        if e.abs() > tol1 {
            let r = (x - w) * (fx - fv);
            let mut q = (x - v) * (fx - fw);
            let mut p = (x - v) * q - (x - w) * r;
            q = 2.0 * (q - r);
            if q > 0.0 {
                p = -p;
            }
            q = q.abs();
            let previous = e;
            e = d;
            if p.abs() < (0.5 * q * previous).abs() && p > q * (a - x) && p < q * (b - x) {
                d = p / q;
                let u = x + d;
                if u - a < tol2 || b - u < tol2 {
                    d = tol1.copysign(xm - x);
                }
                use_golden = false;
            }
        }
        if use_golden {
            e = if x >= xm { a - x } else { b - x };
            d = golden * e;
        }
        let u = if d.abs() >= tol1 {
            x + d
        } else {
            x + tol1.copysign(d)
        };
        let fu = eval(u)?;
        if fu <= fx {
            if u >= x {
                a = x;
            } else {
                b = x;
            }
            v = w;
            fv = fw;
            w = x;
            fw = fx;
            x = u;
            fx = fu;
        } else {
            if u < x {
                a = u;
            } else {
                b = u;
            }
            if fu <= fw || w == x {
                v = w;
                fv = fw;
                w = u;
                fw = fu;
            } else if fu <= fv || v == x || v == w {
                v = u;
                fv = fu;
            }
        }
    };
    Ok(ScalarMinimum {
        x,
        value: fx,
        lower: a,
        upper: b,
        iterations,
        evaluations,
        status,
    })
}

/// A vector-valued residual function `r: R^n -> R^m` for
/// [`levenberg_marquardt`].
pub trait Residuals {
    /// `m`, the number of residuals.
    fn residual_count(&self) -> usize;

    /// Write `r(x)` into `out` (length `m`).
    fn residuals(&mut self, x: &[f64], out: &mut [f64]);

    /// Write the Jacobian `dr_i / dx_j` into `out` (`m x n`).
    ///
    /// The default is forward differences with step
    /// `sqrt(eps) max(|x_j|, 1)`; supply the analytic Jacobian where it is
    /// available, since differences cost `n` evaluations and half the
    /// digits.
    fn jacobian(&mut self, x: &[f64], out: &mut Matrix) {
        let m = self.residual_count();
        let mut base = vec![0.0; m];
        self.residuals(x, &mut base);
        let mut shifted = vec![0.0; m];
        let mut probe = x.to_vec();
        for j in 0..x.len() {
            let h = f64::EPSILON.sqrt() * x[j].abs().max(1.0);
            probe[j] = x[j] + h;
            let actual = probe[j] - x[j];
            self.residuals(&probe, &mut shifted);
            for i in 0..m {
                out[(i, j)] = (shifted[i] - base[i]) / actual;
            }
            probe[j] = x[j];
        }
    }
}

/// Stopping rules and budget for [`levenberg_marquardt`].
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct LevenbergMarquardtOptions {
    /// Stop when `||J^T r||_inf` is at most this.
    pub gradient_tolerance: f64,
    /// Stop when a step is at most this times `||x|| + step_tolerance`.
    pub step_tolerance: f64,
    /// Stop when the cost `1/2 ||r||^2` is at most this (for zero-residual
    /// problems such as equation systems; `0` disables).
    pub cost_tolerance: f64,
    /// Most iterations (accepted or rejected steps).
    pub max_iterations: usize,
}

impl LevenbergMarquardtOptions {
    /// Gradient and step tolerance `tolerance`, no cost tolerance, at most
    /// 500 iterations.
    pub fn new(tolerance: f64) -> Self {
        Self {
            gradient_tolerance: tolerance,
            step_tolerance: tolerance,
            cost_tolerance: 0.0,
            max_iterations: 500,
        }
    }
}

/// Why [`levenberg_marquardt`] stopped.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Termination {
    /// The gradient met its tolerance.
    Gradient,
    /// A step met the step tolerance.
    Step,
    /// The cost met its tolerance.
    Cost,
    /// The iteration budget ran out.
    IterationLimit,
    /// No damping makes the cost decrease: rounding dominates the model.
    NoProgress,
}

/// A local minimiser of `1/2 ||r(x)||^2`.
#[derive(Debug, Clone, PartialEq)]
#[must_use]
pub struct LeastSquaresMinimum {
    /// The minimiser estimate.
    pub x: Vec<f64>,
    /// `1/2 ||r(x)||^2`.
    pub cost: f64,
    /// `||J^T r||_inf` at `x`.
    pub gradient_norm: f64,
    /// Length of the undamped Gauss-Newton step from `x`: near an isolated
    /// minimiser with a well-conditioned Jacobian, an estimate of the
    /// distance to it. Where `jacobian_condition` is large the iteration
    /// converges only linearly and this can underestimate the distance by a
    /// small factor. `None` where the Jacobian is rank-deficient at `x`, so
    /// the minimiser is not isolated and no estimate is given.
    pub error_estimate: Option<f64>,
    /// Estimate of the Jacobian's condition number at `x`.
    pub jacobian_condition: f64,
    /// Iterations performed.
    pub iterations: usize,
    /// Residual evaluations performed (Jacobian evaluations not counted).
    pub evaluations: usize,
    /// Why it stopped.
    pub termination: Termination,
    /// [`Status`] view of `termination`.
    pub status: Status,
}

impl LeastSquaresMinimum {
    /// The minimum if a tolerance was met.
    ///
    /// # Errors
    ///
    /// [`NumericError::NotConverged`] with the gradient norm reached.
    pub fn converged(self) -> NumericResult<Self> {
        if self.status == Status::Converged {
            Ok(self)
        } else {
            Err(NumericError::NotConverged {
                name: "levenberg_marquardt",
                error_estimate: self.error_estimate.unwrap_or(f64::INFINITY),
            })
        }
    }
}

fn norm(v: &[f64]) -> f64 {
    v.iter().map(|t| t * t).sum::<f64>().sqrt()
}

/// Minimise `1/2 ||r(x)||^2` from `x0` by Levenberg-Marquardt.
///
/// # Errors
///
/// Refuses non-finite `x0`, non-finite or negative tolerances, no unknowns
/// and no residuals. A non-finite residual or Jacobian at `x0` or at an
/// accepted point is refused by name; a non-finite residual at a trial
/// point only rejects that step. Fewer residuals than unknowns is accepted
/// (the damping keeps each step defined), but the minimiser is then not
/// isolated and `error_estimate` is `None`.
pub fn levenberg_marquardt<R: Residuals>(
    problem: &mut R,
    x0: &[f64],
    options: LevenbergMarquardtOptions,
) -> NumericResult<LeastSquaresMinimum> {
    all_finite(x0, "initial point")?;
    let gtol = tolerance(options.gradient_tolerance, "gradient tolerance")?;
    let xtol = tolerance(options.step_tolerance, "step tolerance")?;
    let ftol = tolerance(options.cost_tolerance, "cost tolerance")?;
    let n = x0.len();
    let m = problem.residual_count();
    if n == 0 {
        return Err(NumericError::InvalidArgument {
            name: "initial point",
            reason: "must have at least one unknown",
        });
    }
    if m == 0 {
        return Err(NumericError::InvalidArgument {
            name: "residual count",
            reason: "must be at least one",
        });
    }

    let mut evaluations = 0usize;
    let mut evaluate = |problem: &mut R, x: &[f64]| -> Vec<f64> {
        evaluations += 1;
        let mut r = vec![0.0; m];
        problem.residuals(x, &mut r);
        r
    };
    let jacobian = |problem: &mut R, x: &[f64]| -> NumericResult<Matrix> {
        let mut j = Matrix::zeros(m, n);
        problem.jacobian(x, &mut j);
        if (0..m).all(|i| j.row(i).iter().all(|v| v.is_finite())) {
            Ok(j)
        } else {
            Err(NumericError::NonFiniteEvaluation {
                name: "jacobian",
                at: None,
            })
        }
    };

    let mut x = x0.to_vec();
    let mut r = evaluate(problem, &x);
    if r.iter().any(|v| !v.is_finite()) {
        return Err(NumericError::NonFiniteEvaluation {
            name: "residuals",
            at: None,
        });
    }
    let mut cost = 0.5 * r.iter().map(|t| t * t).sum::<f64>();
    let mut jac = jacobian(problem, &x)?;
    let gradient = |jac: &Matrix, r: &[f64]| -> Vec<f64> {
        (0..n)
            .map(|j| (0..m).map(|i| jac[(i, j)] * r[i]).sum())
            .collect()
    };
    let mut g = gradient(&jac, &r);
    let column_norm =
        |jac: &Matrix, j: usize| (0..m).map(|i| jac[(i, j)].powi(2)).sum::<f64>().sqrt();
    let mut scale: Vec<f64> = (0..n)
        .map(|j| {
            let c = column_norm(&jac, j);
            if c > 0.0 {
                c
            } else {
                1.0
            }
        })
        .collect();
    let mut mu = 1e-3f64;
    let mut nu = 2.0f64;
    let mut iterations = 0;

    let termination = loop {
        let g_inf = g.iter().fold(0.0f64, |s, v| s.max(v.abs()));
        if g_inf <= gtol {
            break Termination::Gradient;
        }
        if cost <= ftol {
            break Termination::Cost;
        }
        if iterations >= options.max_iterations {
            break Termination::IterationLimit;
        }
        iterations += 1;

        // [J; sqrt(mu) D] delta = [-r; 0].
        let root_mu = mu.sqrt();
        let augmented = Matrix::from_fn(m + n, n, |i, j| {
            if i < m {
                jac[(i, j)]
            } else if i - m == j {
                root_mu * scale[j]
            } else {
                0.0
            }
        });
        let mut rhs: Vec<f64> = r.iter().map(|v| -v).collect();
        rhs.extend(std::iter::repeat_n(0.0, n));
        let step = Qr::new(&augmented).and_then(|qr| qr.solve_least_squares(&rhs));
        let delta = match step {
            Ok(s) => s.x,
            Err(NumericError::RankDeficient { .. }) => {
                mu *= nu;
                nu *= 2.0;
                if !mu.is_finite() || mu > 1e300 {
                    break Termination::NoProgress;
                }
                continue;
            }
            Err(other) => return Err(other),
        };
        let x_norm = norm(&x);
        if norm(&delta) <= xtol * (x_norm + xtol) {
            break Termination::Step;
        }
        let trial: Vec<f64> = x.iter().zip(&delta).map(|(a, b)| a + b).collect();
        let r_trial = evaluate(problem, &trial);
        let cost_trial = if r_trial.iter().all(|v| v.is_finite()) {
            0.5 * r_trial.iter().map(|t| t * t).sum::<f64>()
        } else {
            f64::INFINITY
        };
        // Predicted reduction of the linear model.
        let jd = jac.mul_vec(&delta)?;
        let model: f64 = r.iter().zip(&jd).map(|(a, b)| (a + b) * (a + b)).sum();
        let predicted = cost - 0.5 * model;
        let actual = cost - cost_trial;
        if predicted > 0.0 && actual > 0.0 {
            let rho = actual / predicted;
            x = trial;
            r = r_trial;
            cost = cost_trial;
            jac = jacobian(problem, &x)?;
            g = gradient(&jac, &r);
            for (j, s) in scale.iter_mut().enumerate() {
                *s = s.max(column_norm(&jac, j));
            }
            mu *= (1.0f64 / 3.0).max(1.0 - (2.0 * rho - 1.0).powi(3));
            nu = 2.0;
        } else {
            mu *= nu;
            nu *= 2.0;
            if !mu.is_finite() || mu > 1e300 {
                break Termination::NoProgress;
            }
        }
    };

    let gradient_norm = g.iter().fold(0.0f64, |s, v| s.max(v.abs()));
    let (error_estimate, jacobian_condition) = match Qr::new(&jac) {
        Ok(qr) => {
            let neg: Vec<f64> = r.iter().map(|v| -v).collect();
            let estimate = qr.solve_least_squares(&neg).ok().map(|s| norm(&s.x));
            let condition = if qr.rank() < n {
                f64::INFINITY
            } else {
                qr.condition_estimate()
            };
            (estimate, condition)
        }
        Err(_) => (None, f64::INFINITY),
    };
    let status = match termination {
        Termination::Gradient | Termination::Step | Termination::Cost => Status::Converged,
        Termination::IterationLimit => Status::BudgetExhausted,
        Termination::NoProgress => Status::RoundoffLimited,
    };
    Ok(LeastSquaresMinimum {
        x,
        cost,
        gradient_norm,
        error_estimate,
        jacobian_condition,
        iterations,
        evaluations,
        termination,
        status,
    })
}
