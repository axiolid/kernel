//! Numerical integration of a scalar function over a finite interval.
//!
//! - [`integrate`]: globally adaptive Gauss-Kronrod (7-point Gauss inside a
//!   15-point Kronrod rule), as QUADPACK `qag` (Piessens et al., 1983) and
//!   OCCT `math_KronrodSingleIntegration`. The panel with the largest error
//!   estimate is bisected until the total estimate meets the tolerance. No
//!   node lies on a panel end, so integrable endpoint singularities
//!   (`x^-1/2`, `ln x`) are handled by refinement toward them.
//! - [`GaussLegendre`]: a fixed `n`-point rule, exact for polynomials of
//!   degree `2n - 1`, for callers that assemble their own integrals (energy
//!   terms over spline spans) and know the integrand's degree.
//!
//! The error estimate is QUADPACK's: per panel `|K15 - G7|`, rescaled by the
//! integrand's variation and floored at the rounding level. It is an
//! estimate, not a proof, and is conservative on smooth integrands; the
//! tests check it against the true error on closed forms, including
//! endpoint singularities.

use std::cmp::Ordering;
use std::collections::BinaryHeap;

use crate::error::{finite, tolerance, NumericError, NumericResult, Status};

/// Kronrod nodes on `[0, 1]`, descending; odd indices are the Gauss nodes.
const XGK: [f64; 8] = [
    0.991_455_371_120_812_6,
    0.949_107_912_342_758_5,
    0.864_864_423_359_769_1,
    0.741_531_185_599_394_4,
    0.586_087_235_467_691_1,
    0.405_845_151_377_397_2,
    0.207_784_955_007_898_5,
    0.0,
];
/// Kronrod weights matching [`XGK`].
const WGK: [f64; 8] = [
    0.022_935_322_010_529_22,
    0.063_092_092_629_978_55,
    0.104_790_010_322_250_2,
    0.140_653_259_715_525_9,
    0.169_004_726_639_267_9,
    0.190_350_578_064_785_4,
    0.204_432_940_075_298_9,
    0.209_482_141_084_727_8,
];
/// 7-point Gauss weights for `XGK[1], XGK[3], XGK[5], XGK[7]`.
const WG: [f64; 4] = [
    0.129_484_966_168_869_7,
    0.279_705_391_489_276_7,
    0.381_830_050_505_118_9,
    0.417_959_183_673_469_4,
];

/// Accuracy request and budget for [`integrate`].
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct IntegrationOptions {
    /// Stop when the error estimate is at most this.
    pub absolute_tolerance: f64,
    /// Or at most this times `|integral|`.
    pub relative_tolerance: f64,
    /// Most panels the interval may be split into.
    pub max_panels: usize,
}

impl IntegrationOptions {
    /// Tolerances with a budget of 2000 panels (30 000 evaluations).
    pub fn new(absolute_tolerance: f64, relative_tolerance: f64) -> Self {
        Self {
            absolute_tolerance,
            relative_tolerance,
            max_panels: 2000,
        }
    }
}

/// The value of an integral and how far it can be trusted.
#[derive(Debug, Clone, Copy, PartialEq)]
#[must_use]
pub struct Integral {
    /// The integral estimate.
    pub value: f64,
    /// Estimate of `|value - exact|`.
    pub error_estimate: f64,
    /// Integrand evaluations performed.
    pub evaluations: usize,
    /// Panels in the final partition.
    pub panels: usize,
    /// Whether the tolerance was met; otherwise the estimate says how
    /// close it came.
    pub status: Status,
}

impl Integral {
    /// The integral if the tolerance was met.
    ///
    /// # Errors
    ///
    /// [`NumericError::NotConverged`] with the error estimate reached.
    pub fn converged(self) -> NumericResult<Self> {
        if self.status == Status::Converged {
            Ok(self)
        } else {
            Err(NumericError::NotConverged {
                name: "integrate",
                error_estimate: self.error_estimate,
            })
        }
    }
}

#[derive(Debug, Clone, Copy)]
struct Panel {
    a: f64,
    b: f64,
    value: f64,
    error: f64,
}

impl PartialEq for Panel {
    fn eq(&self, other: &Self) -> bool {
        self.error.total_cmp(&other.error) == Ordering::Equal
    }
}
impl Eq for Panel {}
impl PartialOrd for Panel {
    fn partial_cmp(&self, other: &Self) -> Option<Ordering> {
        Some(self.cmp(other))
    }
}
impl Ord for Panel {
    fn cmp(&self, other: &Self) -> Ordering {
        self.error.total_cmp(&other.error)
    }
}

/// Integrate `f` over `[lower, upper]` by adaptive Gauss-Kronrod (G7/K15).
///
/// `lower > upper` integrates in reverse (negated value). The integrand is
/// evaluated only strictly inside the interval.
///
/// # Errors
///
/// Refuses non-finite bounds or tolerances, negative tolerances, a
/// tolerance request below the rounding level (absolute `0` with relative
/// under `50 eps`), a zero panel budget, and any non-finite integrand
/// value.
pub fn integrate<F>(
    mut f: F,
    lower: f64,
    upper: f64,
    options: IntegrationOptions,
) -> NumericResult<Integral>
where
    F: FnMut(f64) -> f64,
{
    finite(lower, "lower bound")?;
    finite(upper, "upper bound")?;
    let abs_tol = tolerance(options.absolute_tolerance, "absolute tolerance")?;
    let rel_tol = tolerance(options.relative_tolerance, "relative tolerance")?;
    if abs_tol == 0.0 && rel_tol < 50.0 * f64::EPSILON {
        return Err(NumericError::InvalidArgument {
            name: "tolerances",
            reason: "must allow an error above rounding: absolute > 0 or relative >= 50 eps",
        });
    }
    if options.max_panels == 0 {
        return Err(NumericError::InvalidArgument {
            name: "panel budget",
            reason: "must be at least 1",
        });
    }
    if lower == upper {
        return Ok(Integral {
            value: 0.0,
            error_estimate: 0.0,
            evaluations: 0,
            panels: 0,
            status: Status::Converged,
        });
    }
    let (a, b, sign) = if lower < upper {
        (lower, upper, 1.0)
    } else {
        (upper, lower, -1.0)
    };
    if !(b - a).is_finite() {
        return Err(NumericError::InvalidArgument {
            name: "interval",
            reason: "is wider than f64 can represent",
        });
    }

    let mut evaluations = 0usize;
    let mut heap = BinaryHeap::new();
    let mut settled: Vec<Panel> = Vec::new();
    let first = kronrod(&mut f, a, b, &mut evaluations)?;
    let mut total_value = first.value;
    let mut total_error = first.error;
    heap.push(first);
    let target = |value: f64| abs_tol.max(rel_tol * value.abs());

    let mut status = Status::Converged;
    while total_error > target(total_value) {
        if heap.len() + settled.len() >= options.max_panels {
            status = Status::BudgetExhausted;
            break;
        }
        let Some(worst) = heap.pop() else {
            status = Status::RoundoffLimited;
            break;
        };
        let mid = 0.5 * worst.a + 0.5 * worst.b;
        // A panel too narrow to split further in f64 stays as it is.
        let width = worst.b - worst.a;
        if mid <= worst.a
            || mid >= worst.b
            || width <= 100.0 * f64::EPSILON * worst.a.abs().max(worst.b.abs())
        {
            settled.push(worst);
            continue;
        }
        let left = kronrod(&mut f, worst.a, mid, &mut evaluations)?;
        let right = kronrod(&mut f, mid, worst.b, &mut evaluations)?;
        total_value += left.value + right.value - worst.value;
        total_error += left.error + right.error - worst.error;
        heap.push(left);
        heap.push(right);
        if heap.len() % 64 == 0 {
            // Recompute the running sums to stop cancellation drift.
            let (v, e) = sums(heap.iter().chain(settled.iter()));
            total_value = v;
            total_error = e;
        }
    }
    let (value, error) = sums(heap.iter().chain(settled.iter()));
    if status == Status::Converged && error > target(value) {
        status = Status::RoundoffLimited;
    }
    Ok(Integral {
        value: sign * value,
        error_estimate: error,
        evaluations,
        panels: heap.len() + settled.len(),
        status,
    })
}

fn sums<'a>(panels: impl Iterator<Item = &'a Panel>) -> (f64, f64) {
    let mut value = 0.0;
    let mut compensation = 0.0;
    let mut error = 0.0;
    for p in panels {
        // Kahan summation: many small panels near a singularity.
        let y = p.value - compensation;
        let t = value + y;
        compensation = (t - value) - y;
        value = t;
        error += p.error;
    }
    (value, error)
}

fn kronrod<F: FnMut(f64) -> f64>(
    f: &mut F,
    a: f64,
    b: f64,
    evaluations: &mut usize,
) -> NumericResult<Panel> {
    let center = 0.5 * a + 0.5 * b;
    let half = 0.5 * (b - a);
    let mut eval = |x: f64| -> NumericResult<f64> {
        *evaluations += 1;
        let v = f(x);
        if v.is_finite() {
            Ok(v)
        } else {
            Err(NumericError::NonFiniteEvaluation {
                name: "integrand",
                at: Some(x),
            })
        }
    };
    let fc = eval(center)?;
    let mut kronrod = WGK[7] * fc;
    let mut gauss = WG[3] * fc;
    let mut abs_sum = WGK[7] * fc.abs();
    let mut values = [0.0f64; 15];
    values[7] = fc;
    for j in 0..7 {
        let dx = half * XGK[j];
        let f1 = eval(center - dx)?;
        let f2 = eval(center + dx)?;
        values[j] = f1;
        values[14 - j] = f2;
        kronrod += WGK[j] * (f1 + f2);
        abs_sum += WGK[j] * (f1.abs() + f2.abs());
        if j % 2 == 1 {
            gauss += WG[j / 2] * (f1 + f2);
        }
    }
    let mean = 0.5 * kronrod;
    let mut asc = WGK[7] * (fc - mean).abs();
    for j in 0..7 {
        asc += WGK[j] * ((values[j] - mean).abs() + (values[14 - j] - mean).abs());
    }
    let result = kronrod * half;
    let resabs = abs_sum * half.abs();
    let resasc = asc * half.abs();
    let mut error = ((kronrod - gauss) * half).abs();
    if resasc != 0.0 && error != 0.0 {
        error = resasc * (1.0f64).min((200.0 * error / resasc).powf(1.5));
    }
    if resabs > f64::MIN_POSITIVE / (50.0 * f64::EPSILON) {
        error = error.max(50.0 * f64::EPSILON * resabs);
    }
    Ok(Panel {
        a,
        b,
        value: result,
        error,
    })
}

/// An `n`-point Gauss-Legendre rule on `[-1, 1]`.
///
/// Exact for polynomials of degree up to `2n - 1`; it carries no error
/// estimate, so use it only where the integrand's degree is known (or use
/// [`integrate`]). Nodes come from Newton's method on `P_n` started at
/// Tricomi's asymptotic estimates.
#[derive(Debug, Clone, PartialEq)]
pub struct GaussLegendre {
    nodes: Vec<f64>,
    weights: Vec<f64>,
}

impl GaussLegendre {
    /// Most points supported; beyond it the three-term recurrence's
    /// rounding dominates and an adaptive rule is the better tool.
    pub const MAX_POINTS: usize = 512;

    /// The `n`-point rule.
    ///
    /// # Errors
    ///
    /// Refuses `n == 0` and `n > MAX_POINTS`.
    pub fn new(n: usize) -> NumericResult<Self> {
        if n == 0 || n > Self::MAX_POINTS {
            return Err(NumericError::InvalidArgument {
                name: "point count",
                reason: "must be between 1 and 512",
            });
        }
        let mut nodes = vec![0.0; n];
        let mut weights = vec![0.0; n];
        let nf = n as f64;
        for i in 0..n.div_ceil(2) {
            let mut x = (std::f64::consts::PI * (i as f64 + 0.75) / (nf + 0.5)).cos();
            let mut derivative = 0.0;
            for _ in 0..100 {
                let (p, dp) = legendre(n, x);
                derivative = dp;
                let step = p / dp;
                x -= step;
                if step.abs() <= 1e-16 * x.abs().max(1e-300) {
                    break;
                }
            }
            let (_, dp) = legendre(n, x);
            if dp.is_finite() {
                derivative = dp;
            }
            let w = 2.0 / ((1.0 - x * x) * derivative * derivative);
            nodes[i] = -x;
            nodes[n - 1 - i] = x;
            weights[i] = w;
            weights[n - 1 - i] = w;
        }
        if n % 2 == 1 {
            nodes[n / 2] = 0.0;
        }
        Ok(Self { nodes, weights })
    }

    /// Nodes on `[-1, 1]`, ascending.
    pub fn nodes(&self) -> &[f64] {
        &self.nodes
    }

    /// Weights matching [`nodes`](Self::nodes); they sum to 2.
    pub fn weights(&self) -> &[f64] {
        &self.weights
    }

    /// Apply the rule to `f` on `[lower, upper]`.
    ///
    /// # Errors
    ///
    /// Refuses non-finite bounds and any non-finite value of `f`.
    pub fn apply<F: FnMut(f64) -> f64>(
        &self,
        mut f: F,
        lower: f64,
        upper: f64,
    ) -> NumericResult<f64> {
        finite(lower, "lower bound")?;
        finite(upper, "upper bound")?;
        let center = 0.5 * lower + 0.5 * upper;
        let half = 0.5 * upper - 0.5 * lower;
        let mut sum = 0.0;
        for (x, w) in self.nodes.iter().zip(&self.weights) {
            let t = center + half * x;
            let v = f(t);
            if !v.is_finite() {
                return Err(NumericError::NonFiniteEvaluation {
                    name: "integrand",
                    at: Some(t),
                });
            }
            sum += w * v;
        }
        Ok(sum * half)
    }
}

/// `P_n(x)` and `P_n'(x)` by the three-term recurrence.
fn legendre(n: usize, x: f64) -> (f64, f64) {
    let (mut p0, mut p1) = (1.0, x);
    for k in 2..=n {
        let kf = k as f64;
        let p2 = ((2.0 * kf - 1.0) * x * p1 - (kf - 1.0) * p0) / kf;
        p0 = p1;
        p1 = p2;
    }
    if n == 0 {
        return (1.0, 0.0);
    }
    let dp = n as f64 * (x * p1 - p0) / (x * x - 1.0);
    (p1, dp)
}
