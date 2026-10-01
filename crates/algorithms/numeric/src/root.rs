//! Bracketed scalar root finding with guaranteed convergence.
//!
//! [`find_root`] is Brent's method (inverse quadratic interpolation, the
//! secant step, and bisection), as in Brent, *Algorithms for Minimization
//! without Derivatives* (1973), chapter 4, and OCCT `math_BracketedRoot`.
//! An interpolation step is accepted only if it is less than half the step
//! before last and stays inside the bracket; otherwise the step bisects.
//! Convergence is therefore guaranteed for any continuous function with a
//! sign change, in at most about the square of the bisection count and in
//! practice far fewer.
//!
//! The answer is a bracket, not a point: the function changes sign across
//! `[lower, upper]` (or vanishes exactly at `x`), so for a continuous
//! function a root lies inside and `|x - root| <= upper - lower`.

use crate::error::{finite, tolerance, NumericError, NumericResult, Status};

/// A root located inside a sign-change bracket.
#[derive(Debug, Clone, Copy, PartialEq)]
#[must_use]
pub struct Root {
    /// Best estimate: the bracket end with the smaller `|f|`.
    pub x: f64,
    /// `f(x)`.
    pub value: f64,
    /// Lower end of the final bracket.
    pub lower: f64,
    /// Upper end of the final bracket. The function changes sign across
    /// `[lower, upper]`, or `lower == upper == x` with `value == 0`.
    pub upper: f64,
    /// Iterations performed.
    pub iterations: usize,
    /// Function evaluations performed.
    pub evaluations: usize,
    /// [`Status::Converged`], or [`Status::BudgetExhausted`] if the
    /// iteration cap was hit first; the bracket is valid either way.
    pub status: Status,
}

impl Root {
    /// A bound on `|x - root|`: the final bracket width.
    pub fn error_bound(&self) -> f64 {
        self.upper - self.lower
    }
}

/// Iteration cap. Bisection from any finite bracket reaches the stopping
/// width in under 2200 halvings, and Brent's method rarely needs more steps
/// than bisection.
const MAX_ITERATIONS: usize = 10_000;

/// Find a root of `f` in `[lower, upper]` by Brent's method.
///
/// `f(lower)` and `f(upper)` must differ in sign (or one of them be zero).
/// Iteration stops when the bracket is no wider than
/// `4 * eps * |x| + x_tolerance`, or when its ends are adjacent floats;
/// `x_tolerance = 0` therefore asks for full precision.
///
/// # Errors
///
/// Refuses non-finite bounds or tolerance, `lower > upper`, a negative
/// tolerance, bracket ends of the same sign
/// ([`NumericError::NoSignChange`]), and any non-finite value of `f`.
pub fn find_root<F>(mut f: F, lower: f64, upper: f64, x_tolerance: f64) -> NumericResult<Root>
where
    F: FnMut(f64) -> f64,
{
    finite(lower, "lower bound")?;
    finite(upper, "upper bound")?;
    tolerance(x_tolerance, "x tolerance")?;
    if lower > upper {
        return Err(NumericError::InvalidArgument {
            name: "bracket",
            reason: "must have lower <= upper",
        });
    }
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
    let mut fa = eval(a)?;
    let mut fb = eval(b)?;
    if fa == 0.0 {
        return Ok(exact(a, 0, evaluations));
    }
    if fb == 0.0 {
        return Ok(exact(b, 0, evaluations));
    }
    if fa.signum() == fb.signum() {
        return Err(NumericError::NoSignChange {
            f_lower: fa,
            f_upper: fb,
        });
    }

    let (mut c, mut fc) = (a, fa);
    let mut d = b - a;
    let mut e = d;
    let mut iterations = 0usize;
    loop {
        if fb.signum() == fc.signum() {
            c = a;
            fc = fa;
            d = b - a;
            e = d;
        }
        if fc.abs() < fb.abs() {
            a = b;
            b = c;
            c = a;
            fa = fb;
            fb = fc;
            fc = fa;
        }
        let tol1 = 2.0 * f64::EPSILON * b.abs() + 0.5 * x_tolerance;
        let xm = 0.5 * (c - b);
        let midpoint = b + xm;
        if xm.abs() <= tol1 || fb == 0.0 || midpoint == b || midpoint == c {
            break;
        }
        if iterations >= MAX_ITERATIONS {
            return Ok(finish(
                b,
                fb,
                c,
                iterations,
                evaluations,
                Status::BudgetExhausted,
            ));
        }
        iterations += 1;

        if e.abs() >= tol1 && fa.abs() > fb.abs() {
            let s = fb / fa;
            let (mut p, mut q);
            if a == c {
                p = 2.0 * xm * s;
                q = 1.0 - s;
            } else {
                let qq = fa / fc;
                let r = fb / fc;
                p = s * (2.0 * xm * qq * (qq - r) - (b - a) * (r - 1.0));
                q = (qq - 1.0) * (r - 1.0) * (s - 1.0);
            }
            if p > 0.0 {
                q = -q;
            } else {
                p = -p;
            }
            if 2.0 * p < (3.0 * xm * q - (tol1 * q).abs()).min((e * q).abs()) {
                e = d;
                d = p / q;
            } else {
                d = xm;
                e = d;
            }
        } else {
            d = xm;
            e = d;
        }
        a = b;
        fa = fb;
        b += if d.abs() > tol1 { d } else { tol1.copysign(xm) };
        fb = eval(b)?;
    }
    if fb == 0.0 {
        return Ok(exact(b, iterations, evaluations));
    }
    Ok(finish(b, fb, c, iterations, evaluations, Status::Converged))
}

fn exact(x: f64, iterations: usize, evaluations: usize) -> Root {
    Root {
        x,
        value: 0.0,
        lower: x,
        upper: x,
        iterations,
        evaluations,
        status: Status::Converged,
    }
}

fn finish(b: f64, fb: f64, c: f64, iterations: usize, evaluations: usize, status: Status) -> Root {
    Root {
        x: b,
        value: fb,
        lower: b.min(c),
        upper: b.max(c),
        iterations,
        evaluations,
        status,
    }
}
