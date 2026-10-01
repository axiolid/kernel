//! Refusals and convergence status shared by every routine.

use core::fmt;

/// Why a numeric routine refused its input or could not produce a value.
///
/// Every variant names the offending argument or condition. A routine never
/// substitutes a default, clamps silently, or returns NaN in place of one of
/// these.
#[derive(Debug, Clone, PartialEq)]
#[non_exhaustive]
pub enum NumericError {
    /// An input value was NaN or infinite. `name` is the argument.
    NonFiniteInput {
        /// The argument that was not finite, such as `"lower bound"`.
        name: &'static str,
    },
    /// A caller-supplied function returned NaN or infinity.
    NonFiniteEvaluation {
        /// Which function, such as `"integrand"` or `"residuals"`.
        name: &'static str,
        /// The scalar argument it was evaluated at, when it has one.
        at: Option<f64>,
    },
    /// An argument was finite but outside its valid range.
    InvalidArgument {
        /// The argument, such as `"x tolerance"`.
        name: &'static str,
        /// What it must satisfy.
        reason: &'static str,
    },
    /// Sizes of two arguments disagree.
    DimensionMismatch {
        /// The argument whose size is wrong.
        name: &'static str,
        /// The size it must have.
        expected: usize,
        /// The size it has.
        found: usize,
    },
    /// The function has the same sign at both ends of the bracket, so a
    /// root is not guaranteed between them.
    NoSignChange {
        /// Value at the lower end.
        f_lower: f64,
        /// Value at the upper end.
        f_upper: f64,
    },
    /// Every coefficient of the polynomial is zero: every point is a root.
    ZeroPolynomial,
    /// A square matrix is singular to working precision.
    Singular {
        /// Which matrix, such as `"system matrix"`.
        name: &'static str,
    },
    /// A least-squares matrix has fewer independent columns than unknowns,
    /// so the solution is not unique.
    RankDeficient {
        /// Which matrix.
        name: &'static str,
        /// Numerical rank found.
        rank: usize,
        /// Rank required (the number of unknowns or constraints).
        required: usize,
    },
    /// A matrix given to a symmetric solver is not symmetric.
    NotSymmetric {
        /// Row of the first asymmetric pair found.
        row: usize,
        /// Column of the first asymmetric pair found.
        column: usize,
    },
    /// A symmetric matrix is not positive definite to working precision.
    NotPositiveDefinite {
        /// Index of the first non-positive pivot.
        pivot: usize,
    },
    /// An iterative routine stopped without meeting its tolerance. Returned
    /// only by the `converged()` accessors; the routines themselves report
    /// this case in their result's status with the estimate they reached.
    NotConverged {
        /// Which routine.
        name: &'static str,
        /// The error estimate it reached.
        error_estimate: f64,
    },
}

impl fmt::Display for NumericError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::NonFiniteInput { name } => write!(f, "{name} is not finite"),
            Self::NonFiniteEvaluation { name, at: Some(x) } => {
                write!(f, "{name} returned a non-finite value at {x}")
            }
            Self::NonFiniteEvaluation { name, at: None } => {
                write!(f, "{name} returned a non-finite value")
            }
            Self::InvalidArgument { name, reason } => write!(f, "{name} {reason}"),
            Self::DimensionMismatch {
                name,
                expected,
                found,
            } => write!(f, "{name} has size {found}, expected {expected}"),
            Self::NoSignChange { f_lower, f_upper } => write!(
                f,
                "function values {f_lower} and {f_upper} at the bracket ends do not change sign"
            ),
            Self::ZeroPolynomial => f.write_str("polynomial is identically zero"),
            Self::Singular { name } => write!(f, "{name} is singular to working precision"),
            Self::RankDeficient {
                name,
                rank,
                required,
            } => write!(f, "{name} has numerical rank {rank}, {required} required"),
            Self::NotSymmetric { row, column } => {
                write!(f, "matrix is not symmetric at ({row}, {column})")
            }
            Self::NotPositiveDefinite { pivot } => {
                write!(f, "matrix is not positive definite (pivot {pivot})")
            }
            Self::NotConverged {
                name,
                error_estimate,
            } => write!(
                f,
                "{name} did not converge (error estimate {error_estimate})"
            ),
        }
    }
}

impl std::error::Error for NumericError {}

/// Result of a numeric routine.
pub type NumericResult<T> = Result<T, NumericError>;

/// How an iterative routine stopped.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Status {
    /// The requested tolerance was met.
    Converged,
    /// The iteration, evaluation or subdivision budget ran out first. The
    /// result still carries the best value and its error estimate.
    BudgetExhausted,
    /// Rounding error stopped further progress before the tolerance was
    /// met; the error estimate is the best `f64` can certify here.
    RoundoffLimited,
}

/// Refuse a non-finite scalar by name.
pub(crate) fn finite(value: f64, name: &'static str) -> NumericResult<f64> {
    if value.is_finite() {
        Ok(value)
    } else {
        Err(NumericError::NonFiniteInput { name })
    }
}

/// Refuse a slice containing a non-finite value by name.
pub(crate) fn all_finite(values: &[f64], name: &'static str) -> NumericResult<()> {
    if values.iter().all(|v| v.is_finite()) {
        Ok(())
    } else {
        Err(NumericError::NonFiniteInput { name })
    }
}

/// Refuse a tolerance that is negative or not finite.
pub(crate) fn tolerance(value: f64, name: &'static str) -> NumericResult<f64> {
    finite(value, name)?;
    if value < 0.0 {
        return Err(NumericError::InvalidArgument {
            name,
            reason: "must be non-negative",
        });
    }
    Ok(value)
}
