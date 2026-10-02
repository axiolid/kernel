#![forbid(unsafe_code)]

//! Shared numeric substrate: root finding, quadrature, dense linear algebra
//! and least squares, and minimisation over `f64`.
//!
//! These are the general numeric routines geometry algorithms build on --
//! curve and surface fitting, fairing, projection and inversion -- in one
//! place instead of re-derived inline. The reference designs are OCCT's
//! `math` package, QUADPACK, Brent (1973) and MINPACK.
//!
//! | Need | Routine |
//! | --- | --- |
//! | a root of `f` between two points of opposite sign | [`find_root`] |
//! | every real root of a polynomial in an interval | [`Polynomial::real_roots`] |
//! | `integral f` with an error estimate | [`integrate`] |
//! | fixed Gauss rule for an integrand of known degree | [`GaussLegendre`] |
//! | square system | [`Lu`] |
//! | symmetric positive definite system | [`Cholesky`] |
//! | large sparse symmetric positive definite system | [`conjugate_gradient`] on a [`SparseMatrix`] |
//! | linear least squares, rank | [`Qr`], [`least_squares`] |
//! | least squares with equality constraints | [`constrained_least_squares`] |
//! | nonlinear least squares, nonlinear equations | [`levenberg_marquardt`] |
//! | minimum of `f` on an interval | [`minimize_scalar`] |
//!
//! # Contract
//!
//! - Non-finite input is refused by name ([`NumericError::NonFiniteInput`]),
//!   as is a non-finite value from a caller's function.
//! - Every answer says how far it can be trusted: a bracket that contains
//!   the root, an error estimate, a condition estimate, or an explicitly
//!   unresolved region. Which of these is a proven bound and which is an
//!   estimate is stated on each field.
//! - Iterative routines report a [`Status`]. One that ran out of budget
//!   still returns its best value and estimate; its `converged()` accessor
//!   turns that into [`NumericError::NotConverged`] for callers that only
//!   accept a met tolerance.
//! - Singular, rank-deficient, non-symmetric and indefinite matrices are
//!   refused, never solved into infinities.
//!
//! There are no geometry types here and no dependencies: values are plain
//! `f64` slices, a small row-major [`Matrix`] and a compressed-row
//! [`SparseMatrix`].

pub mod error;
pub mod linalg;
pub mod optimize;
pub mod poly;
pub mod quad;
pub mod root;
pub mod sparse;

pub use error::{NumericError, NumericResult, Status};
pub use linalg::{
    constrained_least_squares, least_squares, Cholesky, ConstrainedSolution, LeastSquaresSolution,
    LinearSolution, Lu, Matrix, Qr,
};
pub use optimize::{
    levenberg_marquardt, minimize_scalar, LeastSquaresMinimum, LevenbergMarquardtOptions,
    Residuals, ScalarMinimum, Termination,
};
pub use poly::{Polynomial, PolynomialRoot, RootKind};
pub use quad::{integrate, GaussLegendre, Integral, IntegrationOptions};
pub use root::{find_root, Root};
pub use sparse::{
    conjugate_gradient, ConjugateGradientOptions, ConjugateGradientSolution, SparseMatrix,
};
