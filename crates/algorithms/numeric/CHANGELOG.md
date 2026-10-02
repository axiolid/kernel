# Changelog

All notable changes to this crate are documented in this file.

The format follows [Keep a Changelog](https://keepachangelog.com/en/1.1.0/),
and this crate adheres to [Semantic Versioning](https://semver.org/).
Pre-1.0: the minor version is the breaking-change slot, per Cargo's own
caret rule for `0.x` versions.

## [Unreleased]

## [0.1.1] - 2026-10-02

### Added

- `SparseMatrix`: compressed sparse rows assembled from triplets, with
  duplicates summed in input order so storage and products do not depend
  on triplet order.
- `conjugate_gradient`: Jacobi-preconditioned conjugate gradients for a
  sparse symmetric positive definite system (#129). Convergence is judged
  on the recomputed true residual, a spent budget is reported as
  `Status::BudgetExhausted` and refused by `converged()`, and an
  asymmetric matrix, a non-positive diagonal or a direction of
  non-positive curvature is refused by name.

## [0.1.0] - 2026-10-02

### Added

- New crate (#136): a shared numeric substrate with no dependencies.
- `find_root`: Brent's method on a sign-change bracket, with a forced
  bisection whenever three steps fail to halve the bracket; returns the
  final bracket as a proven error bound.
- `Polynomial::real_roots`: every real root in an interval, isolated by
  derivative roots, with rigorous rounding bounds on each evaluation.
  Multiple roots and unseparable clusters come back as `Unresolved`
  regions with the most roots they can hold.
- `integrate`: adaptive Gauss-Kronrod (G7/K15) with QUADPACK's error
  estimate and a panel budget; `GaussLegendre` fixed rules up to 512
  points.
- `Matrix`, `Lu`, `Cholesky` and column-pivoted `Qr`, with Hager-Higham
  condition estimates and forward error estimates; `least_squares` and
  `constrained_least_squares` (null-space method).
- `levenberg_marquardt` over a `Residuals` trait (analytic or
  forward-difference Jacobian), and `minimize_scalar` (Brent's localmin).
  Both report a `Status`; `converged()` turns a missed tolerance into an
  error.
- Non-finite input, singular, rank-deficient, asymmetric and indefinite
  matrices, and unbracketed roots are refused by name.
