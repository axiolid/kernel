# Changelog

All notable changes to this crate are documented in this file.

The format follows [Keep a Changelog](https://keepachangelog.com/en/1.1.0/),
and this crate adheres to [Semantic Versioning](https://semver.org/).
Pre-1.0: the minor version is the breaking-change slot, per Cargo's own
caret rule for `0.x` versions.

## [Unreleased]

## [0.3.3] - 2026-09-28

### Changed

- The crates.io page is this crate's own `README.md`, with links to its
  API documentation, its reference page and the source (ADR 0078).

## [0.3.2] - 2026-09-27

### Fixed

- A `SolidOperation::BoundedHalfSpace` boundary must be a 2D curve (#162).
  The graph accepted a 3D curve, which the compiler refuses, so such a
  graph validated and then could never compile; it is now refused when the
  graph is built.

## [0.3.1] - 2026-09-27

### Added

- A `Curve2::QuadraticGraph` or `Curve2::AngleGraph` with finite
  coefficients is accepted as a trim basis (#119).

- A `Curve2::Sinusoid` is a valid trim basis when its three coefficients
  are finite (ADR 0071).
