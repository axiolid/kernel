# Changelog

All notable changes to this crate are documented in this file.

The format follows [Keep a Changelog](https://keepachangelog.com/en/1.1.0/),
and this crate adheres to [Semantic Versioning](https://semver.org/).
Pre-1.0: the minor version is the breaking-change slot, per Cargo's own
caret rule for `0.x` versions.

## [Unreleased]

## [0.3.3] - 2026-10-03

### Added

- `bound` (#232): certified derivative and chord bounds.
  `chord_bound2`/`chord_bound3` bound how far a curve strays from the chord
  of a parameter span (exact sagitta for an arc of at most half a turn,
  `h^2/8 sup|c''|` otherwise, and for a rational B-spline the projective
  form `h^2/8 (|A''| + R|w''|) / w_min`, which carries no first-derivative
  terms); `curve_derivative_bounds2`/`3` bound `|c'|, |c''|, |c'''|` of
  lines, circles, ellipses, sinusoids and B-splines (rational too) from
  closed forms and derivative control polygons; `SurfaceBoundOracle`
  bounds the first and second partials of every elementary surface and of
  B-spline surfaces over a parameter box, with interpolation coefficients
  for linear interpolation over a triangle; `continuity_breaks2`/`3` name
  the knots where a curve may fail to be `C^k`; `certifies_flattening2`/`3`
  name the families whose flattening is certified.

### Changed

- `flatten2` and `flatten3` accept a span only when its certified chord
  bound is within the tolerance too, for every family `bound` covers
  (#232): the midpoint sagitta alone let an ellipse or a spline bulge past
  its chord either side of the midpoint. A B-spline is cut at its corner
  knots (multiplicity at least the degree) first, so a corner is kept as
  an exact vertex. Families without a bound keep the sagitta test.

## [0.3.2] - 2026-09-28

### Changed

- The crates.io page is this crate's own `README.md`, with links to its
  API documentation, its reference page and the source (ADR 0078).

## [0.3.1] - 2026-09-27

### Added

- Evaluation, derivatives and inversion of `Curve3::PairSection`. A
  `Curve2::Lifted` reading a pair section on one of its own B-spline
  surfaces takes that surface's parameters straight from the solve.
- `surface::locate`, `curve::locate2` and `curve::locate3`: parameters of
  a point, iterating where no closed form exists (B-spline surfaces and
  curves: seeded Newton, verified by the round trip). `invert`, `invert2`
  and `invert3` keep their closed-form-only contract.
- Evaluation, derivatives and inversion of `Curve2::Lifted`. An
  `ImplicitSection` on a B-spline carrier is inverted through the
  surface's `locate`.

- Evaluation, derivatives and inversion of `Curve2::Implicit` and
  `Curve3::ImplicitSection` (ADR 0077). `invert2` and `invert3` now also
  cover `QuadraticGraph`, `AngleGraph`, `RuledSection` and `TorusSection`,
  reading the angle off the point and trying whole turns.

- Evaluation, first and second derivatives of `Curve2::QuadraticGraph`,
  `Curve3::RuledSection`, `Curve2::AngleGraph` and `Curve3::TorusSection`
  (#119, ADR 0076); a parameter outside the graph's
  spans is refused, not extrapolated.

- `Curve2::Sinusoid` evaluation: point, first and second derivative, a
  one-turn domain, and exact inversion (the parameter is the point's first
  coordinate, then its height is checked) (ADR 0071).
