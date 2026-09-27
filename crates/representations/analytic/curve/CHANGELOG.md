# Changelog

All notable changes to this crate are documented in this file.

The format follows [Keep a Changelog](https://keepachangelog.com/en/1.1.0/),
and this crate adheres to [Semantic Versioning](https://semver.org/).
Pre-1.0: the minor version is the breaking-change slot, per Cargo's own
caret rule for `0.x` versions.

## [Unreleased]

## [0.3.1] - 2026-09-27

### Added

- `Field2::value` and `SeriesField2::value` evaluate the value alone,
  without the jet, and agree with `jet` to the last bit.
- Bridge cells in `ImplicitCurve2` (`ImplicitCell::bridge`, a cubic into a
  point where two branches cross, bounded by its Bezier control values),
  `ImplicitCell::part`, `reversed` and `solved_range`,
  `ImplicitCurve2::solve_cell`, and `Field2::scale_at` (the size of the
  terms that make up the value at a point, which its rounding scales
  with).

- `Curve3::PairSection` (`PairSection3`, `PairNode`): the section of two
  B-spline surfaces, carried by nodes on both and defined between them by
  the surfaces themselves (where both meet on the plane across the chord).
  It has `solve` (the parameters on both surfaces and the point), `rates`,
  `second_rates`, `sub`, `reversed`, `parameter_of` and `side`, and
  `pair_section::solve4` for the 4x4 systems behind it (ADR 0077).
- `BSplineSurface` is defined here, and `axiolid_surface` re-exports it
  unchanged, so a traced curve can carry a B-spline carrier
  (`Carrier::Spline`). It has `BSplineSurface::jet` (point and first and
  second partials, rational) and `domain`.
- `Field2` is now an enum of `SeriesField2` (the former struct: powers and
  harmonics) and `PatchField2` (piecewise Bernstein polynomials on a grid,
  bounded by their coefficients over any box and continued past the grid by
  their edge polynomials). `ImplicitCurve2::clipped` cuts a curve to a box.
- `Curve2::Lifted(LiftedCurve2)`: a space curve read in an analytic
  surface's parameters, sharing the curve's parameter. It is the pcurve,
  on the analytic face, of a section only a B-spline can carry.

- `ImplicitCurve2::sub`, `rotated`, `reversed`, `shifted`, `closure` and
  `turning_points`. `implicit::{bound, bound_simple, partial}` give
  interval bounds and partial derivatives of a `Field2` over parameter
  boxes.

- `Curve2::Implicit(ImplicitCurve2)` and
  `Curve3::ImplicitSection(ImplicitSection3)` (#119, ADR 0077): a stretch of
  a `Field2`'s zero set in monotone cells, where each point is the field's
  unique root in its cell's bracket, and the same curve on its analytic
  `Carrier` (plane, ruled surface, sphere or torus).

- `Curve2::QuadraticGraph(QuadraticGraph2)` and
  `Curve3::RuledSection(RuledSection3)` (#119, ADR 0076): one root branch of
  `a(t) v^2 + b(t) v + c(t) = 0` with degree-2 trigonometric coefficients
  (`Trig2`, `Branch`), and the same curve lifted onto a cylinder,
  elliptical cylinder or cone (`RuledCarrier`). The exact pcurve and edge
  of a quadric's cut across a ruled surface.
- `Curve2::AngleGraph(AngleGraph2)` and `Curve3::TorusSection(TorusSection3)`
  (#119, ADR 0076): the solution `u(t)` of `a(t) cos u + b(t) sin u = c(t)`,
  and the same curve on a torus (`TorusCarrier`). The exact pcurve and edge
  of a plane's or sphere's cut across a torus.

- `Curve2::Sinusoid(Sinusoid2)`: the graph `v = mean + a cos(t) + b sin(t)`,
  the exact pcurve of a plane's cut across a cylinder in its (angle, height)
  parameters (ADR 0071). The parameter is the first coordinate. Additive:
  `Curve2` is `#[non_exhaustive]`.
