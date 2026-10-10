# Changelog

All notable changes to this crate are documented in this file.

The format follows [Keep a Changelog](https://keepachangelog.com/en/1.1.0/),
and this crate adheres to [Semantic Versioning](https://semver.org/).
Pre-1.0: the minor version is the breaking-change slot, per Cargo's own
caret rule for `0.x` versions.

## [Unreleased]

### Added

- Curve paths (#290, ADR 0082 amendment): `path::CurvePath`, pieces of
  curves laid end to end, each a `PathPiece` -- the span `[start, end]` of
  a `PathCurve` (an owned `Curve2` or `Curve3`) in its station measure,
  traversed forwards or backwards, optionally carried by a rigid placement
  that is or is not exact. The neutral form of a composite, trimmed or
  segmented curve relation, which the curve-evaluation contract can see:
  it holds and composes (`reversed`, `placed`, `extend`, `length`,
  `PathPiece::frame_is_exact`) and measures nothing. `PathCurve` and
  `PathPiece` are `#[non_exhaustive]`, so a later piece kind is additive.
- Offset pieces of a curve path (#289, ADR 0082 amendment):
  `PathCurve::Offset(Box<PathOffset>)`, appended last; `PathOffset { base,
  law }`, an offset of one base piece measured in its own length;
  `OffsetLaw::Planar` (an offset curve 2D, left positive),
  `OffsetLaw::Directed` (an offset curve 3D, along `normalise(V x T)`) and
  `OffsetLaw::Linear` (offsets at stations, linear between the piece's
  ends, in an `OffsetFrame`), with `PathOffsets` and `OffsetFrame`. All
  `#[non_exhaustive]`. `PathCurve::is_line` is also true for an offset of a
  line by a constant law, and `PathPiece::frame_is_exact` also needs the
  offset's base placed exactly.

## [0.3.5] - 2026-10-09

### Added

- Rotation about a held rail (#279, ADR 0081 amendment):
  `CantForm::AboutRail { rail, elevation }`, a pivot piece whose
  elevation is derived from the cant law at the same plan distance,
  `e = e0 + s D / 2 = e0 + s (b / 2) sin(psi)`, `s = +1` about the right
  rail and `-1` about the left (`RailSide`, `RailSide::pivot_sign`). A
  Viennese bend, a polynomial in the bank angle, now rotates about its
  low rail exactly; under `BankConvention::VerticalRise` the held rail
  head stands `e0` above the profile at every station. Built with
  `CantPiece::about_rail`; `CantLaw::has_rail_pieces` finds one.
- `Banked3::cant_rate_at`: `dD/dd`, the law's own rate or
  `b cos(psi) psi'` for an angle piece.
- `BankError::RailInCant`: a held-rail piece in a cant law is refused by
  name.

### Changed

- `Banked3::pivot_at` reads a held-rail piece from the cant law: value
  `e0 + s D / 2` and rate `s D' / 2`, with the cant law's refusals at
  that station. Height pieces read as before; an angle piece in the
  pivot law is still `BankError::AngleInPivot`. `CantLaw::value_at` and
  `rate_at` have no value on a held-rail piece.

## [0.3.4] - 2026-10-09

### Added

- `SeamSide` (#263, ADR 0082 amendment): which piece a position exactly on
  a seam of a composite curve is read from -- `Outgoing`, the piece that
  starts there and what every evaluator reads (the default), or
  `Incoming`, the piece that ends there. `#[non_exhaustive]`; shared by
  `axiolid-model`'s stations and `axiolid-evaluate`'s station readers.

## [0.3.3] - 2026-10-03

### Added

- `Curve3::Banked(Banked3)` (#240, ADR 0081): an elevated centreline
  carrying a roll law, evaluated by plan distance. `CantLaw` is a run of
  `CantPiece`s from plan distance zero, each a `CantForm` over its own
  `xi = s / length`: a polynomial in `xi` (constant, linear, Bloss, and
  Helmert as two quadratic pieces through `CantPiece::helmert`), a
  half-cosine, a sine transition, or a Viennese bend, which gives the bank
  angle itself (`CantValue::Angle`). A second `CantLaw` of height pieces
  is the pivot: the rotation point's elevation above the profile. The
  rail-head distance `b` gives `psi = asin(D / b)` (`bank_angle`), and a
  cant beyond it is refused by name (`BankError`). `BankConvention` names
  how a cant is read on a grade, with no default: `TangentRotation` rolls
  the section about the 3D tangent by `psi` (the rail heads rise
  `D cos theta`), `VerticalRise` by `asin(D / (b cos theta))` (they rise
  exactly `D`).
- `ElevationLaw::CircularArc { height, grade, radius }` (#238): a vertical
  circular arc in plan distance, the circle itself rather than a parabola.
  With `t0 = atan(grade)` and signed `R` (positive sag), `sin t = sin t0 +
  d/R`, `z = height + R (cos t0 - cos t)` and `grade = tan t`, evaluated
  in a form that does not cancel for small `d/R`. `height_at` and
  `grade_at` answer it in closed form and report `None` where
  `|sin t0 + d/R| >= 1`.
- `ElevationLaw::Intrinsic { height, grade, curvature }` (#238): a profile
  given by curvature against its own arc length (a linear law is the
  clothoid between grades). Stored exactly; `height_at` and `grade_at`
  report `None` for it, since its height is a quadrature, read by
  `axiolid-evaluate`.
- `ElevationLaw::circular_arc`, `ElevationLaw::intrinsic` and
  `ElevationLaw::piece_at`, the innermost piece of a composed law at a
  distance, rebased to its own start.
- `Curve2::Chain(Chain2)` (#239): a plane curve parameterised by
  cumulative arc length, its pieces placed rigidly end to end from a start
  frame. `ChainPiece2::Intrinsic { curvature, length }` is a curvature-law
  piece; `ChainPiece2::Parametric { curve, start, length }` reads any
  `Curve2`, written in the piece's local frame (through the origin with
  tangent `+x` at `start`), by arc length over `length`. Each piece starts
  at the previous end point along the previous end tangent, so a chain is
  tangent-continuous at its joins: a layout such as line, cubic parabola,
  arc is one exact value. `Chain2::is_well_formed`, `length`, `joins` and
  `piece_at` read the data; evaluation lives in `axiolid-evaluate`.
  `Elevated3` takes a chain as its plan, since its parameter is plan
  distance (ADR 0060).

## [0.3.2] - 2026-09-28

### Changed

- The crates.io page is this crate's own `README.md`, with links to its
  API documentation, its reference page and the source (ADR 0078).

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
