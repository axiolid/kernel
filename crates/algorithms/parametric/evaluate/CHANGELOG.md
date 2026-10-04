# Changelog

All notable changes to this crate are documented in this file.

The format follows [Keep a Changelog](https://keepachangelog.com/en/1.1.0/),
and this crate adheres to [Semantic Versioning](https://semver.org/).
Pre-1.0: the minor version is the breaking-change slot, per Cargo's own
caret rule for `0.x` versions.

## [Unreleased]

### Added

- `station::SectionFrame::oriented` (#246, ADR 0082 amendment): turns a
  section frame by an explicit axis and reference direction given as
  components in that frame, `(tangent, lateral, up)` -- the axis becomes
  the exact up, the reference direction is orthonormalised against it
  (Gram-Schmidt, axis primary), the point stays -- and refuses a zero,
  non-finite or parallel pair by name (`station::ORIENTATION_TOLERANCE`).
- `bound::chord_bound2` bounds implicit curves (`Curve2::Implicit`, ADR
  0077) one cell at a time (#249), and `continuity_breaks2` names their
  cell joins. A regular cell's solved parameter is bounded through the
  implicit function theorem with interval bounds of the field's partials
  over the stretch's box (first and second order, the smaller kept); a
  bridge cell into a crossing by its cubic's Bezier control points. The
  deviation of a boolean whose exact result has traced pcurves (a hole
  through an I-beam's root fillets) can now be certified.

## [0.3.5] - 2026-10-03

### Added

- `station` (#241, ADR 0082): `station_section2` and `station_section3`
  resolve a distance along a curve -- plan distance on an elevated or
  banked curve, arc length (through `arc_parameter`) on any other -- to a
  `SectionFrame`: the point, the unit tangent, the lateral axis to the
  left and `up = tangent x lateral`. A 2D curve is framed in `z = 0` with
  `+Z` up, a banked curve by its rolled section, every other 3D curve by
  the reference-up frame against `+Z` that `ReferenceCurveEvaluator`
  returns. `SectionFrame::place` applies offsets, `plan` gives the upright
  frame, `frame` the provider's layout (`x` tangent, `y` up, `z` right).
  `station_length2`/`station_length3` give the length a station is
  checked against; a negative, non-finite or too-long distance and a
  vertical tangent are refused by name.

## [0.3.4] - 2026-10-03

### Added

- `banked` (#240, ADR 0081): `banked_point`, `banked_derivative`,
  `banked_tangent` and `banked_section` evaluate a `Curve3::Banked` by plan
  distance. A `BankedSection` carries the rotation point, the unit tangent
  of the banked curve (its grade includes the pivot's rate), the lateral
  and up axes rolled about it, the cant, the nominal bank angle, the roll
  its convention gives, the pivot and the grade, with `frame()` and
  `rail_heads()`. Every refusal of the cant laws is passed on by name.
- `evaluate3`, `derivative3` and `domain3` take a banked curve, parameterised
  by plan distance over its cant law's span, so `flatten3` flattens it;
  `certifies_flattening3` reports it uncertified (midpoint sagitta, #232).
- `ReferenceCurveEvaluator` measures a banked curve in plan distance, and
  its `frame_at` returns the curve's section frame. An evaluator built
  against a reference up other than `+Z` refuses a banked curve, whose cant
  is measured against `+Z`.
- `elevation` (#238): `elevation_height` and `elevation_grade` read every
  elevation law, alone or inside a piecewise one. Closed forms go through
  `ElevationLaw::height_at`/`grade_at`; an intrinsic (clothoid) profile is
  integrated with `intrinsic_point` and its plan distance inverted by a
  bracketed Newton solve to `INVERSION_TOLERANCE * max(1, d)`, after
  certifying that the profile stays below vertical over the bracket, so a
  profile that turns vertical is refused rather than read on a branch
  where plan distance runs backwards. Against a 50-digit reference a
  150 m clothoid agrees to 4e-15 m in height.
- `elevation_chord_bound`: a certified bound on how far an elevation law
  strays from the chord of its heights over a span, `h^2/8 sup|z''|` with
  `z'' = k / cos^3 t`, for polynomials, circular arcs and intrinsic
  profiles; `None` across a piecewise seam, outside a law's domain, or for
  a family it does not bound.

- `arc_parameter` (#239): `arc_length2`/`arc_length3` (signed arc length
  between two parameters of any evaluable curve) and
  `parameter_at_arc_length2`/`parameter_at_arc_length3` (the parameter at
  a signed arc length from a start parameter). Adaptive 8-point
  Gauss-Legendre of the speed, split at knots and vertices, and a
  bracketed Newton root find, both to `ARC_LENGTH_TOLERANCE` (`1e-12`)
  relative to `max(1, length)`; lines, intrinsic curves and chains are
  exact. A length past the end of a bounded curve (named with the length
  available), a start outside the domain, non-finite input and a
  quadrature past `MAX_PANELS` are refused by name. Against the binomial
  series of a cubic parabola the inverse lands within `1e-13` relative.
- `chain` (#239): `chain_point` and `chain_tangent` evaluate a
  `Curve2::Chain` by arc length; `evaluate2`, `derivative2` and `domain2`
  take a chain, and an `Elevated3` over a chain evaluates by plan
  distance (`elevated_point`, `elevated_tangent`, `ReferenceCurveEvaluator`,
  and so `Curve3::Banked` over it). A piece longer than its curve, a
  parametric piece off its local frame by more than
  `PIECE_FRAME_TOLERANCE`, a malformed chain and a distance outside the
  chain are refused by name.
- Chains in the certified flattening (#232): `chord_bound2` bounds a span
  inside one piece (an intrinsic piece by `h^2/8` times a bound on its
  curvature law, a parametric piece by its curve's own chord bound between
  the parameters the span reads), `continuity_breaks2` names the joins,
  and `certifies_flattening2` accepts a chain whose intrinsic pieces have
  bounded laws and whose parametric pieces are certified families with no
  corner.

### Changed

- `elevated_point`, `elevated_tangent` and `banked_derivative` read the
  elevation through `elevation_height`/`elevation_grade`, so elevated and
  banked curves carry circular-arc and intrinsic profiles.

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
