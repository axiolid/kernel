# Changelog

All notable changes to this crate are documented in this file.

The format follows [Keep a Changelog](https://keepachangelog.com/en/1.1.0/),
and this crate adheres to [Semantic Versioning](https://semver.org/).
Pre-1.0: the minor version is the breaking-change slot, per Cargo's own
caret rule for `0.x` versions.

## [Unreleased]

### Added

- `surface::invert` (and so `locate`) and `surface::jet` take an
  elliptical cylinder (#287): the angle is that of the section
  coordinates over the semi-axes, the inverse of `evaluate`'s affine
  image of the circle. They were refused as unsupported. `project` still
  refuses it: the closest point on an ellipse has no closed form.

- `ReferenceCurveEvaluator` reads curve paths (#290, ADR 0082 amendment):
  the contract's `path_*` queries build `CompositeBasis::from_path` and
  answer with `CompositeBasis::section_on` for the same side, so a frame
  through the contract on a curve relation equals the same station lowered
  as geometry, bitwise against `+Z`; against another reference up the
  reference-up frame of the section's point and tangent (a banked piece
  refused there by name). A native parameter along a path is refused by
  name; the convention and exactness are the composite's.
- `CompositeBasis::from_path` (the path's curves borrowed, each piece
  checked against its curve, then measured as `CompositeBasis::new` does)
  and `CompositeBasis::path`; `StationPiece::from_path_piece`,
  `StationCurve::from_path_curve`, and `From<StationPiece>` /
  `From<StationCurve>` for `PathPiece` / `PathCurve` (the curve copied).
- Offsets as station bases (#289, ADR 0082 amendment): `station::offset`,
  with `StationOffset` (an offset of one base piece by an `OffsetLaw`,
  measured in its own length: arc length, or its own plan length beside an
  elevated or banked curve), `StationCurve::Offset` (appended last),
  `StationPiece::offset` and `offset_pieces` (one offset piece per span of
  the base between seams, so every seam and joint of the base, and every
  break of a by-distances law, is a joint of the offset). A constant offset
  of a line is a line (frame exact), of a circle a circle of radius
  `r - a` (length in closed form); every other offset is read numerically
  to `OFFSET_TOLERANCE` and never claimed exact. A collapse, a cusp or
  reversal (`CUSP_TOLERANCE`), a self-crossing in a horizontal plane, an
  offset of an offset and a base span across a seam are refused by name,
  as is an offset across a corner of its base, where its two sides do not
  meet. `CompositeBasis` (and so `from_path` and the reference provider's
  `path_*` queries) reads offset pieces, building each one's measure table
  once.

## [0.3.9] - 2026-10-09

### Added

- `ReferenceCurveEvaluator` reads a seam side (#286, ADR 0082
  amendment): `point_at_on`, `tangent_at_on` and `frame_at_on` answer a
  measure ON a seam by the station rule (within `ARC_LENGTH_TOLERANCE *
  max(1, s)`, for both sides) with `station_section3_on`, so a placement
  framed through the contract agrees with the same station lowered as
  geometry; the frame is built against the evaluator's reference up, and
  a banked curve's against `+Z` only, as side-lessly. Off a seam both
  sides are the side-less answer. A distance is located in the provider's
  convention, as is an elevated or banked curve's native parameter; the
  incoming side at a polyline's or a B-spline's native parameter is
  refused by a typed `UnsupportedInput`.
- Banked curves rotating about a held rail (#279, ADR 0081 amendment):
  a pivot piece `CantForm::AboutRail` is evaluated from the cant law,
  point, tangent and section alike. `banked_second_derivative` adds
  `e'' = s D'' / 2`, for an angle piece
  `s (b / 2) (cos(psi) psi'' - sin(psi) psi'^2)`.
  `banked_derivative_bounds` and `banked_chord_bound` take half the
  cant's bounds there, for a Viennese bend `(b / 2) P_1`,
  `(b / 2)(P_2 + P_1^2)` and `(b / 2)(P_3 + 3 P_1 P_2 + P_1^3)` with
  `P_k` the bend's exact derivative suprema over the span (its ends and
  critical points), so `flatten3` and the sweep certification carry
  over. `banked_breaks` and `grade_corners3` name the cant's seams inside
  a held-rail piece as the point path's. A cant law with a held-rail
  piece is refused by name.

## [0.3.8] - 2026-10-09

### Added

- Composite station bases (#285, ADR 0082 amendment), in
  `station::composite`: `CompositeBasis` measures stations along
  `StationPiece`s -- spans of atomic 2D or 3D curves (`StationCurve`)
  between two distances in the curve's station measure, reversed or
  carried by a rigid placement -- laid end to end. The distance runs
  through the pieces in their common convention (plan distance when
  every piece is an elevated or banked curve, arc length when none is;
  mixed pieces refused by name); consecutive pieces must meet within
  `JOINT_TOLERANCE` (relative), a gap and an undeclared reversed piece
  refused by name; every joint is a seam read by the #263 rule
  (`section_on` with a `SeamSide`); `seams` / `exact_seams` list the
  joints and each piece's own seams; `frame_is_exact_at` claims an exact
  frame only where every piece up to the one read is an exactly placed
  line; `pieces_between` clips the pieces a trim of the composite keeps.
- `SectionFrame::carried`: the frame a station on a curve placed by a
  rigid motion has -- the source's moved where the motion keeps `+Z`
  (within the new `KEEPS_UP_TOLERANCE`), the placed curve's own
  reference-up frame of the moved point and tangent where it tilts `+Z`
  on an arc-length-measured source; an elevated or banked source so
  placed is refused by a typed `UnsupportedInput`.
- `station::DistanceConvention`, the curve-evaluation contract's
  convention re-exported so a caller can name a basis's measure.

## [0.3.7] - 2026-10-09

### Added

- Placing at a station (#264, ADR 0082 amendment):
  `SectionFrame::placement` is the rigid motion taking local `x`, `y`,
  `z` onto the tangent, the left lateral and up and the origin onto the
  point (a linear placement's reading, not the provider layout of
  `SectionFrame::frame`); `SectionFrame::moved` carries a frame by a
  rigid motion; `station_frame_is_exact2` / `station_frame_is_exact3`
  say whether a station frame is exact (a line only: every other family
  reads its distance by the arc-length inverse or its point by
  quadrature).

- Seams of a station's basis curve (#263, ADR 0082 amendment), in
  `station::seam`: `station_seams2` / `station_seams3` list where two
  pieces of a curve meet -- a polyline's vertices at the running sum of
  its segment lengths, a B-spline's corner knots (multiplicity at least
  its degree; by quadrature, flagged not exact), an intrinsic curve's and
  a chain's law seams and joins (smooth), an elevated curve's plan and
  profile seams, a banked curve's cant and pivot seams -- as
  `StationSeam { distance, parameter, smooth, exact }`, without
  evaluating the curve; `exact_station_seams2` / `exact_station_seams3`
  refuse an inexact seam with a typed `UnsupportedInput`.
- `station_section2_on` / `station_section3_on`: a station within the
  arc-length tolerance of a seam that is not smooth is read at the seam
  from the piece the `SeamSide` names, the incoming one as the curve
  truncated at the seam (a polyline's previous segment, a B-spline's
  previous span, a profile and cant and pivot law cut there).
- `Mitre`: the plane a run of sections crossing a seam stands in, normal
  to the bisector of the two tangents; `Mitre::between` is `None` within
  `SEAM_TANGENT_TOLERANCE` and refuses a near reversal (`MITRE_TOLERANCE`
  on the cosine of half the turn) and a seam whose sides do not share
  their point by name; `Mitre::place` projects a section point placed in
  each side's frame along that side's tangent onto the plane.

### Changed

- `station_section2` / `station_section3` read a station on a seam (within
  `ARC_LENGTH_TOLERANCE * max(1, s)`) at the seam itself, from the piece
  starting there. On a polyline or a B-spline the arc-length inverse could
  previously land a hair either side of the vertex, and so either frame.

## [0.3.6] - 2026-10-04

### Added

- Elevated curves through the generic 3D curve functions (#252):
  `curve::evaluate3`, `derivative3`, `second_derivative3` and `domain3`
  read a `Curve3::Elevated` by PLAN distance (ADR 0082, like `Banked3`):
  the derivative is `(p', grade)`, not unit, the second `(k n, z'')`,
  the domain `[0, L]` over the plan's length (unbounded for a line plan,
  a full turn for a circle); `second_derivative3` also reads a
  `Curve3::Banked` (`banked::banked_second_derivative`, the pivot's `e''`
  added).
- `elevated` module (#252): `elevated_derivative`,
  `elevated_second_derivative`, `elevated_span`, and the certified
  `elevated_chord_bound` -- `sqrt(P^2 + Z^2)`, the plan's Taylor chord
  bound `h^2/8 sup |k|` and the profile's `elevation_chord_bound` at
  equal plan distance, orthogonal parts composed in quadrature (the
  module documents the derivation and why the plan's own `chord_bound2`
  does not compose) -- with `banked_chord_bound`,
  `elevated_derivative_bounds` and `banked_derivative_bounds`
  (`sup |c'|, |c''|, |c'''|`, plan and profile combined per law, a
  chain's parametric piece through its own derivative suprema),
  `elevated_breaks`/`banked_breaks`, and `grade_corners3`, the seams
  across which the grade jumps, read from the laws on either side.
  `bound::chord_bound3`, `curve_derivative_bounds3`, `continuity_breaks3`
  and `certifies_flattening3` dispatch to them, so `flatten3` certifies an
  elevated curve over a certified plan with a closed-form profile.

### Changed

- `elevation`: an intrinsic reading carries its profile arc length, and
  `frenet`'s pointwise curvature-law value is shared crate-wide (#252);
  no public signature changed.

## [0.3.5] - 2026-10-03

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
