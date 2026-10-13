# Changelog

All notable changes to this crate are documented in this file.

The format follows [Keep a Changelog](https://keepachangelog.com/en/1.1.0/),
and this crate adheres to [Semantic Versioning](https://semver.org/).
Pre-1.0: the minor version is the breaking-change slot, per Cargo's own
caret rule for `0.x` versions.

## [Unreleased]

### Added

- Joints near a window (#305, ADR 0082 amendment):
  `CurveEvaluator::path_joint_distances_within(path, around, window,
  tolerance)`, a defaulted trait method answering, as `(joint, interval)`
  in order, only the interior joints whose certified interval meets
  `[around - window, around + window]`, each certified as
  `path_joint_distances` certifies it; a joint left out truly lies outside
  the window. A provider certifies only the pieces before and within the
  window: one piece whose length cannot be certified after it no longer
  refuses every placement on the path, while one before it, or within
  reach of it, still refuses by name. The default delegates to
  `path_joint_distances` and keeps the joints that meet the window, so a
  provider that does not certify joints refuses it naming
  `PATH_JOINT_DISTANCES_UNSUPPORTED`. `window_span` checks the window
  (a centre that is not finite, a window negative or not finite:
  `InvalidInput`) and `DistanceInterval::meets` tests an interval against
  it. The conformance suite asks four lines turning at right angles and a
  widening offset then a line for windows on, between, over and past
  their joints: every returned joint must be certified and meet its
  window, every joint truly within a window must be returned, each must
  agree with the provider's whole-path answer, and a bad window or
  tolerance must be refused; a provider may refuse it as unsupported.

## [0.3.6] - 2026-10-12

### Added

- Open curve-path ends (#298, ADR 0082 amendment): the trait documentation
  states that a provider reading paths resolves a piece's open end with
  its own length of the curve, and `OPEN_PATH_END_UNSUPPORTED` names the
  typed refusal of a provider that reads stated spans only. The
  conformance suite reads an open path (the whole of a widening offset of
  the whole of a line, then a line to its end) forwards and reversed,
  refuses a distance past its end and certifies its joint at `sqrt(101)`;
  a provider may refuse it as unsupported, typed.

## [0.3.5] - 2026-10-10

### Added

- Certified joint distances of a curve path (#294, ADR 0082 amendment):
  `CurveEvaluator::path_joint_distances(path, tolerance)`, one
  `DistanceInterval` per interior joint that contains the joint's true
  distance (the sum of the true lengths of the pieces before it) and is no
  wider than `tolerance`; lengths read from stored data sum to point
  intervals. A defaulted trait method, so no provider breaks: the default
  refuses by a typed `UnsupportedInput` naming the new
  `PATH_JOINT_DISTANCES_UNSUPPORTED`, never an estimate handed out as a
  certificate. `DistanceInterval` (`lower`, `upper`, `#[non_exhaustive]`,
  built by `new` and `point`) is new.
- The conformance suite checks the query: a provider either refuses it as
  unsupported, or answers one interval containing the joint of a line and
  a quarter arc (exactly at 10) and of an offset of a line whose lateral
  offset widens from 0 to 1 over 10 (`sqrt(101)`, a quadrature for a
  provider), no wider than the tolerance, and refuses a non-positive or
  non-finite tolerance and a path whose pieces do not meet.

## [0.3.4] - 2026-10-10

### Added

- Point, tangent and frame queries over a curve path (#290, ADR 0082
  amendment): `CurveEvaluator::path_point_at_on`, `path_tangent_at_on`
  and `path_frame_at_on` read a `CurvePath` (re-exported from
  `axiolid-curve` with `PathPiece` and `PathCurve`) as a composite station
  basis is read -- the distance end to end, every joint a seam read from
  the piece the `SeamSide` names -- and `path_point_at`, `path_tangent_at`
  and `path_frame_at` are their outgoing reading; `path_distance_convention`
  and `path_frame_is_exact_at` say which distance runs through the path and
  whether a frame on it is exact. All are defaulted, so no provider
  breaks: the sided queries refuse by a typed `GeomError::UnsupportedInput`
  naming the new `CURVE_PATH_UNSUPPORTED`, the plain ones delegate to them,
  the convention is `Unsupported` and nothing is claimed exact.
- Conformance: a line and a quarter arc meeting at a right angle, forwards
  and reversed, on the joint and within the seam tolerance either side --
  each side reads its own piece, the plain queries the outgoing one, a
  distance off the path and a path with a gap are refused, the convention
  is arc length and nothing on or after the arc is exact -- or every path
  query refused as unsupported.
- Re-exports `PathOffset`, `PathOffsets`, `OffsetLaw` and `OffsetFrame`
  from `axiolid-curve` (#289): a curve path may carry offset pieces, which
  a provider that does not read them refuses by name.

## [0.3.3] - 2026-10-09

### Added

- A seam side on point, tangent and frame queries (#286, ADR 0082
  amendment): `CurveEvaluator::point_at_on`, `tangent_at_on` and
  `frame_at_on` take a `SeamSide` (re-exported from `axiolid-curve`) and
  read a measure on a seam -- a polyline's vertex, a grade break, a cant
  jump -- from the piece that ends there (`Incoming`) or starts there
  (`Outgoing`), by the station seam rule. They are defaulted, so no
  provider breaks: `Outgoing` delegates to the side-less method, and
  every other side is refused by a typed `GeomError::UnsupportedInput`
  naming the new `SEAM_SIDE_UNSUPPORTED`, never answered with the
  outgoing frame.
- Conformance: both sides agree with the side-less answer off a seam and
  keep its refusals; at a polyline corner and a grade break each side
  reads its own piece, also within the seam tolerance for a provider that
  reads sides; a side a provider does not read is refused as unsupported.

### Changed

- `CurveEvaluator::frame_at` documents the axes the reference provider
  has always returned (#242): `x` the tangent, `y` up, `z` to the right.
  The text used to call `z` up. No behaviour changes; `axiolid-evaluate`
  pins the axes on a level line, a grade and a banked curve.

## [0.3.2] - 2026-10-03

### Changed

- `CurveEvaluator::frame_at` documents that a curve carrying its own roll
  (`Curve3::Banked`, #240) is framed by it: its section frame, rolled about
  the tangent under the curve's bank convention (ADR 0081).

## [0.3.1] - 2026-09-28

### Changed

- The crates.io page is this crate's own `README.md`, with links to its
  API documentation, its reference page and the source (ADR 0078).

