# Changelog

All notable changes to this crate are documented in this file.

The format follows [Keep a Changelog](https://keepachangelog.com/en/1.1.0/),
and this crate adheres to [Semantic Versioning](https://semver.org/).
Pre-1.0: the minor version is the breaking-change slot, per Cargo's own
caret rule for `0.x` versions.

## [Unreleased]

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

