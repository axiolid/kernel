# Changelog

All notable changes to this crate are documented in this file.

The format follows [Keep a Changelog](https://keepachangelog.com/en/1.1.0/),
and this crate adheres to [Semantic Versioning](https://semver.org/).
Pre-1.0: the minor version is the breaking-change slot, per Cargo's own
caret rule for `0.x` versions.

## [Unreleased]

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

