# Changelog

All notable changes to this crate are documented in this file.

The format follows [Keep a Changelog](https://keepachangelog.com/en/1.1.0/),
and this crate adheres to [Semantic Versioning](https://semver.org/).
Pre-1.0: the minor version is the breaking-change slot, per Cargo's own
caret rule for `0.x` versions.

## [Unreleased]

### Added

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
