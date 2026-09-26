# Changelog

All notable changes to this crate are documented in this file.

The format follows [Keep a Changelog](https://keepachangelog.com/en/1.1.0/),
and this crate adheres to [Semantic Versioning](https://semver.org/).
Pre-1.0: the minor version is the breaking-change slot, per Cargo's own
caret rule for `0.x` versions.

## [Unreleased]

### Added

- `Curve2::QuadraticGraph(QuadraticGraph2)` and
  `Curve3::RuledSection(RuledSection3)` (#119, ADR 0076): one root branch of
  `a(t) v^2 + b(t) v + c(t) = 0` with degree-2 trigonometric coefficients
  (`Trig2`, `Branch`), and the same curve lifted onto a cylinder,
  elliptical cylinder or cone (`RuledCarrier`). The exact pcurve and edge
  of a quadric's cut across a ruled surface.

- `Curve2::Sinusoid(Sinusoid2)`: the graph `v = mean + a cos(t) + b sin(t)`,
  the exact pcurve of a plane's cut across a cylinder in its (angle, height)
  parameters (ADR 0071). The parameter is the first coordinate. Additive:
  `Curve2` is `#[non_exhaustive]`.
