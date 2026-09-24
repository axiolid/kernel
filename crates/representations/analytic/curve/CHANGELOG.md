# Changelog

All notable changes to this crate are documented in this file.

The format follows [Keep a Changelog](https://keepachangelog.com/en/1.1.0/),
and this crate adheres to [Semantic Versioning](https://semver.org/).
Pre-1.0: the minor version is the breaking-change slot, per Cargo's own
caret rule for `0.x` versions.

## [Unreleased]

### Added

- `Curve2::Sinusoid(Sinusoid2)`: the graph `v = mean + a cos(t) + b sin(t)`,
  the exact pcurve of a plane's cut across a cylinder in its (angle, height)
  parameters (ADR 0071). The parameter is the first coordinate. Additive:
  `Curve2` is `#[non_exhaustive]`.
