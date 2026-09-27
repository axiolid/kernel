# Changelog

All notable changes to this crate are documented in this file.

The format follows [Keep a Changelog](https://keepachangelog.com/en/1.1.0/),
and this crate adheres to [Semantic Versioning](https://semver.org/).
Pre-1.0: the minor version is the breaking-change slot, per Cargo's own
caret rule for `0.x` versions.

## [Unreleased]

## [0.3.1] - 2026-09-27

### Changed

- `BSplineSurface` is defined in `axiolid-curve` (so a curve traced on a
  B-spline surface can carry its carrier, ADR 0077) and re-exported here
  unchanged: same fields, same derives. Requires `axiolid-curve` 0.3.1.
