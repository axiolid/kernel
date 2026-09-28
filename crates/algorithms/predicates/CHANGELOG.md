# Changelog

All notable changes to this crate are documented in this file.

The format follows [Keep a Changelog](https://keepachangelog.com/en/1.1.0/),
and this crate adheres to [Semantic Versioning](https://semver.org/).
Pre-1.0: the minor version is the breaking-change slot, per Cargo's own
caret rule for `0.x` versions.

## [Unreleased]

## [0.3.2] - 2026-09-28

### Changed

- The crates.io page is this crate's own `README.md`, with links to its
  API documentation, its reference page and the source (ADR 0078).

## [0.3.1] - 2026-09-27

### Fixed

- `incircle` and `insphere`: the exact fallback rounded the coordinate
  differences to `f64` before its exact expansion arithmetic, so on nearly
  cocircular (cospherical) points whose differences do not fit an `f64` --
  exactly where the filter hands over -- it could return the wrong sign.
  Delaunay flips driven by it cycled for ever (#190). The differences are
  now exact two-term expansions and every product after them is an
  expansion product; checked against an exact dyadic determinant.
