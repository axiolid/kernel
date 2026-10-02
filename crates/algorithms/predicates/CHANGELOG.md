# Changelog

All notable changes to this crate are documented in this file.

The format follows [Keep a Changelog](https://keepachangelog.com/en/1.1.0/),
and this crate adheres to [Semantic Versioning](https://semver.org/).
Pre-1.0: the minor version is the breaking-change slot, per Cargo's own
caret rule for `0.x` versions.

## [Unreleased]

## [0.3.3] - 2026-10-02

### Added

- `in_diametral_sphere` and `in_diametral_sphere_filter` (#126): is a point
  inside the smallest sphere through a triangle? For a point in the
  triangle's plane this is the coplanar in-circle test in 3D, which a 3D
  Delaunay triangulation needs for points in the plane of a hull face. A
  running-error filter, then exact expansions; `Uncertain` outside the
  exactly evaluable range (non-zero coordinates beyond `[2^-100, 2^100]`).

### Fixed

- `insphere_filter` bounded its rounding error by the rounded 3x3 minors
  instead of by the absolute values of the elementary products (Shewchuk's
  permanent). When a minor cancelled, the bound was far too small, and the
  filter certified a non-zero sign for exactly cospherical points -- five
  lattice points of one sphere, four of them on a thin tetrahedron, were
  reported strictly outside (#126). The permanent is now Shewchuk's.

### Changed

- `insphere`'s exact path first tries `i128` arithmetic when every
  coordinate difference is exact and all lie on one dyadic grid spanning
  fewer than 20 bits (integer lattices, for instance), and falls back to the
  expansions otherwise. Exactly cospherical lattice points reach the exact
  path on every call; this keeps it from allocating there.

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
