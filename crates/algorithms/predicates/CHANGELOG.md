# Changelog

All notable changes to this crate are documented in this file.

The format follows [Keep a Changelog](https://keepachangelog.com/en/1.1.0/),
and this crate adheres to [Semantic Versioning](https://semver.org/).
Pre-1.0: the minor version is the breaking-change slot, per Cargo's own
caret rule for `0.x` versions.

## [Unreleased]

### Fixed

- `orient2d` returned a wrong sign, certified as exact, for four classes of
  input, in every published release; it now returns the exact sign for every
  finite input (#404):
  - Coordinate differences of about 1e154 or more: the products overflowed
    (a debug build panicked, a release build returned a wrong sign).
  - Coordinates of very different magnitudes (a ratio of about 1e32 or
    more, e.g. `(1e-100, 1e-100)`, `(1.5e-100, 1e-100)`, `(-4.5e50, -1)`
    gave Zero for a Negative triple): the expansion path lost error terms.
  - Differences near 1e-155: the filter's products underflowed and it
    proved a wrong sign (a Negative for a truly Positive triple).
  - A NaN or infinite coordinate was reported as a certain Zero; it now
    returns `Uncertain`, as `orient3d` does.

  The exact path now evaluates the determinant over the original coordinates,
  as an exact sum of expansions of the coordinate differences and their
  rounding tails when those are in range (the usual exactly collinear case),
  otherwise in a fixed-size dyadic accumulator. `orient2d_filter` declines
  differences outside [1e-90, 1e90] as `Uncertain`. The cost is real for
  exactly collinear triples: per call, 24 ns on a full-mantissa millimetre
  grid whose differences are exact (12 ns before) and 50 ns on one crossing
  binades, where 45% of the triples have a rounding difference (12 ns
  before); measured with the committed `predicates` bench. Non-degenerate
  input is unaffected (about 15 ns).

### Added

- The package now ships `LICENSE` and `LICENSE-COMMERCIAL` (symlinks to the
  repository root), so the licence texts travel with the archive (#358).

## [0.3.4] - 2026-10-28

### Changed

- The licence is now `MPL-2.0 OR LicenseRef-Commercial` (was `MPL-2.0`):
  a commercial licence is available from the copyright holder
  (`LICENSE-COMMERCIAL`, ADR 0089, #358).
- The minimum supported Rust version is now 1.99 (was 1.88), and the
  workspace builds on Rust 1.99.1 (#349).

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
