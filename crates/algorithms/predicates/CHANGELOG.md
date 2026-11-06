# Changelog

All notable changes to this crate are documented in this file.

The format follows [Keep a Changelog](https://keepachangelog.com/en/1.1.0/),
and this crate adheres to [Semantic Versioning](https://semver.org/).
Pre-1.0: the minor version is the breaking-change slot, per Cargo's own
caret rule for `0.x` versions.

## [Unreleased]

### Fixed

- `incircle` and `insphere` returned a wrong sign, certified as exact, for
  three classes of input, in every published release; each now returns the
  exact sign for every finite input (#430):
  - Coordinate differences of about 1e77 or more (`incircle`) or 1e61
    (`insphere`): the exact path's degree-4 and degree-5 products
    overflowed, and a point clearly inside or outside the circle was
    reported as a certain Zero (`(s, 0)`, `(0, s)`, `(-s, 0)` against
    `(0, 0)` and `(2s, 2s)` for s = 1e80, 1e100, 1e160).
  - Very small differences (s = 1e-100 and 1e-160 failed, 1e-80 did
    not): the products underflowed and the same Zero was reported.
  - A NaN or infinite coordinate: it returned a certain sign; it now
    returns `Uncertain`, as `orient3d` does.

  The filters decline when their permanent is outside [1e-250, 1e300]. The
  expansion stage of the exact path runs only where every difference lies in
  [1e-60, 1e60] (`incircle`) or [1e-45, 1e45] (`insphere`) and every nonzero
  rounding tail is at least as large as the lower bound; anything else goes
  to a fixed-size dyadic evaluation over the original coordinates (133 and
  166 words, the accumulator `orient2d` and `orient3d` share), which forms
  no difference and no rounded product. Cost, measured with the committed
  `predicates` bench (three runs each, ns per call; before the change, the
  "huge" scene is the wrong Zero of the overflowing path): random points
  15.6 to 16.1 before, 16.1 to 16.9 after; exactly cocircular lattice points
  (2^-10 grid, every call escalates) 911 to 951 before, 923 to 1027 after;
  the same points scaled to about 1e80 1196 to 1209 before, 5330 to 5440
  after, now exact.
- `StaticFilter::orient3d` certified a wrong sign when the declared range was
  below about 1e-103, in every published release: the precomputed error bound
  underflowed to zero, so every nonzero rounded determinant was accepted
  (#429; 540 wrong signs in 3.2M fuzzed cases). `StaticFilter::new` now returns
  `None` when a derived bound is below the normal range or not finite, the
  same refusal it already gave an overflowing bound. The module no longer
  claims that the static filter never certifies a sign the dynamic one would
  reject: it certifies only exact signs, and the dynamic `orient2d` filter
  declines more (outside [1e-90, 1e90]) since #404.
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
