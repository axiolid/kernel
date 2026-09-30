# Changelog

All notable changes to this crate are documented in this file.

The format follows [Keep a Changelog](https://keepachangelog.com/en/1.1.0/),
and this crate adheres to [Semantic Versioning](https://semver.org/).
Pre-1.0: the minor version is the breaking-change slot, per Cargo's own
caret rule for `0.x` versions.

## [Unreleased]

## [0.3.3] - 2026-09-30

### Added

- Wachspress and discrete harmonic coordinates (#143, ledger row H5):
  `wachspress_coordinates2` and `discrete_harmonic_coordinates2`, both for a
  strictly convex polygon (refused by name otherwise, `BarycentricError::NotConvex`),
  positive and smooth inside, boundary-linear on the edges, and refused as
  `BarycentricError::Undefined` at a pole (a denominator vanishing on the
  line through a non-adjacent edge, or discrete harmonic's own collinear
  fan angle). And `mean_value_coordinates3` (Floater, Kos and Reimers 2005),
  extending mean-value coordinates to a closed triangle mesh in space, one
  weight per vertex, exact at vertices and on a face's own plane inside its
  triangle.

## [0.3.2] - 2026-09-28

### Changed

- The crates.io page is this crate's own `README.md`, with links to its
  API documentation, its reference page and the source (ADR 0078).

## [0.3.1] - 2026-09-27

### Added

- Barycentric coordinates (#143): `triangle_barycentric2`,
  `triangle_barycentric3` (for the point's projection onto the triangle's
  plane) and `tetrahedron_barycentric`, exact at corners; and
  `mean_value_coordinates2` for simple polygons, convex or not, which
  interpolate the boundary linearly and reproduce points inside. Shapes
  thinner than the linear tolerance, non-simple polygons and points where
  mean-value weights cancel are refused with `BarycentricError`.
