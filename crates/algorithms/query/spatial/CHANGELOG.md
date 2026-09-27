# Changelog

All notable changes to this crate are documented in this file.

The format follows [Keep a Changelog](https://keepachangelog.com/en/1.1.0/),
and this crate adheres to [Semantic Versioning](https://semver.org/).
Pre-1.0: the minor version is the breaking-change slot, per Cargo's own
caret rule for `0.x` versions.

## [Unreleased]

## [0.3.1] - 2026-09-27

### Added

- Barycentric coordinates (#143): `triangle_barycentric2`,
  `triangle_barycentric3` (for the point's projection onto the triangle's
  plane) and `tetrahedron_barycentric`, exact at corners; and
  `mean_value_coordinates2` for simple polygons, convex or not, which
  interpolate the boundary linearly and reproduce points inside. Shapes
  thinner than the linear tolerance, non-simple polygons and points where
  mean-value weights cancel are refused with `BarycentricError`.
