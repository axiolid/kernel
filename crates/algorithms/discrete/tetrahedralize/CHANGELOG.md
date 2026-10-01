# Changelog

All notable changes to this crate are documented in this file.

The format follows [Keep a Changelog](https://keepachangelog.com/en/1.1.0/),
and this crate adheres to [Semantic Versioning](https://semver.org/).
Pre-1.0: the minor version is the breaking-change slot, per Cargo's own
caret rule for `0.x` versions.

## [Unreleased]

### Added

- `Delaunay3` (#126): exact incremental 3D Delaunay tetrahedralization.
  Bowyer-Watson insertion in Hilbert-curve order, located by a stochastic
  visibility walk, with an infinite vertex for the convex hull.
  `from_points`, `insert`, `tetrahedra` (positively oriented, with
  adjacency), `hull_triangles` (outward), `vertex_of`, `dimension`.
  `orient3d` and `insphere` decide every step exactly; cospherical points,
  points in the plane of a hull face (decided by `in_diametral_sphere`) and
  other ties are broken by a symbolic perturbation of the lifted points,
  so the result is unique for a point set whatever the insertion order.
  Exact duplicates are merged (`Insertion::Duplicate`). Refused by name
  (`Delaunay3Error`): fewer than three dimensions, non-finite coordinates,
  and non-zero coordinates outside `[2^-100, 2^100]`.
