# Changelog

All notable changes to this crate are documented in this file.

The format follows [Keep a Changelog](https://keepachangelog.com/en/1.1.0/),
and this crate adheres to [Semantic Versioning](https://semver.org/).
Pre-1.0: the minor version is the breaking-change slot, per Cargo's own
caret rule for `0.x` versions.

## [Unreleased]

## [0.3.4] - 2026-09-28

### Changed

- The crates.io page is this crate's own `README.md`, with links to its
  API documentation, its reference page and the source (ADR 0078).

## [0.3.3] - 2026-09-27

### Added

- `topology` (#144): every connected component of a two-manifold triangle
  mesh, classified exactly from its connectivity -- counts, Euler
  characteristic, boundary loops, orientability and consistent winding, and
  the surface (`SurfaceKind::Orientable { genus }` or
  `NonOrientable { crosscaps }`). A closed orientable component also gets a
  basis of its first homology: `2g` simple closed edge loops, by the
  tree-cotree construction. Meshes with an edge on three or more triangles,
  or a vertex whose triangles form several fans, are refused
  (`TopologyError::NonManifold`). `genus` is unchanged.

## [0.3.2] - 2026-09-27

### Added

- `detect_planes` (#131): the planar regions of a triangle mesh, grown
  over shared edges within an angle and a distance (`PlaneTolerance`) and
  then certified. Each `DetectedPlane` gives its triangles, a point and
  unit normal, a `deviation` that is a proven upper bound on every member
  corner's distance from that plane (regions are peeled until it is
  within the requested distance), whether the region is exactly coplanar
  (by `orient3d`), and its area. Largest first; deterministic.

## [0.3.1] - 2026-09-27

### Added

- `intersection_volume`, `difference_volume` and `enclosed_volume` (#183):
  the volume two closed triangle meshes share, the volume of the first
  outside the second, and a mesh's own volume, each as a `VolumeInterval`
  certified to contain the true value and never negative. No boolean is
  built: every face pair whose shadows overlap contributes the integral of
  the lower of the two planes over the overlap, with sides of lines and
  the lower plane decided exactly and the arithmetic rounded outward.
  Coplanar faces and touching bodies are exact cases. Open, non-manifold,
  self-intersecting or non-finite meshes are refused with the operand
  named (`OverlapError`); either winding is accepted.
