# Changelog

All notable changes to this crate are documented in this file.

The format follows [Keep a Changelog](https://keepachangelog.com/en/1.1.0/),
and this crate adheres to [Semantic Versioning](https://semver.org/).
Pre-1.0: the minor version is the breaking-change slot, per Cargo's own
caret rule for `0.x` versions.

## [Unreleased]

### Added

- `remesh::remesh` (#149): isotropic remeshing towards a target edge length `L` in the Botsch–Kobbelt scheme CGAL's `isotropic_remeshing` follows. Each iteration splits edges longer than `4/3 L`, collapses edges shorter than `4/5 L` where no edge longer than `4/3 L` results, flips edges to bring valences towards 6 (4 on a boundary), and moves free vertices towards the area-weighted centroid of their triangles in the tangent plane before projecting them onto the input. Works on a `HalfedgeMesh`; deterministic.
- Boundary edges and, by default, edges sharper than 60 degrees (`RemeshOptions::feature_angle`) are protected: never flipped, split on their own segment, collapsed only where their line runs straight, so feature corners stay bit-identical and the Euler characteristic, boundary loops and orientation are kept. Free vertices are projected only onto their own patch of the input (the pieces between protected edges), with a bounding-volume hierarchy per patch.
- A collapse, flip or move that would turn a triangle's normal by a right angle or more, or leave it flatter than a height-to-longest-edge ratio of `1e-3`, is skipped, never forced; so is a collapse or flip that would lower the smallest angle of the triangles it rewrites below 15 degrees.
- `RemeshReport` measures the result: edits made, the fraction of edges within `[4/5 L, 4/3 L]`, mean valence deviation before and after, the smallest angle, and the largest distance of an output vertex from the input.
- `RemeshError` refuses non-manifold, inconsistently wound, ragged or out-of-range input (wrapping `HalfedgeBuildError`), non-finite positions, degenerate input triangles, a non-positive target, a feature angle outside `[0, pi]` and an output over the triangle budget.

## [0.3.1] - 2026-09-28

### Changed

- The crates.io page is this crate's own `README.md`, with links to its
  API documentation, its reference page and the source (ADR 0078).

## [0.3.0] - 2026-09-23

### Fixed

- A refinement that creates no vertex returns the input's channels and normals. It previously reported them `Preserved` and returned a mesh without them.
