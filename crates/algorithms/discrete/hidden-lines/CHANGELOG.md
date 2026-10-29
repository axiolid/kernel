# Changelog

All notable changes to this crate are documented in this file.

The format follows [Keep a Changelog](https://keepachangelog.com/en/1.1.0/),
and this crate adheres to [Semantic Versioning](https://semver.org/).
Pre-1.0: the minor version is the breaking-change slot, per Cargo's own
caret rule for `0.x` versions.

## [Unreleased]

## [0.1.2] - 2026-10-29

### Changed

- Candidate triangles for an edge, and the cap triangles of a cut solid, come from an `axiolid-spatial` BVH over view-space boxes instead of a linear scan. The predicate, the candidate set, its ascending order and `candidate_pairs` are unchanged (#138). Wall clock from `scenario/hidden-lines` (`cargo bench -p axiolid-benchmark --bench scenario -- hidden-lines`, which times the BVH and the linear scan, `ScanHiddenLines`, on each case; a loaded shared machine, so read these as ranges), linear scan / BVH: isometric 20x20 boxes 412 / 389 ms (a reviewer's run: 459 / 426), cut plan 20x20 174 / 140 ms (117 / 100; very noisy), isometric 50x50 3.34 / 2.52 s (3.62 / 2.88), cut plan 50x50 1.59 / 0.76 s (1.67 / 0.82). The gain is small at 400 boxes, because exact evaluation of the candidates dominates there, and 25-30% (isometric) to about 50% (cut plan) at 2,500.

## [0.1.1] - 2026-10-28

### Changed

- The near-plane caps of closed solids occlude (#138, ADR 0090): a
  `Solid` that the near plane cuts no longer fails as `UnsupportedInput`.
  A point with `W < near` (farther than the plane) is hidden when the point
  of the plane above it lies in the solid, by exact ray parity with a
  symbolic perturbation of the image point. A face of a cut solid in the
  near plane is drawn, not refused.

## [0.1.0] - 2026-10-28

### Added

- `ScalarHiddenLines`, the exact scalar provider for `MeshHiddenLines`
  (ADR 0090, #138): boundary, silhouette and crease edges of triangle
  meshes split into visible and hidden pieces by exact interval and dyadic
  decisions, with endpoints that are the exact points rounded once. A near
  clipping plane that cuts a `Solid` mesh is refused as `UnsupportedInput`
  naming it.
