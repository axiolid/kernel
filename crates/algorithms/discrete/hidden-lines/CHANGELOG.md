# Changelog

All notable changes to this crate are documented in this file.

The format follows [Keep a Changelog](https://keepachangelog.com/en/1.1.0/),
and this crate adheres to [Semantic Versioning](https://semver.org/).
Pre-1.0: the minor version is the breaking-change slot, per Cargo's own
caret rule for `0.x` versions.

## [Unreleased]

### Changed

- The near-plane caps of closed solids occlude (#138, ADR 0090): a
  `Solid` that the near plane cuts no longer fails as `UnsupportedInput`.
  A point nearer than the plane is hidden when the point of the plane above
  it lies in the solid, by exact ray parity with a symbolic perturbation of
  the image point. A cut solid with a triangle that meets the near plane
  along the image of a drawn edge is refused as `Degenerate`, naming it.

## [0.1.0] - 2026-10-28

### Added

- `ScalarHiddenLines`, the exact scalar provider for `MeshHiddenLines`
  (ADR 0090, #138): boundary, silhouette and crease edges of triangle
  meshes split into visible and hidden pieces by exact interval and dyadic
  decisions, with endpoints that are the exact points rounded once. A near
  clipping plane that cuts a `Solid` mesh is refused as `UnsupportedInput`
  naming it.
