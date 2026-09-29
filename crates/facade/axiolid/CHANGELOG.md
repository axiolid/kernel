# Changelog

All notable changes to this crate are documented in this file.

The format follows [Keep a Changelog](https://keepachangelog.com/en/1.1.0/),
and this crate adheres to [Semantic Versioning](https://semver.org/).
Pre-1.0: the minor version is the breaking-change slot, per Cargo's own
caret rule for `0.x` versions.

## [Unreleased]

## [0.4.0] - 2026-09-28

### Changed

- **Breaking:** requires `axiolid-ray-mesh` 0.4, so `axiolid::ray_mesh`
  (re-exported under the ray features) carries the new
  `RayMeshError::TriangleIndexOutOfRange`, and `RayIndex` queries refuse an
  out-of-range candidate instead of skipping it. An exhaustive `match` on
  `RayMeshError` needs the new arm.
- The crates.io page is this crate's own `README.md`, with links to its
  API documentation, its reference page and the source (ADR 0078).

