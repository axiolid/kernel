# Changelog

All notable changes to this crate are documented in this file.

The format follows [Keep a Changelog](https://keepachangelog.com/en/1.1.0/),
and this crate adheres to [Semantic Versioning](https://semver.org/).
Pre-1.0: the minor version is the breaking-change slot, per Cargo's own
caret rule for `0.x` versions.

## [Unreleased]

## [0.3.3] - 2026-11-03

### Fixed

- `self_intersections` no longer reports disjoint, nearly coplanar
  triangles as intersecting after a rotation has rounded their vertices off
  their common plane. The triangle-triangle test was exact only in its side
  tests and fell back to a float `n_p x n_q` direction, which is noise for
  nearly parallel planes. Every decision is now a certified `orient3d` or
  `orient2d` sign: an edge of one triangle against the closed other triangle,
  and for exactly coplanar pairs a separating-edge test on the exact
  projection (#399).
- Coplanar triangles that overlap without a proper edge crossing or a vertex
  strictly inside the other (equal triangles, a vertex on an edge with shared
  area) are now reported; touching along an edge or at a vertex still is not
  (#399).
- Zero-area triangles are never reported as intersecting; the previous
  test reported a needle piercing another triangle. `diagnose` still flags
  such a triangle as `DegenerateElement` (#399).

### Added

- The package now ships `LICENSE` and `LICENSE-COMMERCIAL` (symlinks to the
  repository root), so the licence texts travel with the archive (#358).

## [0.3.2] - 2026-10-28

### Changed

- The licence is now `MPL-2.0 OR LicenseRef-Commercial` (was `MPL-2.0`):
  a commercial licence is available from the copyright holder
  (`LICENSE-COMMERCIAL`, ADR 0089, #358).
- The minimum supported Rust version is now 1.99 (was 1.88), and the
  workspace builds on Rust 1.99.1 (#349).

## [0.3.1] - 2026-09-28

### Changed

- The crates.io page is this crate's own `README.md`, with links to its
  API documentation, its reference page and the source (ADR 0078).

## [0.3.0] - 2026-09-23

### Added

- Repairs carry corner-indexed channels: weld leaves them untouched (they index values, not positions, so a seam is lossless), dropping and flipping triangles move their entries (#112).

### Fixed

- Repairs keep attribute channels and normals in step with the geometry they rewrite (#114). Weld compacts per-vertex channels and normals; a seam drops the channel by name, a hard edge switches normals to corner-indexed. Dropping or flipping triangles moves corner-indexed normals with them.
- `RepairReport::attribute_fates` names every input channel's fate.
