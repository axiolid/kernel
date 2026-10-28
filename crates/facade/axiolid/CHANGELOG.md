# Changelog

All notable changes to this crate are documented in this file.

The format follows [Keep a Changelog](https://keepachangelog.com/en/1.1.0/),
and this crate adheres to [Semantic Versioning](https://semver.org/).
Pre-1.0: the minor version is the breaking-change slot, per Cargo's own
caret rule for `0.x` versions.

## [Unreleased]

## [0.4.3] - 2026-10-28

### Added

- `mesh-hidden-lines` and `dispatch-mesh-hidden-lines` features, with the
  `axiolid::mesh_hidden_lines` module: the orthographic hidden-line
  drawing contract and its registry (ADR 0090, #138). Opt-in; the bundles
  do not include it.
- `hidden-lines-provider` feature: `axiolid::mesh_hidden_lines::ScalarHiddenLines`,
  the exact scalar provider (`axiolid-hidden-lines`, #138).

## [0.4.2] - 2026-10-28

### Changed

- The licence is now `MPL-2.0 OR LicenseRef-Commercial` (was `MPL-2.0`):
  a commercial licence is available from the copyright holder
  (`LICENSE-COMMERCIAL`, ADR 0089, #358).
- The minimum supported Rust version is now 1.99 (was 1.88), and the
  workspace builds on Rust 1.99.1 (#349).

## [0.4.1] - 2026-10-10

### Added

- A `serde` feature (#267, ADR 0085) forwarding to the `serde` feature of
  every representation package already in the build; with `model` it
  brings the geometry graph's versioned JSON and CBOR wire format
  (`axiolid::model::wire`). It enables no representation by itself.

## [0.4.0] - 2026-09-28

### Changed

- **Breaking:** requires `axiolid-ray-mesh` 0.4, so `axiolid::ray_mesh`
  (re-exported under the ray features) carries the new
  `RayMeshError::TriangleIndexOutOfRange`, and `RayIndex` queries refuse an
  out-of-range candidate instead of skipping it. An exhaustive `match` on
  `RayMeshError` needs the new arm.
- The crates.io page is this crate's own `README.md`, with links to its
  API documentation, its reference page and the source (ADR 0078).

