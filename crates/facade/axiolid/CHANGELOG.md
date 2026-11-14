# Changelog

All notable changes to this crate are documented in this file.

The format follows [Keep a Changelog](https://keepachangelog.com/en/1.1.0/),
and this crate adheres to [Semantic Versioning](https://semver.org/).
Pre-1.0: the minor version is the breaking-change slot, per Cargo's own
caret rule for `0.x` versions.

## [Unreleased]

## [0.5.1] - 2026-11-14

### Added

- `Application::hidden_lines` draws perspective views: the portable
  hidden-line provider claims them (ADR 0107, #387). The perspective claim
  has its own conformance suite, run once per process and only by a call
  whose view is perspective; a failure refuses the call as the case
  `perspective claim` of the hidden-line report. Visible faces still refuse a
  perspective view (`Unsupported`).

## [0.5.0] - 2026-11-13

### Added

- Hidden-line drawings and visible faces join the `discrete` (contract) and
  `application` (dispatch and exact scalar providers) bundles (#371,
  ADR 0105). `Application` gains `hidden_lines`, `audit_hidden_lines` and
  `visible_faces`; `ApplicationBuilder` gains `with_portable_hidden_lines` and
  `with_portable_visible_faces`, and `Application::portable()` selects both.
  The descriptor advertises `MESH_HIDDEN_LINES` and `MESH_VISIBLE_FACES`
  (`Profile2d`, `ToleranceBounded`). Each provider's conformance suite runs
  once per process, on the first call that reaches it; a failed verdict
  refuses the call.
- `IntegrationDescriptor::enabled_features` reports `mesh-hidden-lines`,
  `dispatch-mesh-hidden-lines` and `hidden-lines-provider`.

### Changed

- **Breaking:** `ApplicationErrorSource` and `ProviderConformanceReport` are
  `#[non_exhaustive]`. `ApplicationErrorSource` gains `HiddenLines` and
  `ProviderConformance`; `ProviderConformanceReport` gains `MeshHiddenLines`
  and `MeshVisibleFaces`. A `match` on either needs a wildcard arm (#371,
  ADR 0105).

- **Breaking:** requires `axiolid-overlay` 0.4, so `axiolid::overlay` carries
  `OverlayError` as `#[non_exhaustive]` with the new `HolesOverlap` and
  `PolygonsOverlap`, and `Region::new` refuses rings that cross into each
  other. A `match` on `OverlayError` needs a wildcard arm (#448, ADR 0104).

## [0.4.7] - 2026-11-03

### Added

- The package now ships `LICENSE` and `LICENSE-COMMERCIAL` (symlinks to the
  repository root), so the licence texts travel with the archive (#358).

- `mesh_hidden_lines` re-exports `ScalarVisibleFaces` under the
  `hidden-lines-provider` feature (ADR 0094, #383).

## [0.4.6] - 2026-10-29

### Added

- `mesh_section` exports `CutRegion`, `CutRegionLimits` and the cut-region
  refusal names, so `Application::section_mesh_request` can return the cut of
  a solid as closed regions (ADR 0092, #372).

## [0.4.5] - 2026-10-29

### Added

- `Application::section_mesh_request` sections open meshes when the request
  has `SectionRole::Sheet` (ADR 0091, #359).

## [0.4.4] - 2026-10-28

### Added

- `Application::section_mesh_request` sections a closed oriented mesh and
  reports faces lying in the plane as regions with a side (ADR 0091, #359).

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

