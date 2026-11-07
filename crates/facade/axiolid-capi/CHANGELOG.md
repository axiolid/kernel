# Changelog

All notable changes to this crate are documented in this file.

The format follows [Keep a Changelog](https://keepachangelog.com/en/1.1.0/),
and this crate adheres to [Semantic Versioning](https://semver.org/).
Pre-1.0: the minor version is the breaking-change slot, per Cargo's own
caret rule for `0.x` versions.

## [Unreleased]

### Changed

- **Breaking for Rust callers only.** The public C-mirroring enums
  `AxiolidStatus`, `AxiolidOperation`, `AxiolidRepresentation`,
  `AxiolidExactness` and `AxiolidGeometryKind` gain variants (`MeshHiddenLines`,
  `MeshVisibleFaces`, `HiddenLines`) and are now `#[non_exhaustive]`, so later
  C-ABI-additive variants are no longer Rust breaks. The C ABI is unchanged
  and additive (`abi_patch` 1).

### Added

- Visible faces through the C ABI (#386): `axiolid_v0_4_visible_faces` with
  `AxiolidVisibleFaceLimits`; the result is read with
  `axiolid_v0_4_visible_faces_counts`, `_copy` (regions, rings, points and
  triangle indices as flat offset-indexed tables, nothing written unless
  every capacity covers its count) and `_evidence`. `AxiolidFaceSource`
  distinguishes group regions from near-plane caps, `AxiolidGeometryKind`
  gains `VisibleFaces`, and the capability table lists `MeshVisibleFaces`.
  Arguments, refusals and mesh ids behave as for hidden lines. Existing
  symbols and layouts are unchanged.
- Hidden lines through the C ABI (#473): `axiolid_v0_4_hidden_lines` with
  `AxiolidView`, `AxiolidViewMesh` (roles, group modes including planar),
  `AxiolidEdgeSelection` (including group boundaries) and
  `AxiolidHiddenLineLimits`; the result is read with
  `axiolid_v0_4_hidden_lines_count`, `_copy` and `_evidence`.
  `axiolid_v0_4_context_last_error_mesh` returns the id of the input mesh a
  refusal names. `AxiolidOperation` gains `MeshHiddenLines` and
  `MeshVisibleFaces`, `AxiolidGeometryKind` gains `HiddenLines`, and the
  capability table lists `MeshHiddenLines`. All arguments are validated at the
  boundary. `abi_patch` is 1; existing symbols and layouts are unchanged.

## [0.3.4] - 2026-11-07

### Added

- The package now ships `LICENSE` and `LICENSE-COMMERCIAL` (symlinks to the
  repository root), so the licence texts travel with the archive (#358).

## [0.3.3] - 2026-10-28

### Changed

- The package metadata names rheonda UG (haftungsbeschränkt) as author, and commercial
  licences are available from rheonda UG (haftungsbeschränkt) (`LICENSE-COMMERCIAL`, ADR 0089).

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

