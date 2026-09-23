# Changelog

All notable changes to this crate are documented in this file.

The format follows [Keep a Changelog](https://keepachangelog.com/en/1.1.0/),
and this crate adheres to [Semantic Versioning](https://semver.org/).
Pre-1.0: the minor version is the breaking-change slot, per Cargo's own
caret rule for `0.x` versions.

## [Unreleased]

### Added

- The pairwise boolean carries attribute channels (#116). Each result triangle's source triangle is tracked through the CSG core, so a corner that is a source corner copies its value and a corner on a cut is derived in its source triangle under the channel's `Blend`. Output is corner-indexed; faces from a tool without the channel are `UNMAPPED`. Fates: `Preserved` when nothing was derived, `Interpolated` otherwise, `Dropped(NotBlendable)` for a `Blend::None` channel a cut needs. The analytic box path and the grouped/tree batch paths still drop channels (`ProviderLimitation`).

### Fixed

- `dedupe_edge` pushed per-face data indexed by a vertex id when pinching a vertex, desynchronising face normals and provenance from the faces (inherited from upstream; Manifold pushes only per-vertex data there). The `ProviderLimitation` comment claiming the boolean returns positions only, stale since ADR 0047, is corrected.
