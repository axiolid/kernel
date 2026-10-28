# Changelog

All notable changes to this crate are documented in this file.

The format follows [Keep a Changelog](https://keepachangelog.com/en/1.1.0/),
and this crate adheres to [Semantic Versioning](https://semver.org/).
Pre-1.0: the minor version is the breaking-change slot, per Cargo's own
caret rule for `0.x` versions.

## [Unreleased]

## [0.3.3] - 2026-10-28

### Added

- `SectionRequest`, `PlaneSection` and the provided `MeshPlaneSection` methods
  `section_request` and `request_scratch_requirement`: faces lying in the
  cutting plane come back as `InPlaneRegion`s with an `InPlaneSide`, cut curves
  are transversal `SectionChain`s (open ones end on a region), and touching
  edges are counted in `PlaneSectionEvidence`. A provider that does not
  implement the method refuses with `Unsupported`. The refusal names are
  `REFUSED_OVERLAPPING_IN_PLANE_FACES`, `REFUSED_IN_PLANE_FACES_ON_BOTH_SIDES`
  and `REFUSED_SADDLE_VERTEX`. Additive (ADR 0091, #359).
- `conformance::RegionConformanceSuite` runs hand-derived fixtures for the
  method, among them a cube through its top face, a frame with a hole and a
  parapet with a post (#359).

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

