# Changelog

All notable changes to this crate are documented in this file.

The format follows [Keep a Changelog](https://keepachangelog.com/en/1.1.0/),
and this crate adheres to [Semantic Versioning](https://semver.org/).
Pre-1.0: the minor version is the breaking-change slot, per Cargo's own
caret rule for `0.x` versions.

## [Unreleased]

## [0.3.7] - 2026-11-03

### Added

- `PlaneSectionEvidence::in_plane_crossings` and
  `with_in_plane_crossings`: the distinct points, other than source vertices,
  where in-plane source edges cross (ADR 0093, #366).
- Fixtures for stacked boxes, tubes and an L-step block, whose faces coincide
  in the plane, and for `Both` regions as cut regions.

### Changed

- Every `InPlaneSide::Both` region is a cut region, equal to it bit for bit
  (ADR 0093, #366). The region conformance suite checks it, and the registry
  checks it instead of refusing cut regions beside a `Both` region. A
  provider that reports a `Both` region and refuses only the cut request now
  fails the conformance suite; it must either refuse the plain request too or
  return the cut regions.
- `in_plane_pairs` and `RegionLimits::max_in_plane_pairs` count the candidate
  pairs examined for the in-plane faces (segment, piece-triangle and hole
  containment pairs of a merge), not only overlap tests. The count grows for
  the same mesh (68,104 to 105,600 on the 8,192-triangle bench grid), so a
  `max_in_plane_pairs` sized to the old counts can now give
  `BudgetExceeded`.
- Documentation: the classes of the plane (`Below`, `Above`, `Both`), that
  regions have pairwise disjoint interiors, and that a provider may still
  refuse overlapping in-plane faces by name.

## [0.3.6] - 2026-10-29

### Changed

- Documentation of `CutRegion`: the invariant that rings never cross, and the
  exactness section, now say that the decisions are made on exact source
  points (source vertices in the plane and exact crossings of source edges),
  and that the returned points are those computed in binary64, so rings that
  touch may cross by rounding error (ADR 0092, #377).

## [0.3.5] - 2026-10-29

### Added

- `SectionRequest::with_cut_regions(CutRegionLimits)` and
  `PlaneSection::cut_regions`: the cut of a solid, where its material lies on
  both sides of the plane, as closed `CutRegion`s (outer ring, holes) whose
  rings run through the chain and region vertices bit for bit. Refusals
  `REFUSED_PINCHED_CUT_REGION` and `REFUSED_CROSSING_CUT_RINGS`; evidence
  fields `cut_regions`, `cut_region_vertices` and `cut_pairs`. Additive
  (ADR 0092, #372).
- `RegionConformanceSuite` runs every case again with cut regions requested,
  adds the L-step, wall at a sill, wall at a head, stepped slab, pocket floor
  and nested tubes fixtures, and reports `cut_regions_claimed`;
  `RegionConformanceFailure::WrongCutRegions`. A provider that answers the
  first solid case without cut regions predates them and is not claimed; one
  that returns them and later drops them fails.

### Fixed

- The documentation said a closed chain runs counter-clockwise; a closed chain
  of a solid keeps the material on its left, so one around a hole runs
  clockwise.

## [0.3.4] - 2026-10-29

### Added

- `SectionRole` (`Solid`, `Sheet`) and `SectionRequest::role` with
  `with_role`: a sheet is an open surface (terrain, a surface body), sectioned
  with boundary edges allowed. `ChainEnd::MeshBoundary` and the evidence fields
  `role` and `boundary_edges` (`PlaneSectionEvidence::with_role`). The default
  stays `Solid`, today's semantics. Additive (ADR 0091, #359).
- `RegionConformanceSuite` adds the terrain cut by a vertical plane and
  through vertices, a sheet lying in the plane in both windings, an annulus
  sheet and a saddle (#359). `RegionConformanceReport::sheets_claimed` is
  `false` when the provider refuses the first sheet case with `Unsupported`;
  the sheet cases are then skipped.

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

