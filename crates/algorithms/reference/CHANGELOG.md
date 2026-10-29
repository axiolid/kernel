# Changelog

All notable changes to this crate are documented in this file.

The format follows [Keep a Changelog](https://keepachangelog.com/en/1.1.0/),
and this crate adheres to [Semantic Versioning](https://semver.org/).
Pre-1.0: the minor version is the breaking-change slot, per Cargo's own
caret rule for `0.x` versions.

## [Unreleased]

## [0.3.14] - 2026-10-29

### Fixed

- Cut rings are decided on exact source points: crossing, contact,
  containment and outer or hole no longer use the returned coordinates, which
  the projection to a turned frame and the edge interpolation round. Shells
  that touch in the source are no longer refused as crossing under a turned
  frame or an oblique section. An overlap that is in the input is still
  refused. The returned coordinates are unchanged (ADR 0092, #377).

### Changed

- Depends on `axiolid-exact` for the exact points of the cut.

## [0.3.13] - 2026-10-29

### Added

- `ScalarSection::section_request` traces the cut of a solid as closed regions
  when asked (`SectionRequest::with_cut_regions`): region edges that the cut
  lies across are found with exact source signs, linked to the cut segments
  into rings, and nested exactly on the returned coordinates. Pinches and
  crossing rings are refused by name, contacts between rings are decided
  exactly (a ring crossing another through a vertex, a repeated vertex and a
  spike are refused); a sheet request is invalid input (ADR 0092, #372).

## [0.3.12] - 2026-10-29

### Added

- `ScalarSection::section_request` sections `SectionRole::Sheet` meshes: a
  chain ends on the mesh boundary at a crossed boundary edge or a boundary
  vertex (`ChainEnd::MeshBoundary`), a boundary edge in the plane next to an
  off-plane triangle is a touching edge, a sheet in the plane is a region, and
  a sheet folded flat onto itself is refused as in-plane faces of both sides
  (ADR 0091, #359). The request scratch bound is 1536 bytes per element.

## [0.3.11] - 2026-10-28

### Added

- `ScalarSection::section_request`: faces lying in the plane become regions
  with a side and holes, found by exact orientation predicates on the source
  coordinates; cut curves are transversal chains, open ones ending on a region,
  oriented and simplified exactly (no tolerance);
  touching edges are counted. Overlapping in-plane faces and saddle vertices
  are refused by name (ADR 0091, #359).

## [0.3.10] - 2026-10-28

### Changed

- The licence is now `MPL-2.0 OR LicenseRef-Commercial` (was `MPL-2.0`):
  a commercial licence is available from the copyright holder
  (`LICENSE-COMMERCIAL`, ADR 0089, #358).
- The minimum supported Rust version is now 1.99 (was 1.88), and the
  workspace builds on Rust 1.99.1 (#349).

## [0.3.9] - 2026-10-18

### Fixed

- Requires `axiolid-evaluate` 0.3.17, whose certified frame laws and
  intrinsic positions it re-exports and `axiolid-mesh-compile` 0.3.28 uses
  through it. 0.3.8 accepted any 0.3.12 or later, so a lockfile holding an
  older `axiolid-evaluate` did not build `axiolid-mesh-compile` 0.3.28.

## [0.3.8] - 2026-10-12

### Fixed

- `clash::interference` counted every narrow-phase pair of well-formed
  triangles that do not meet (`TriangleTriangleRelation::Disjoint`) as a
  `degenerate_skips`, through a catch-all arm, so a clash report over-stated
  how much of its input was untestable (two unit boxes 5 apart at a 10 m
  tolerance reported 120 of 144 pairs skipped). Disjoint pairs are now a
  verdict, and the match is exhaustive so a new relation cannot fall
  through again.

## [0.3.7] - 2026-10-04

### Added

- `elevated`: `axiolid-evaluate`'s derivatives, certified chord and
  derivative bounds and grade corners of elevated and banked curves
  (#252), re-exported on the same footing as `curve` and `surface`.

## [0.3.6] - 2026-10-03

### Added

- `station`: `axiolid-evaluate`'s station resolution (#241), re-exported
  on the same footing as `curve` and `surface`.

## [0.3.5] - 2026-10-03

### Added

- Re-exports `axiolid_evaluate::arc_parameter` and
  `axiolid_evaluate::chain` (#239).

## [0.3.4] - 2026-10-03

### Added

- `bound`, re-exported from `axiolid-evaluate` like `curve` and `surface`
  (#232).

### Fixed

- `tessellate_primitive` refuses a cylinder or a cone whose chord budget
  needs more than 4096 segments with `BudgetExceeded`, as it already did
  for spheres and tori, instead of clamping silently to a coarser mesh
  (#232). Their bound is now documented: both are curved one way only, so
  the ring's sagitta is the whole distance from any point of the exact
  surface, caps included, to the mesh, and it gets the whole budget.

## [0.3.3] - 2026-10-02

### Fixed

- `tessellate_primitive` keeps every point of a sphere and a torus within
  the chord budget of the mesh, not only its vertices (#231). Both
  directions of each used to get the whole budget, so their sagittas added
  inside a triangle (a 0.3 m sphere at 1 mm lay 1.99 mm from its mesh, a
  torus 1.90 mm); each now gets half, and the sphere's stack count is
  rounded up rather than down, which had let an odd segment count leave
  the polar step alone above its share. A budget that needs more than
  4096 segments is refused with `BudgetExceeded` instead of clamped.

## [0.3.2] - 2026-09-28

### Changed

- The crates.io page is this crate's own `README.md`, with links to its
  API documentation, its reference page and the source (ADR 0078).

## [0.3.1] - 2026-09-27

### Added

- `tessellate_primitive` meshes `Primitive::Torus` and `Primitive::Wedge`
  (#142). The torus is a grid of planar trapezoids sized by the chord
  budget, round the axis for the outer equator and round the tube for the
  tube; horn and spindle tori, and non-positive or non-finite radii, are
  refused by name. The wedge's faces are planar and shared corners of a
  collapsed top are merged, so a ridge or apex wedge is still a closed,
  outward-wound solid; a reversed or non-finite top range is refused.
