# Changelog

All notable changes to this crate are documented in this file.

The format follows [Keep a Changelog](https://keepachangelog.com/en/1.1.0/),
and this crate adheres to [Semantic Versioning](https://semver.org/).
Pre-1.0: the minor version is the breaking-change slot, per Cargo's own
caret rule for `0.x` versions.

## [Unreleased]

## [0.3.17] - 2026-11-13

### Added

- `bound::ElevatedPositionTable` is re-exported from `axiolid-evaluate`
  0.3.21 (#393): certified points of an elevated or banked curve over an
  intrinsic plan, which `axiolid-mesh-compile` reads for its swept disks.
  Over a chain plan it also reads the tangent and the second derivative
  (`jet3_on`, `jet_errors`, #438).
- `bound::ChainPositionTable` and `bound::SeamSide` are re-exported from
  `axiolid-evaluate` 0.3.21 (#438): certified points, tangents and
  curvature vectors of an arc-length chain.

### Changed

- Requires `axiolid-evaluate` 0.3.21.
- `section_request` with cut regions no longer refuses a cut whose boundary
  meets itself at a vertex (#373, ADR 0092 amendment): the edges are paired
  there by exact angle and a ring that returns to a vertex is split, so
  corners that touch give two regions and a hole touching its outer ring is a
  hole that touches it. A ring that touches another at every corner is
  placed by the exact midpoint of one of its edges. Still refused by name: an
  edge in and an edge out along one ray or a chain the rings would part
  (`REFUSED_PINCHED_CUT_REGION`), edges that do not alternate or two of one
  direction along one ray (`REFUSED_CROSSING_CUT_RINGS`); a `Both` ring that repeats a
  vertex stays refused.

## [0.3.16] - 2026-11-03

### Added

- The package now ships `LICENSE` and `LICENSE-COMMERCIAL` (symlinks to the
  repository root), so the licence texts travel with the archive (#358).

### Changed

- Requires `axiolid-evaluate` 0.3.20, whose station frame laws bound their
  signed bend (#361).
- The exact point store, the segment arrangement and the ring tracing of the
  section's in-plane faces moved into the new crate
  `axiolid-exact-arrangement`, which the visible faces of
  `axiolid-hidden-lines` use as well (ADR 0094, #383). The section's output
  and refusals are unchanged. A state that cannot occur in that code is now
  `GeomError::Degenerate`, no longer `BackendContractViolation`. Wall clock of
  `cargo bench -p axiolid-reference --bench section_in_plane`, best of nine
  runs on a loaded machine, before and after (ms per section): grid 51.9 and
  54.1, footprint on a grid 51.4 and 51.0, skew grids 44.9 and 46.7, box on a
  slab 0.015 and 0.018. Single runs of one build spread from 53 to 192 ms, so
  this is unchanged within noise.

### Fixed

- `triangulate_simple` and `retriangulate_face` refuse a NaN or infinite
  coordinate with `GeomError::InvalidInput`. `orient2d` now reports such a
  coordinate as uncertain instead of a certain zero (#404), and these
  entry points would otherwise panic on it. `coplanar_overlap` keeps its
  signature and panics on a non-finite coordinate, documented under
  `# Panics`.

## [0.3.15] - 2026-11-03

### Changed

- `ScalarSection::section_request` merges the in-plane faces on an exact
  arrangement instead of refusing overlaps: faces of one side that overlap or
  only share a boundary segment form one region, and where faces of both
  sides coincide the region is `InPlaneSide::Both`. It no longer returns
  `REFUSED_OVERLAPPING_IN_PLANE_FACES` or
  `REFUSED_IN_PLANE_FACES_ON_BOTH_SIDES`. Every `Both` region is a cut region
  when cut regions are requested; a `Both` ring that visits a vertex twice is
  refused as pinched, and a region edge that bounds the cut with in-plane
  faces on the cut's side of it as crossing cut rings. It reports
  `in_plane_crossings`, and `in_plane_pairs` counts segment, piece-triangle
  and hole containment pairs (ADR 0093, #366).
- The request scratch bound is 2560 bytes per element (was 1536), and the
  crossings of in-plane edges are charged against the memory budget as they
  are created.

### Added

- Bench `section_in_plane`: the in-plane faces of a section request in wall
  clock.

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
