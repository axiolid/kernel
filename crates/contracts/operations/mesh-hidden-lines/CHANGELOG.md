# Changelog

All notable changes to this crate are documented in this file.

The format follows [Keep a Changelog](https://keepachangelog.com/en/1.1.0/),
and this crate adheres to [Semantic Versioning](https://semver.org/).
Pre-1.0: the minor version is the breaking-change slot, per Cargo's own
caret rule for `0.x` versions.

## [Unreleased]

## [0.1.5] - 2026-11-07

### Added

- `Projection` (`Orthographic`, the default, and `Perspective`), the
  `OrthographicView::projection` field, `OrthographicView::perspective(frame,
  near)` (refuses a `near` that is not finite and below zero) and a defaulted
  `supports_projection` on `MeshHiddenLines` and `MeshVisibleFaces` that is
  `true` for `Orthographic` only. Registries route by it, so a provider that
  predates the field is never given a perspective view. The new
  `conformance::perspective::PerspectiveConformanceSuite` infers the claim
  from `supports_projection`: an unclaiming provider must refuse a perspective
  request by name, and a claim fails closed until perspective fixtures exist
  (ADR 0107, #387). No provider draws perspective yet.
- `FaceGroups::Planar` and `HiddenLineMesh::with_planar_groups`: the groups are
  the edge-connected components of triangles joined across an edge of exactly
  two triangles with parallel, same-facing normals, decided exactly, each named
  by its smallest triangle (ADR 0106, #384). Conformance gains the planar
  isometric cube, a flat sheet and a fold.
- `EdgeKind::GroupBoundary` and `EdgeSelection::group_boundary` (default
  `false`) with `EdgeSelection::with_group_boundary`: an edge between two
  triangles of different face groups that is neither a silhouette nor a
  selected crease. With the selection off the groups stay ignored and the
  output is unchanged (ADR 0106, #385). Conformance gains the planar
  isometric cube and a two-group square sheet.

## [0.1.4] - 2026-11-07

### Added

- `VisibleFaceEvidence::snap_displacement`, `rerouted_edges`,
  `collapsed_regions` and `collapsed_holes` (#410, ADR 0100): how far snap
  rounding moved region boundaries, as a `DisplacementBound`, and how many
  edges it rerouted and regions and holes it collapsed. `input_mesh` starts
  them at zero. The conformance suite requires zero in every exact case.

### Changed

- `MeshVisibleFaces` regions are the exact regions snap-rounded to the
  binary64 grid (ADR 0100): a region whose area collapses under rounding is
  dropped and counted, never refused; rings are weakly simple; and
  `VisibleRegion::triangles` lists the triangles frontmost somewhere in an
  exact region it was snapped from. Where no edge is rerouted the vertices
  are still the exact points rounded once.

## [0.1.3] - 2026-11-03

### Added

- The `MeshVisibleFaces` operation contract (ADR 0094, #383): visible face
  regions of meshes in the same orthographic view, with `FaceGroup`s per
  triangle, `FaceSource` (group or near-plane cap), `VisibleFaceRequest`,
  `VisibleFaceLimits`, `VisibleRegion`, evidence and
  `FACES_CAPABILITY_ID`.
- `conformance::faces::FaceConformanceSuite` for `MeshVisibleFaces`
  providers: twelve exact cases (overlapping boxes, a hole under a box, an
  elevation, a cut column on a slab, a sheet that is never capped,
  interpenetrating sheets, a coplanar tie, an edge-on wall, two groups on
  one face, a crop, a far clip and a floating cluster), an area case for
  non-binary64 coordinates, the refusal cases, and determinism and
  orientation checks (#383).
- `HiddenLineMesh::groups` and `HiddenLineMesh::with_groups`; the field is
  read by `MeshVisibleFaces` only.

## [0.1.2] - 2026-11-03

### Added

- The package now ships `LICENSE` and `LICENSE-COMMERCIAL` (symlinks to the
  repository root), so the licence texts travel with the archive (#358).
- `RefusedMesh`: a refused input mesh and its `GeomError`, as listed by the
  registry's audit so a caller learns every refusal in one call (#370).

## [0.1.1] - 2026-10-28

### Added

- Conformance cases for near-plane caps (#138): a solid cut by the near
  plane hides what lies below its cap, a sheet is never capped, a shared
  diagonal needs the perturbation, a plane at the top of a solid cuts
  nothing, a cut solid with a face in the plane is drawn, and a fold under
  an equator edge needs the orientation of the perturbation.

## [0.1.0] - 2026-10-28

### Added

- The `MeshHiddenLines` operation contract (ADR 0090, #138): an
  orthographic view with optional near and far clipping, a 2D crop region,
  boundary, silhouette and crease edge selection, explicit limits,
  visible and hidden segments naming their source edge, evidence counts,
  `HiddenLineError` (a `GeomError` plus the mesh at fault), and a
  conformance suite with self-contained fixtures.
