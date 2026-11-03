# Changelog

All notable changes to this crate are documented in this file.

The format follows [Keep a Changelog](https://keepachangelog.com/en/1.1.0/),
and this crate adheres to [Semantic Versioning](https://semver.org/).
Pre-1.0: the minor version is the breaking-change slot, per Cargo's own
caret rule for `0.x` versions.

## [Unreleased]

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
