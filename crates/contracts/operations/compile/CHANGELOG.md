# Changelog

All notable changes to this crate are documented in this file.

The format follows [Keep a Changelog](https://keepachangelog.com/en/1.1.0/),
and this crate adheres to [Semantic Versioning](https://semver.org/).
Pre-1.0: the minor version is the breaking-change slot, per Cargo's own
caret rule for `0.x` versions.

## [Unreleased]

### Changed

- The minimum supported Rust version is now 1.99 (was 1.88), and the
  workspace builds on Rust 1.99.1 (#349).

## [0.3.4] - 2026-10-14

### Added

- `PolylineCompiler::compile_polyline` (#311): a curve node compiled to a
  3D polyline, a separate operation rather than a second output kind of
  `compile_mesh`, whose signature stays as it is. `PolylineOutcome`
  carries the vertices, whether the curve closes, and a `ChordDeviation`
  reported the way a mesh's deviation is: the requested chord budget and
  a certified upper bound on the distance from every point of the exact
  curve to the polyline, with `meets_requested`. The bound is never
  absent: a curve a compiler cannot bound is refused by name.
- A direct dependency on `axiolid-core` for `Point3` and `Scalar`; it was
  already in this crate's closure through `axiolid-mesh`.

## [0.3.3] - 2026-10-09

### Added

- `MeshClosure::OpenSolid` (#265): the source declared a solid, but the
  compiled mesh is not closed (some edge is not shared by exactly two
  triangles running it in opposite directions). A compiler reports it
  instead of claiming `Solid` for a mesh with boundary edges;
  `CompileOutcome::solid_mesh` refuses it by name. `MeshClosure` is
  `#[non_exhaustive]`, so the new variant is additive.

## [0.3.2] - 2026-09-28

### Changed

- The crates.io page is this crate's own `README.md`, with links to its
  API documentation, its reference page and the source (ADR 0078).

## [0.3.1] - 2026-09-24

### Added

- `MeshClosure` and `CompileOutcome::closure` (#161): whether a compiled mesh
  bounds a solid (`Solid`), is a surface model with area but no volume
  (`Surface`), or was not reported (`Unknown`, the default for `untracked`
  and `tracked`, so existing compilers build unchanged).
  `CompileOutcome::solid_mesh` returns the mesh only for `Solid`, so volume
  readers refuse a surface model instead of measuring a closed shell the
  source never declared a solid. `with_closure` sets it.

## [0.3.0] - 2026-09-23

### Added

- `CompileOutcome` and the provided method `MeshCompiler::compile_mesh_reported` (#115): a compiled mesh with the fate of each attribute channel. The default wraps `compile_mesh` and reports `attribute_fates: None` ("not tracked", not "nothing dropped"), so existing implementations compile unchanged.
