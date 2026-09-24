# Changelog

All notable changes to this crate are documented in this file.

The format follows [Keep a Changelog](https://keepachangelog.com/en/1.1.0/),
and this crate adheres to [Semantic Versioning](https://semver.org/).
Pre-1.0: the minor version is the breaking-change slot, per Cargo's own
caret rule for `0.x` versions.

## [Unreleased]

## [0.3.1] - 2026-09-24

### Added

- `PolygonMesh` faces that are not plain triangles compile (#160): n-gons,
  concave faces and faces with holes (IFC4 `IfcIndexedPolygonalFaceWithVoids`)
  are triangulated in their own plane, keeping the authored positions and
  winding. Plain triangles keep their exact corner order as before. A face
  whose corners leave its plane by more than the linear tolerance, that has
  no area, or whose rings cross is refused with an error naming its index.
- B-reps with shells but no solid tessellate (#161): every shell is
  tessellated as authored and the result is reported as
  `MeshClosure::Surface` through `compile_mesh_reported`, even when the
  shell is closed. Collections are `Solid` only if every member is, and a
  boolean with a surface operand is refused. Authored meshes report
  `Solid` exactly when they are closed, consistently wound two-manifolds.

### Changed

- A `PolygonMesh` with non-triangular faces used to fail with
  `GeomError::Unsupported`; it now compiles. A B-rep with no solid and no
  shell is refused as "neither a solid nor a shell" instead of "has no
  solid".

## [0.3.0] - 2026-09-23

### Fixed

- `Instance` and `Collection` nodes no longer drop attribute channels and normals (#115). Both used to rebuild the mesh from positions and indices only, so a textured item lost its `uv` channel as soon as a product had a second item or was instanced — silently.
  - `Instance` carries channels through unchanged and now transforms normals by the inverse transpose instead of dropping them. Under a mirroring transform, corner-indexed channels and normals swap corners with the triangle.
  - `Collection` merges channels by name. A channel on only some members becomes corner-indexed, with the other members' triangles `UNMAPPED`; a channel on every member as per-vertex stays per-vertex. Members defining one name with a different width or blend drop it as `DropReason::IncompatibleChannels`.
  - Booleans keep the provider's channel fates, composed onto what each operand already went through, instead of discarding the evidence.

### Added

- `ReferenceMeshCompiler` implements `MeshCompiler::compile_mesh_reported`: the compiled mesh plus each channel's fate on its way to the root (worst over parallel members, sequential through booleans).
