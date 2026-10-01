# Changelog

All notable changes to this crate are documented in this file.

The format follows [Keep a Changelog](https://keepachangelog.com/en/1.1.0/),
and this crate adheres to [Semantic Versioning](https://semver.org/).
Pre-1.0: the minor version is the breaking-change slot, per Cargo's own
caret rule for `0.x` versions.

## [Unreleased]

### Added

- `HalfedgeMesh` (#140): an editable halfedge surface mesh with O(1) adjacency, laid out like CGAL's `Surface_mesh` (edge `e` owns halfedges `2e` and `2e + 1`, so `opposite` stores nothing). Built from a `TriMesh` (`from_tri_mesh`) or polygonal faces (`from_faces`) keeping the input numbering, converted back with `to_tri_mesh`, which reproduces the input index buffer exactly.
- Navigation: `next`, `prev`, `opposite`, `source`, `target`, `face`, `edge`, `vertex_halfedge`, `face_halfedge`, `find_halfedge`; counter-clockwise circulators `outgoing_halfedges`, `incoming_halfedges`, `vertex_vertices`, `vertex_faces`; face and hole loops `face_halfedges`, `face_vertices`, `face_faces`, `loop_halfedges`, `boundary_halfedges`, `boundary_loops`. A boundary vertex stores, and circulates from, its boundary halfedge.
- Local edits that refuse by name and keep every invariant and the Euler characteristic: `flip_edge`, `split_edge` (re-triangulating adjacent triangles), `collapse_edge` (link condition, with tetrahedron, lone-triangle and pillow guards), `split_face` (centre fan), `split_face_diagonal`; plus `fill_hole` and `compact`, which renumbers densely and returns a `HalfedgeRemap`.
- `HalfedgeMesh::validate` checks every structural invariant and reports a `HalfedgeInvariantError`.
- `HalfedgeBuildError` names the input a halfedge mesh refuses: `NonManifoldEdge`, `NonManifoldVertex` (a bowtie or two closed fans at a point), `InconsistentOrientation`, `DegenerateFace`, `FaceTooSmall`, `IndexOutOfRange`, `IncompleteTriangle`. `HalfedgeEditError` names refused edits (`LinkCondition`, `WouldDegenerate`, `EdgeExists`, `BoundaryEdge`, `NotATriangle`, ...).
- Typed ids `VertexId`, `HalfedgeId`, `EdgeId`, `FaceId`.

## [0.3.1] - 2026-09-28

### Changed

- The crates.io page is this crate's own `README.md`, with links to its
  API documentation, its reference page and the source (ADR 0078).

## [0.3.0] - 2026-09-23

### Added

- `AttributeFate::then`: the fate of a channel through two sequential steps (dropped wins and keeps the first reason; any interpolation interpolates).
- Corner-indexed attribute channels (#112): `AttributeChannel::corner_indices`, one entry per triangle corner, mirroring `NormalAttribute::indices`. Source formats store texture coordinates this way; positions stay shared, so UV seams no longer force a choice between splitting vertices (breaking closure) and smearing values.
- `AttributeChannel::corner_indexed`, `is_corner_indexed`, `value_count`, `at_corner` (reads either addressing, `None` for an unmapped corner), and `AttributeChannel::UNMAPPED` for triangles that carry no value.
- `validate_structure` checks corner channels: whole tuples, one entry per corner, entries in range, and each triangle fully mapped or fully unmapped. New `MeshValidationError` variants name the channel.
- `DropReason::ConflictingValues`: merged vertices carried different values, so a per-vertex channel could not keep both (#114).
- `DropReason::IncompatibleChannels`: inputs being combined define one channel name with a different width or blend (#115).

### Changed

- **Breaking** (minor slot pre-1.0, ADR 0067): `AttributeChannel` gains the public field `corner_indices`, so struct-literal construction must add `corner_indices: None`; `AttributeChannel::new` is unaffected. `DropReason` is now `#[non_exhaustive]`, so an exhaustive `match` on it needs a wildcard arm. No caller in this workspace or in openbim does either.
