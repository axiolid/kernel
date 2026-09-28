# Changelog

All notable changes to this crate are documented in this file.

The format follows [Keep a Changelog](https://keepachangelog.com/en/1.1.0/),
and this crate adheres to [Semantic Versioning](https://semver.org/).
Pre-1.0: the minor version is the breaking-change slot, per Cargo's own
caret rule for `0.x` versions.

## [Unreleased]

### Changed

- **Breaking:** a candidate triangle index at or beyond the mesh's
  triangle count is refused with the new
  `RayMeshError::TriangleIndexOutOfRange` instead of being skipped by
  `nearest_hit_among`; a broad phase built over a different mesh would
  otherwise report "no hit" for triangles it never tested.
  `triangle_hit` refuses the same index instead of panicking in the
  mesh view.
