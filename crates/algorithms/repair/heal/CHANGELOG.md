# Changelog

All notable changes to this crate are documented in this file.

The format follows [Keep a Changelog](https://keepachangelog.com/en/1.1.0/),
and this crate adheres to [Semantic Versioning](https://semver.org/).
Pre-1.0: the minor version is the breaking-change slot, per Cargo's own
caret rule for `0.x` versions.

## [Unreleased]

### Added

- Repairs carry corner-indexed channels: weld leaves them untouched (they index values, not positions, so a seam is lossless), dropping and flipping triangles move their entries (#112).

### Fixed

- Repairs keep attribute channels and normals in step with the geometry they rewrite (#114). Weld compacts per-vertex channels and normals; a seam drops the channel by name, a hard edge switches normals to corner-indexed. Dropping or flipping triangles moves corner-indexed normals with them.
- `RepairReport::attribute_fates` names every input channel's fate.
