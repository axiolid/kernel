# Changelog

All notable changes to this crate are documented in this file.

The format follows [Keep a Changelog](https://keepachangelog.com/en/1.1.0/),
and this crate adheres to [Semantic Versioning](https://semver.org/).
Pre-1.0: the minor version is the breaking-change slot, per Cargo's own
caret rule for `0.x` versions.

## [Unreleased]

### Added

- Corner-indexed attribute channels (#112): `AttributeChannel::corner_indices`, one entry per triangle corner, mirroring `NormalAttribute::indices`. Source formats store texture coordinates this way; positions stay shared, so UV seams no longer force a choice between splitting vertices (breaking closure) and smearing values.
- `AttributeChannel::corner_indexed`, `is_corner_indexed`, `value_count`, `at_corner` (reads either addressing, `None` for an unmapped corner), and `AttributeChannel::UNMAPPED` for triangles that carry no value.
- `validate_structure` checks corner channels: whole tuples, one entry per corner, entries in range, and each triangle fully mapped or fully unmapped. New `MeshValidationError` variants name the channel.
- `DropReason::ConflictingValues`: merged vertices carried different values, so a per-vertex channel could not keep both (#114).

### Changed

- **Breaking** (minor slot pre-1.0, ADR 0067): `AttributeChannel` gains the public field `corner_indices`, so struct-literal construction must add `corner_indices: None`; `AttributeChannel::new` is unaffected. `DropReason` is now `#[non_exhaustive]`, so an exhaustive `match` on it needs a wildcard arm. No caller in this workspace or in openbim does either.
