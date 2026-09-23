# Changelog

All notable changes to this crate are documented in this file.

The format follows [Keep a Changelog](https://keepachangelog.com/en/1.1.0/),
and this crate adheres to [Semantic Versioning](https://semver.org/).
Pre-1.0: the minor version is the breaking-change slot, per Cargo's own
caret rule for `0.x` versions.

## [Unreleased]

### Added

- `DropReason::ConflictingValues`: merged vertices carried different values, so a per-vertex channel could not keep both (#114).
- `DropReason` is now `#[non_exhaustive]`, so future reasons are not breaking.
