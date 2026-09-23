# Changelog

All notable changes to this crate are documented in this file.

The format follows [Keep a Changelog](https://keepachangelog.com/en/1.1.0/),
and this crate adheres to [Semantic Versioning](https://semver.org/).
Pre-1.0: the minor version is the breaking-change slot, per Cargo's own
caret rule for `0.x` versions.

## [Unreleased]

### Added

- `merge_fates`: compose per-channel fates across sequential steps (#116).

### Fixed

- Composed evidence reported only the last step's attribute fates: `BooleanEvidence::absorb` (the `subtract_many`/`union_many` defaults) and `symmetric_difference_via_composition` now compose them, so a channel a middle step derived or dropped is reported that way.
