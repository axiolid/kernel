# Changelog

All notable changes to this crate are documented in this file.

The format follows [Keep a Changelog](https://keepachangelog.com/en/1.1.0/),
and this crate adheres to [Semantic Versioning](https://semver.org/).
Pre-1.0: the minor version is the breaking-change slot, per Cargo's own
caret rule for `0.x` versions.

## [Unreleased]

### Changed

- The licence is now `MPL-2.0 OR LicenseRef-Commercial` (was `MPL-2.0`):
  a commercial licence is available from the copyright holder
  (`LICENSE-COMMERCIAL`, ADR 0089, #358).
- The minimum supported Rust version is now 1.99 (was 1.88), and the
  workspace builds on Rust 1.99.1 (#349).

## [0.1.0] - 2026-10-02

### Added

- `tangent_circles` (#159, ledger row B14): every circle tangent to three
  given points, lines or circles, in any combination (Apollonius' problem
  and its degenerate cases), returning every real solution. Each of the
  three constraints may independently require internal or external
  tangency for a circle constraint, or either side for a line constraint,
  enumerated as sign choices; solutions are deduplicated and sorted for a
  deterministic order. Refused by name: non-finite input, a non-positive
  circle radius, a zero line direction, and a configuration whose
  elimination system is singular for every sign choice (three coincident
  or otherwise inseparable constraints).
