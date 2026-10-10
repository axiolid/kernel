# Changelog

All notable changes to this crate are documented in this file.

The format follows [Keep a Changelog](https://keepachangelog.com/en/1.1.0/),
and this crate adheres to [Semantic Versioning](https://semver.org/).
Pre-1.0: the minor version is the breaking-change slot, per Cargo's own
caret rule for `0.x` versions.

## [Unreleased]

### Added

- An optional `serde` feature, off by default, deriving `Serialize` and
  `Deserialize` for `Line` and `Polyline`: the values the geometry graph's versioned
  wire format carries (#267, ADR 0085). Enums are externally tagged by
  their variant names and unknown fields are refused. No default-build
  change.

## [0.3.1] - 2026-09-28

### Changed

- The crates.io page is this crate's own `README.md`, with links to its
  API documentation, its reference page and the source (ADR 0078).

