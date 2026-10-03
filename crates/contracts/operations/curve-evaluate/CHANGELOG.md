# Changelog

All notable changes to this crate are documented in this file.

The format follows [Keep a Changelog](https://keepachangelog.com/en/1.1.0/),
and this crate adheres to [Semantic Versioning](https://semver.org/).
Pre-1.0: the minor version is the breaking-change slot, per Cargo's own
caret rule for `0.x` versions.

## [Unreleased]

## [0.3.2] - 2026-10-03

### Changed

- `CurveEvaluator::frame_at` documents that a curve carrying its own roll
  (`Curve3::Banked`, #240) is framed by it: its section frame, rolled about
  the tangent under the curve's bank convention (ADR 0081).

## [0.3.1] - 2026-09-28

### Changed

- The crates.io page is this crate's own `README.md`, with links to its
  API documentation, its reference page and the source (ADR 0078).

