# Changelog

All notable changes to this crate are documented in this file.

The format follows [Keep a Changelog](https://keepachangelog.com/en/1.1.0/),
and this crate adheres to [Semantic Versioning](https://semver.org/).
Pre-1.0: the minor version is the breaking-change slot, per Cargo's own
caret rule for `0.x` versions.

## [Unreleased]

## [0.3.2] - 2026-09-28

### Changed

- The crates.io page is this crate's own `README.md`, with links to its
  API documentation, its reference page and the source (ADR 0078).

## [0.3.1] - 2026-09-27

### Added

- `Primitive::Torus` (#142): a ring torus about local +z, by major and
  minor radius. Horn and spindle tori are not solids and are refused by
  the tessellator.
- `Primitive::Wedge` (#142): OCCT's `MakeWedge` general form with the
  height along local +z -- a base rectangle at z = 0 and a top rectangle,
  narrowed or shifted, at z = height. The top may collapse to a ridge or
  an apex.
