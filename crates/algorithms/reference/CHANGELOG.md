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

- `tessellate_primitive` meshes `Primitive::Torus` and `Primitive::Wedge`
  (#142). The torus is a grid of planar trapezoids sized by the chord
  budget, round the axis for the outer equator and round the tube for the
  tube; horn and spindle tori, and non-positive or non-finite radii, are
  refused by name. The wedge's faces are planar and shared corners of a
  collapsed top are merged, so a ridge or apex wedge is still a closed,
  outward-wound solid; a reversed or non-finite top range is refused.
