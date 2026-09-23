# Changelog

All notable changes to this crate are documented in this file.

The format follows [Keep a Changelog](https://keepachangelog.com/en/1.1.0/),
and this crate adheres to [Semantic Versioning](https://semver.org/).
Pre-1.0: the minor version is the breaking-change slot, per Cargo's own
caret rule for `0.x` versions.

## [Unreleased]

## [0.3.0] - 2026-09-23

### Added

- Exact extrusion of rounded rectangles (`IfcRoundedRectangleProfileDef`)
  and of hollow rectangles with outer and inner corner radii
  (`IfcRectangleHollowProfileDef`). Each corner is an exact quarter arc that
  extrudes to a cylinder wall. `section_lower::rectangle_contour` builds the
  contour through the same router the structural sections use
  ([#111](https://github.com/axiolid/kernel/issues/111)).
- Exact full-turn revolution of a filled rounded rectangle; each corner
  sweeps a torus.
- Rounded and hollow rectangles are accepted as the basis of a derived
  profile, with the same similarity check every other contour goes through.

### Changed

- Invalid rectangle radii are refused instead of clamped: negative,
  non-finite, wider than the half-extent, an inner radius on a filled
  rectangle, and a hollow section whose corners leave no wall.
