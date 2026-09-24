# Changelog

All notable changes to this crate are documented in this file.

The format follows [Keep a Changelog](https://keepachangelog.com/en/1.1.0/),
and this crate adheres to [Semantic Versioning](https://semver.org/).
Pre-1.0: the minor version is the breaking-change slot, per Cargo's own
caret rule for `0.x` versions.

## [Unreleased]

### Fixed

- `half_space::bounded_half_space_in_frame` now places the boundary at the
  authored frame's origin, projected onto the clip plane (#164). It used to
  take only the frame's axes and anchor the boundary at the clip plane's
  origin, so a boundary frame offset within the plane cut the wrong region
  with no error: the mesh stayed closed and correctly wound. The offset along
  the normal is still dropped, so the sweep starts on the clip plane and
  polarity and depth are unchanged. `ReferenceMeshCompiler` passes
  `BoundedHalfSpace.placement.translation` as that origin, so compiled
  bounded half-spaces now honour the placement's translation as well as its
  rotation.

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
