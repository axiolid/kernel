# Changelog

All notable changes to this crate are documented in this file.

The format follows [Keep a Changelog](https://keepachangelog.com/en/1.1.0/),
and this crate adheres to [Semantic Versioning](https://semver.org/).
Pre-1.0: the minor version is the breaking-change slot, per Cargo's own
caret rule for `0.x` versions.

## [Unreleased]

## [0.3.7] - 2026-10-04

### Added

- `elevated`: `axiolid-evaluate`'s derivatives, certified chord and
  derivative bounds and grade corners of elevated and banked curves
  (#252), re-exported on the same footing as `curve` and `surface`.

## [0.3.6] - 2026-10-03

### Added

- `station`: `axiolid-evaluate`'s station resolution (#241), re-exported
  on the same footing as `curve` and `surface`.

## [0.3.5] - 2026-10-03

### Added

- Re-exports `axiolid_evaluate::arc_parameter` and
  `axiolid_evaluate::chain` (#239).

## [0.3.4] - 2026-10-03

### Added

- `bound`, re-exported from `axiolid-evaluate` like `curve` and `surface`
  (#232).

### Fixed

- `tessellate_primitive` refuses a cylinder or a cone whose chord budget
  needs more than 4096 segments with `BudgetExceeded`, as it already did
  for spheres and tori, instead of clamping silently to a coarser mesh
  (#232). Their bound is now documented: both are curved one way only, so
  the ring's sagitta is the whole distance from any point of the exact
  surface, caps included, to the mesh, and it gets the whole budget.

## [0.3.3] - 2026-10-02

### Fixed

- `tessellate_primitive` keeps every point of a sphere and a torus within
  the chord budget of the mesh, not only its vertices (#231). Both
  directions of each used to get the whole budget, so their sagittas added
  inside a triangle (a 0.3 m sphere at 1 mm lay 1.99 mm from its mesh, a
  torus 1.90 mm); each now gets half, and the sphere's stack count is
  rounded up rather than down, which had let an odd segment count leave
  the polar step alone above its share. A budget that needs more than
  4096 segments is refused with `BudgetExceeded` instead of clamped.

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
