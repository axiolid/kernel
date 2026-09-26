# Changelog

All notable changes to this crate are documented in this file.

The format follows [Keep a Changelog](https://keepachangelog.com/en/1.1.0/),
and this crate adheres to [Semantic Versioning](https://semver.org/).
Pre-1.0: the minor version is the breaking-change slot, per Cargo's own
caret rule for `0.x` versions.

## [Unreleased]

### Added

- `FaceDomain` (#167): an exact face's parameter domain, answering whether
  a point lies in the face with a certificate or not at all; a point outside
  the domain's box in a coordinate that does not wrap is decided at once.
- `boundary_distance` and `boundary_clearance` (#125, ADR 0074): an
  interval certain to contain the distance between the boundaries of two
  exact B-reps, with witness points on both, and a comparison with a limit
  that answers `Clearance::Below`, `Above` or `Indeterminate`. Branch and
  bound over face patches and edge spans; lower bounds from Lipschitz
  spheres and exact projection ranges, upper bounds only from points
  certified on the boundary.

- `exact_properties` measures curved faces (#125, ADR 0073): cylinders,
  cones, spheres, tori, elliptical cylinders, B-spline faces, and planar
  faces bounded by arcs or ellipses. Each face is integrated over its own
  parameter domain by Green's theorem round its pcurves, with adaptive
  Gauss-Kronrod quadrature held to a relative error of 1e-13; nothing is
  tessellated. Seams, poles and apexes, and torus faces bounded by meridians
  are handled; a boundary that encloses nothing in the surface's parameters
  is refused.

### Changed

- **Breaking:** `ExactMeasureError` is `#[non_exhaustive]` and gains
  `ParameterDomain`, `Evaluation` and `NotConverged`. `NonPlanarFace` now
  means a surface family the module cannot integrate at all.
- The `exact` feature now also enables `axiolid-evaluate` and
  `axiolid-curve`.

### Fixed

- `FaceDomain` and the certified distance read a point on a periodic face
  with a negative angle as outside it (#167): the whole-period shifts tried
  had the wrong sign, so `-0.2` was tried at `-0.2 - 2 pi`, not at
  `2 pi - 0.2`. Faces reaching a pole hid it, since the pole adds its own
  crossing.

- A planar face with a hole reported the hole's area added to its own: the
  fan summed triangle magnitudes. Areas are now summed as vectors, so a hole
  subtracts (a 4 x 4 plate with a 2 x 2 hole read 20 per cap, not 12).

- `exact_properties` honours face, shell-use and bound orientation. It
  read loop winding alone, which is only right for faces used forward; a
  `Reversed` cap off the plane `z = 0` added its volume instead of
  subtracting it (a unit cube at `2 <= z <= 3` measured 7/3). Every
  solid tested before sat on `z = 0`, where the error vanishes.
