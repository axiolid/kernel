# Changelog

All notable changes to this crate are documented in this file.

The format follows [Keep a Changelog](https://keepachangelog.com/en/1.1.0/),
and this crate adheres to [Semantic Versioning](https://semver.org/).
Pre-1.0: the minor version is the breaking-change slot, per Cargo's own
caret rule for `0.x` versions.

## [Unreleased]

### Added

- `exact_surface_intersection` derives ruled quadric sections (#119,
  ADR 0076): a cylinder or elliptical cylinder against a plane, sphere,
  cylinder or elliptical cylinder, and a cone against a plane, when no
  line, circle or ellipse applies -- pipe tees, off-axis sphere/cylinder
  junctions, oblique cone cuts, parabolas and hyperbolas. Which spans of
  angle carry the curve is decided by exact root isolation of the
  discriminant. `ExactIntersectionCurve` gains `spans` and `Derivation`
  gains `RuledQuadricSection`.

- Exact intersection of analytic curves (#119): `exact_curve_surface_intersection`
  (line, circle or ellipse against plane, cylinder, elliptical cylinder, cone,
  sphere or torus) and `exact_curve_curve_intersection2`/`3` (lines, circles,
  ellipses). The equations become integer polynomials in the line parameter or
  the conic's half-angle parameter, solved with `axiolid-exact`: each hit, its
  multiplicity (2 = touching) and "the curve lies on the surface" are exact
  decisions, and a hit at a conic's half-angle singularity (`theta = pi`) is
  reported as `Antipode`. A cone counts only its modelled nappe; a line on the
  cone through its apex is refused as `PartialOverlap`. Points are rounded to
  `f64` for output only. B-spline operands stay on the certified tier.

### Fixed

- `exact_surface_intersection` decides tangent, parallel and perpendicular
  cases exactly (#119). In `f64` a plane exactly tangent to a sphere along
  a normal like (3, 2, 6) came out as a circle of radius about 1e-7, an
  exactly perpendicular oblique plane cut a cylinder in a near-circular
  ellipse, and an exactly parallel one produced a huge ellipse instead of
  rulings. Only output coordinates are rounded now.
