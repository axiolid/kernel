# Changelog

All notable changes to this crate are documented in this file.

The format follows [Keep a Changelog](https://keepachangelog.com/en/1.1.0/),
and this crate adheres to [Semantic Versioning](https://semver.org/).
Pre-1.0: the minor version is the breaking-change slot, per Cargo's own
caret rule for `0.x` versions.

## [Unreleased]

### Added

- `spline_pair_intersection`: every component of the section of two
  B-spline surfaces within parameter windows, as `PairSection3` curves.
  Bezier sub-patch pairs are split until their normal cones are apart, so
  no closed loop hides in one; sub-patch edge crossings seed every
  component, and each is followed on both surfaces to the window's edge,
  where its last node is solved exactly. `exact_surface_intersection`
  returns it for two B-splines (`Derivation::PairTrace`).
- `exact_curve_surface_intersection` for a B-spline curve against a
  B-spline surface: hull pruning over the Bezier pieces of both, then
  Newton.

### Changed

- The trace separates the singular points it meets. An isolated point
  where the surfaces touch without crossing is dropped, with the boxes
  about it that rounding cannot decide. A trace that finds only such
  points answers `NotRegularCurve`, as the closed forms do for touching.
  Branches crossing at a saddle are still refused by name: within about
  `sqrt(rounding / curvature)` of the crossing (micrometres at metre
  scale), where the branches run cannot be decided in doubles. Whole-turn
  windows start an irrational fraction of a radian past `-pi`, so a
  symmetric section's special points never sit on the window's edge.

### Added

- `section_curve_curve_intersection2` (B8, #119): lines, conics,
  sinusoids, quadratic and angle graphs and implicit curves meet each
  other, in any pairing. One piece is traced over its span and the other's
  defining field has its roots isolated along it; roots on a graph's other
  branch are filtered out.
- `section_curve_curve_intersection3` (B8): space curves of the section
  families meet where one crosses a surface the other lies on, kept within
  a tolerance.
- `extrema::minimum_distance` (B7, #119): the smallest distance between
  points, curve spans and surface patches of any analytic or B-spline
  family and traced sections, bracketed and certified, with witness points.
  It is a branch and bound whose lower bounds are exact image boxes and
  projections on the joining direction, in mean-value form so they close
  quadratically at a closest point.

- Sections of B-spline surfaces by planes, quadrics and tori (#119,
  ADR 0077). The analytic equation is read on the spline's rational Bezier
  patches as Bernstein polynomials (knot insertion, then Bernstein
  products) and traced with the same certified subdivision, clipped to the
  spline's domain. `implicit_surface_intersection` and
  `trace_section_pcurves` take a B-spline carrier.
- `exact_curve_surface_intersection` takes a B-spline surface for lines,
  circles and ellipses. The curve is written as the meeting of two
  analytic surfaces; the first is traced on the spline and the second's
  roots are found along it.

- A plane through a cone's apex cuts rays along the rulings of the
  modelled nappe (`Derivation::ConeApexRulings`, spans `[0, inf)`). It was
  refused before. A plane flatter than the cone meets it only at the apex
  (`NotRegularCurve`).
- Curve/surface intersection for the section families (B9, #119):
  `exact_curve_surface_intersection` takes an `ImplicitSection`, and
  `section_curve_surface_intersection` takes any ruled, torus or traced
  section over a span. Roots are isolated along the curve's cells with
  interval bounds. Parameters are `ExactCurveParameter::Certified`.
- `trace_section_pcurves`, `extract_stretch` and `implicit_view` give a
  section's implicit pcurve on any analytic face, cut out between given
  points, and a ruled or torus section as an implicit curve.

- Traced sections (#119, ADR 0077): `exact_surface_intersection` now
  builds a torus against a cylinder, elliptical cylinder, cone or torus off
  its axis (`Derivation::ImplicitTrace`). The section is found in the
  torus's parameters as every component of the other surface's equation,
  by certified subdivision into monotone cells. There is no marching, so no
  loop is missed for want of a small step. A singular point (surfaces
  touching where branches cross) is refused as `NotRegularCurve`.
  `implicit_surface_intersection` exposes it with an explicit window for
  any analytic pair, and `section_field_of` the field itself.

- `exact_surface_intersection` derives ruled quadric sections (#119,
  ADR 0076): a cylinder or elliptical cylinder against a plane, sphere,
  cylinder, elliptical cylinder or cone, and a cone against a plane, sphere
  or cone (on the modelled nappes, decided exactly), when no
  line, circle or ellipse applies -- pipe tees, off-axis sphere/cylinder
  junctions, oblique cone cuts, parabolas and hyperbolas. Which spans of
  angle carry the curve is decided by exact root isolation of the
  discriminant. `ExactIntersectionCurve` gains `spans` and `Derivation`
  gains `RuledQuadricSection`.
- `exact_surface_intersection` derives a ring torus's section by a plane or
  sphere off its axis (#119, ADR 0076) -- a pipe bend meeting a wall, a ball
  against a ring -- as `u(v)` solving `A(v) cos u + B(v) sin u = C(v)`,
  with spans and wrap points decided exactly
  (`Derivation::TorusAngleSection`).

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
