# Changelog

All notable changes to this crate are documented in this file.

The format follows [Keep a Changelog](https://keepachangelog.com/en/1.1.0/),
and this crate adheres to [Semantic Versioning](https://semver.org/).
Pre-1.0: the minor version is the breaking-change slot, per Cargo's own
caret rule for `0.x` versions.

## [Unreleased]

### Added

- `clip_arc_prism_exact` (#120): an arc prism cut by a half-space whose
  plane passes between its caps, the "column under a sloped roof" case.
  Cylindrical walls stay `Cylinder` faces trimmed by an exact `Ellipse3`
  edge with a `Sinusoid2` pcurve (ADR 0071); planar walls get sloped edges;
  the cut cap is unnamed. A plane parallel to the axis is refused by name.

- `boolean_prisms_exact_solids` and `boolean_arc_prisms_exact_solids`
  (#120): coaxial booleans whose result falls apart into separate pieces
  return one solid per piece, ordered by lowest vertex (x, then y), each
  audited on its own. An empty result is an empty list. The single-solid
  functions keep refusing a disconnected result, so callers that expect
  one solid are not silently handed the first piece.

### Changed

- `boolean_arc_prisms_exact` runs on the exact arc overlay (ADR 0070) and
  builds results it used to refuse: a result with interior holes becomes a
  solid with through-passages (#120), and a result starting above `z = 0`
  is extruded from its own base height. Disconnected results go through
  the `_solids` variants.
- Stepped coaxial booleans are built, not refused (#120, ADR 0072):
  `boolean_prisms_exact`, `boolean_arc_prisms_exact` and their `_solids`
  variants return a union of prisms with different spans, a difference
  whose tool stops inside the subject (notch, counterbore, blind pocket,
  slot through the middle heights) as exact solids with their ledge faces.
  Walls are named after the operand edge they lie on, caps and ledges
  after the operand cap that made them. A result enclosing a cavity, and
  pieces touching only along an edge, are refused by name.
- `clip_arc_prism_exact` builds a plane that crosses a cap inside the
  section: the part of the old cap that survives keeps its name, next to
  the unnamed cut.
- `boolean_stepped` docs: the bands are the lighter alternative to the
  stepped solid; their volumes are checked against it.

## [0.3.1] - 2026-09-24

### Fixed

- `extrude` (and so `extrude_profile` and the reference mesh compiler) wound
  a solid inside-out when the extrusion direction pointed below the profile
  plane (`direction.z < 0`), e.g. an opening body extruded downward from its
  lintel (#166). The signed volume was `-area * depth`, so any boolean using
  the solid refused it as inside-out. Such a solid is now outward-oriented
  with the same magnitude, for outer and hole loops alike.
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
