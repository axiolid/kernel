# Changelog

All notable changes to this crate are documented in this file.

The format follows [Keep a Changelog](https://keepachangelog.com/en/1.1.0/),
and this crate adheres to [Semantic Versioning](https://semver.org/).
Pre-1.0: the minor version is the breaking-change slot, per Cargo's own
caret rule for `0.x` versions.

## [Unreleased]

### Added

- Exact partial-turn revolution (#172). `revolve_profile_exact` now builds
  a revolution through `0 < |angle| < 2 pi` instead of refusing it: every
  profile the full turn revolves (rectangles, circles, sections, contours
  with arcs and holes, derived, centre-line and composite profiles) becomes
  a closed exact B-rep whose walls are cylinder, cone, plane and torus
  patches trimmed to the swept angle, capped at both ends by the profile
  itself; each hole is a tunnel between the caps. Unlike the full turn, a
  partial turn may touch the axis: a vertex on it sweeps no arc (a cone
  closes at its apex, a planar wall becomes a sector) and a straight
  segment on it sweeps no wall, its edge shared by the two caps. The angle
  follows the right-hand rule about the axis direction, as the mesh path
  does; a section on either side of the axis is accepted. The new
  `revolve_partial::revolve_section_partial` builds one lowered section.
  Volume (`theta R A`) and area (`theta integral(r ds) + 2 A`) match
  Pappus at 90, 180, 270 degrees and small angles
  (`tests/revolve_partial.rs`); `scripts/probe_partial_revolve_mutants.py`
  lists the faults the tests must catch.

### Changed

- `revolve_profile_exact` refuses a turn beyond a full turn by name
  (`"exact revolution beyond a full turn"`) and a zero angle as invalid
  input; both were previously reported as the partial-turn refusal, which
  no longer exists. A section crossing the axis, an arc whose circle
  reaches the axis, and an ellipse remain refused by name.

## [0.3.7] - 2026-09-30

### Fixed

- `boolean_polyhedra_exact` completes long chains of grid-aligned
  subtraction (#199). The depth-2 Menger sponge, 147 differences from a
  unit cube, refused at subtraction 82 with "every probe direction met a
  vertex or edge exactly"; it now finishes closed, with no face enclosing
  zero area, and with exactly the volume of its inputs. Two causes, both
  fixed at the root rather than by dropping fragments:
  - A split point is now the double nearest the exact crossing, not the
    result of `a + (b - a) * t` rounded twice. The old formula cut the
    plane `z = 1/3` at `0.33333333333333326`, so two splits of one exact
    point landed ULPs apart and a later split through the pair emitted a
    ring enclosing no area. Correct rounding depends on the exact point
    alone, and a crossing that is a double comes back exactly.
  - A fragment one ULP wide has no double strictly inside it, so its f64
    centroid rounded onto its own edge and was classified as the wrong
    point, leaving a hole in the shell. Such fragments are now classified
    at an exact dyadic interior point, with every predicate evaluated
    exactly. A fragment that encloses no area at all is refused by name
    ("a split fragment encloses no area") instead of failing later in ray
    classification.
- Coplanar contact on slanted planes. A fragment was classified at its f64
  centroid even when rounding put that point off the fragment's plane, so
  a face lying in the other solid's surface read as inside or outside it:
  a slanted tetrahedron united with itself came back with four faces whose
  windings disagreed, a silently broken solid. The classification point
  must now be certified in the fragment's plane (else the exact probe is
  used), and whether coplanar normals agree is an exact sign rather than an
  f64 dot product.
- The contact-matrix thin-overlap sweep (#200) passes at every overlap down
  to 1e-15: union, intersection and difference are closed and measure the
  exact slab. The boolean was closed all along; the test measured it with
  the 1 um `Tolerance::METRE`, which calls a real 1e-12-wide face
  degenerate, drops it and reports the gap as a hole. It now measures with
  `Tolerance::ZERO`, and the assertions are relative to the slab.

## [0.3.6] - 2026-09-28

### Changed

- The crates.io page is this crate's own `README.md`, with links to its
  API documentation, its reference page and the source (ADR 0078).

## [0.3.5] - 2026-09-28

### Fixed

- Structural sections mesh (#193): `profile_rings` flattens every
  `Profile::Section` family -- I, asymmetric I, L, T, U, C, Z, trapezium --
  from the exact contour `section_contour` builds, fillets and toe radii
  chorded within the budget, straight edges exact. It refused them with
  `Unsupported { ProfileTriangulation }` before. Rectangles with corner
  radii mesh with them, through `rectangle_contour`: the mesh path used to
  drop the radii and mesh a sharp box.
- Circles are chorded from half a step off the axes, not from angle 0
  (#194): the same chords, turned, so a chord's middle, inside the circle,
  sits at every quarter turn. A circular void tangent to a face along an
  axis direction, as openings are, leaves a sliver of material under it
  instead of a chord point on the face, which pinched the solid.

## [0.3.4] - 2026-09-27

### Added

- `bounding::minimum_enclosing_sphere` (#118): the least sphere holding a
  point set, by Welzl's algorithm with exact in/out decisions (diametral,
  least sphere through three points, circumsphere through four), the
  centre enclosed from exact dyadic values and the radius rounded up, with
  `SphereEvidence::error` bounding centre offset and radius excess.
- `bounding::oriented_bounding_box` (#118): a box holding every point,
  exactly (`|(p - centre) . axes[i]| <= half_extents[i]` for the returned
  doubles), never larger than the axis-aligned box. Tries the axis-aligned
  box, the principal axes, and each world axis, principal axis and exact
  hull face normal with the exact minimum-area rectangle across it. It does
  not claim the global minimum volume.
- Depends on `axiolid-exact` for the exact decisions.

## [0.3.3] - 2026-09-27

### Added

- Exact full-turn revolution of circles, hollow circles, hollow rectangles
  and every section with holes (#111). Each hole revolves on its own and
  joins the solid as a void shell. Verified by Pappus on every case.
- Composite profiles are unioned by `profile_lower::composite_regions` over
  one `ArcArrangement` (#111): members may carry arcs and their own
  openings, a member's opening stays open unless another member fills it,
  and members that do not touch become separate solids of one `ExactBRep`,
  for extrusion and revolution alike. Member vertices within the linear
  tolerance are welded first, so members authored to meet do meet.
- `section_lower::circle_contour`: a circle profile (and its bore) as exact
  quarter-arc contours. A circle moved off the origin by a derived profile
  now lowers to that contour instead of being refused.

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

### Changed

- Coaxial booleans whose result encloses a cavity return one solid with a
  void shell (#120) instead of refusing: `boolean_prisms_exact`,
  `boolean_arc_prisms_exact` and their `_solids` variants. A cavity in a
  result of several pieces is still refused by name.

### Fixed

- A sharp rectangle revolved through `revolve_rectangle` wrote its cap
  holes' pcurve intervals forwards although their uses run the edge
  backwards (ADR 0024), so each hole's pcurve ran against its edge: the
  geometric audit put it 8 off the edge and the solid measured 289 instead
  of 188.5. Found by the general boolean (#167); present before #125.

- Exact revolutions are no longer built inside out (#125).
  `revolve_profile_exact` and the contour revolution put their surface
  frames at `(x, y, z) = (X, Z, Y)`, which is left-handed; every loop is
  built anticlockwise in its parameters, so every `Forward` face pointed
  into the solid. The topological and geometric audits compare faces with
  each other and passed it; `exact_properties` measured `-2 pi R A` for
  every Pappus fixture. The frames are now `(X, -Z, Y)`.

## [0.3.2] - 2026-09-25

### Fixed

- `sweep::swept_disk` carries its section frame along the path by
  rotation-minimising frames (double reflection) instead of one fixed axis
  seeded from the first segment (#169). A pipe whose later leg ran along that
  axis was refused ("sweep reference direction must not be parallel to the
  directrix"), and a leg NEARLY along it projected the fixed axis to a
  residue of arbitrary direction, so the ring rotated between stations and
  the tube's volume collapsed with no error: in a real Revit rebar model 878
  bent bars were refused and 6,644 compiled more than 0.5 % short, up to
  64 %. The section is a circle, so the choice of perpendicular only rotates
  it about its own axis. `fixed_reference_sweep` is unchanged: there the
  reference is the author's.

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
