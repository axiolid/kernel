# Changelog

All notable changes to this crate are documented in this file.

The format follows [Keep a Changelog](https://keepachangelog.com/en/1.1.0/),
and this crate adheres to [Semantic Versioning](https://semver.org/).
Pre-1.0: the minor version is the breaking-change slot, per Cargo's own
caret rule for `0.x` versions.

## [Unreleased]

### Added

- `profile::ring_touches` and `profile::RingTouch` (#265): every vertex
  lying inside another ring edge, which `triangulate_with` under
  `PinchPolicy::Accept` inserts into that edge, found by the same exact
  validation and listed once per vertex and edge. A caller whose rings
  share their edges with neighbouring patches (the faces of a closed
  shell) splits the neighbour's copy of each such edge at the same vertex,
  so the patches do not meet at a T-junction.

### Fixed

- **Exact extrusions whose direction points against the profile normal
  build (#275).** `extrude_profile_exact` refused any offset with
  `offset.z <= tolerance` as `"non-forward planar extrusion"`, so an
  opening cut down from a slab's top (`ExtrudedDirection (0, 0, -1)`) had
  no exact result and a boolean with it no certified mesh deviation. Such
  an extrusion is the forward prism along `(o.x, o.y, -o.z)` mirrored in
  the profile plane; it is now built that way, through
  `ExactBRep::transformed`, for every profile family and for oblique
  directions. Negating `z` is exact, so the profile-plane cap keeps the
  profile's coordinates bit for bit, and the reflection flips every face,
  so the solid stays outward with the forward prism's positive volume.
  Only a direction within tolerance of the profile plane
  (`|o.z| <= tolerance`) is refused, now named `"extrusion direction in
  the profile plane"`. The mesh path already built these solids; both now
  agree. Mutation probe: `scripts/probe_downward_extrusion_mutants.py`.

## [0.3.15] - 2026-10-04

### Added

- **Rings that touch at single points triangulate on 2D and surface
  paths (#262).** The certified clipper of #253 refuses rings touching at
  one vertex, which is right for a solid (the extrusion shares one wall
  edge between four faces) but refused valid regions a consumer re-feeds
  from overlay output: two rooms meeting at a corner, an L-shaped room
  wrapped round a column corner, a corridor eroded to a point, the same
  pinch reported through even-odd fill, a room minus door zones touching
  its walls. `profile::triangulate_with(rings, PinchPolicy)` makes the
  choice explicit: `PinchPolicy::Refuse` is `triangulate` (extrusion and
  loft caps, unchanged), `PinchPolicy::Accept` triangulates the region
  (ADR 0083 amendment). Two edges may then meet at one point that is a
  vertex of at least one of them; the vertex is inserted into the other
  edge, coincident vertices become one (triangles use a point's first
  index in `outer ++ holes`), each ring is split at its repeated vertices
  into simple loops turned by their nesting depth (inside an odd number of
  loops is inside), the loops' edges must alternate leaving and arriving
  round each shared vertex (else the rings cross there, refused by name),
  and each wedge there becomes one node, so the region falls into its
  connected parts -- two lobes of a pinch separate, a hole touching the
  outer ring merges with it, holes touching in a ring close off a part of
  their own -- each bridged and clipped as before. The certificate is the
  same over the split loops, with each visit of a pinch point counted as
  its own vertex and `n + 2h - 2c` triangles for `c` parts. Rings touching
  nowhere triangulate identically under both policies; crossing and
  overlapping rings stay refused under both. Tests
  (`tests/profile_pinches.rs`): the consumer's cases (rooms at a corner,
  an L round a column, an eroded corridor, even-odd figure eights, a room
  minus door zones), footprints touching at a corner, holes touching the
  outer ring at a corner, inside an edge, at every corner and across the
  region, holes touching each other at a corner, inside a side, at their
  extreme corner and in a ring that closes off a part with its own hole,
  and an outer ring pinched round a hole, each from every start vertex and
  either way round, tile exactly (strictly counter-clockwise, exact area,
  no vertex on a triangle, untwinned edges covering every ring edge once)
  and are refused by `triangulate` and `extrude_profile` by name; 300
  random grid layouts of checkerboard squares and touching diamonds tile
  exactly, 200 off the grid tile, and layouts shrunk apart triangulate
  identically under both policies; crossing at shared vertices and
  overlapping rings are refused by name. Mutation probe:
  `scripts/probe_ring_triangulation_mutants.py` (52/52 killed).

### Fixed

- **Profiles with several holes extrude to closed solids (#253).**
  `profile::triangulate` delegated to `earcut`, which drops nodes where
  the bridged ring runs straight on and lets only reflex nodes block an
  ear. For two holes side by side in one horizontal band (a 4 x 4 outer
  ring, 1 x 1 holes at x in [-1.5, -0.5] and [0.5, 1.5]) it returned 12
  triangles of the right area instead of 14: one edge ran along the
  band's bottom line past both holes' inner corners, a T-junction, so the
  plain extrusion had 8 boundary edges and failed its volume and closure
  checks. Triangulation is now the crate's own ear clipper, every
  decision an exact `orient2d` sign (ADR 0083, superseding ADR 0015 for
  this crate): holes are bridged in order of decreasing largest x from
  their rightmost vertex, each bridge accepted only when it enters the
  polygon at both ends and touches no ring edge or earlier bridge, a
  vertex on it included; an ear is refused if any other vertex lies in
  its closed triangle. The output is certified before it is returned --
  every triangle strictly counter-clockwise, every ring edge used once
  from inside, every other edge once in each direction -- and refused
  with `Degenerate` otherwise. Rings may come either way round; every
  vertex is now a triangle corner, `outer ++ holes` unchanged. Ring sets
  that bound no polygon with holes are refused with `InvalidInput`
  naming the ring: a non-finite or repeated vertex, a ring folding back
  on or crossing itself, holes overlapping or touching each other or the
  outer ring (including at one vertex, which earcut accepted and which
  extrudes to a non-manifold edge), a hole outside the outer ring or
  inside another hole. The crate no longer depends on `earcut`. Tests
  (`tests/profile_holes.rs`): holes side by side, stacked, three on one
  ray, a hole vertex exactly on another's ray, five holes 1/1024 from the
  outer ring, and three flattened round holes in a row each triangulate
  with n + 2h - 2 counter-clockwise triangles, every ring edge once and
  every other edge twinned, cover exactly outer minus holes, and extrude
  to a closed, outward solid of volume area x depth; so do a triangle
  and an L with a straight or reflex vertex, from every starting vertex,
  as outer ring and as hole, a vertex in line with another hole's side
  just past its end, and 400 seeded
  random layouts of rectangles, triangles, L shapes, diamonds and
  rectangles with collinear midpoints on a quarter grid (compared with
  `==`) and 200 off-grid layouts of regular 3- to 40-gons; each invalid
  ring set above is refused by name; the certificate refuses earcut's
  output for the issue's profile. Mutation probe:
  `scripts/probe_ring_triangulation_mutants.py` (31/31 killed).

## [0.3.14] - 2026-10-04

### Fixed

- **Parametric profiles with decimal sizes lower at `Tolerance::ZERO`
  (#250).** `contour_lower::contour_to_arc_ring` demanded bit-equal
  segment joints at `Tolerance::ZERO`, which a `Line2`/`Circle2` contour
  cannot give for non-dyadic sizes: a line ends at a rounded
  `origin + direction`, an arc at the `cos`/`sin` of its sweep. So every
  I section with IPE or HEA sizes (sharp or with root fillets), and the
  other families with decimal sizes, was refused before any boolean ran
  ("contour segments leave a gap of 2.6e-18"). The section router already
  computes each corner and tangent point once and hands it to both
  segments; the ring takes one vertex per joint (the leaving segment's
  start, a line's stored origin bit for bit), so it is closed by
  construction. The joint check now allows the larger of the tolerance
  and the rounding of the two evaluations meeting there (eight machine
  epsilons of the magnitudes they are computed from); a contour open by
  more is still refused at `Tolerance::ZERO`. Rings are bit-identical to
  before wherever lowering succeeded, so results at a positive tolerance
  are unchanged. Tests: I (IPE 300/200, HEA 200, HEB 340; sharp, root
  fillets, toe radii), asymmetric I, T, U, L, Z, C, trapezium, rounded
  and hollow rectangles and an annulus lower at ZERO with each line's
  shared corner as its ring vertex, and extrude at ZERO to their
  closed-form areas, fillet terms `(1 - pi/4) r^2` included; tapered
  I/U/L/T lower and extrude at ZERO; joints far from the origin close to
  their own rounding; contours open by 1e-13 or 1e-9 are refused at ZERO.
  In `axiolid-mesh-compile`, an IPE 300 beam (with and without fillets)
  minus a round web hole clear of the fillets compiles at ZERO under
  exact placements with an empty report and the closed-form volume.
  Mutation probe: `scripts/probe_section_zero_tolerance_mutants.py`
  (9/9 killed).

## [0.3.13] - 2026-10-03

### Fixed

- **Restored: swept disks mitre their corners again (#245).** Since 0.3.11
  (#232) `pipe::swept_disk_along_pieces` refused a non-tangent joint as
  "undefined"; that was a regression (0.3.11 and 0.3.12). A disk swept
  round a corner between two straight segments is mitred at half angle,
  as `IfcSweptDiskSolid` defines it: both tubes are cut by the plane that
  bisects the two tangents, and the two pieces share one ring on it whose
  every vertex lies on both exact cylinders (the outgoing frame is the
  least rotation of the incoming one, which agrees with the reflection in
  the mitre plane). The walls are planar trapezoids, exact faces of the
  ring prism; the ring is chorded to `(c/2) cos(theta/2)` for the
  sharpest mitre, because the mitre stretches the section by
  `1 / cos(theta/2)`, so the exact tube stays within the chord budget
  (derivation in the module notes). Hollow disks mitre their bore the same
  way. The near-tangent tolerance `joint_tolerance` is unchanged. Refused
  by name: a segment whose mitres reach past its length
  (`L <= r |g_perp|`, `r tan(theta/2)` for one mitre: the tube would cut
  through itself), a reversal, and a corner beside an arc (no ring lies
  on both the cylinder's and the torus's cut).
- A disk radius equal to the fillet or bend radius is refused with its own
  reason (#245): the format rule permits a fillet radius equal to the disk
  radius, but the bend is then a horn torus whose inner wall pinches to a
  point on the bend's axis, which no closed two-manifold mesh bounds
  without meeting itself there.

## [0.3.12] - 2026-10-03

### Changed

- A `Curve2::Chain` profile segment is named in deviation reports
  (`arc-length chain profile segment with an unbounded piece`) and in
  contour-lowering refusals (#239); a chain the flattener certifies is
  bounded by the chord budget like any other certified family.

## [0.3.11] - 2026-10-03

### Added

- `pipe::swept_disk_along_pieces`, with `pipe::PathPiece` and
  `pipe::joint_tolerance` (#232): a disk, optionally hollow, swept along a
  chain of straight segments and circular arcs, with every point of the
  exact tube within the chord budget of the mesh. Half the budget chords
  the disk; each piece is then bounded on its own by the #231 span bound
  (exact along an arc, where the stations are rotations about its axis,
  and zero along a segment), its end frames exact. Consecutive pieces
  share one station, so the tube is watertight and wound one way; the
  shared station's measured distance from the outgoing piece's own start
  station is added to that piece's first span. A joint turning by more
  than `2 asin(c / (8 r))` is a corner and is refused by name, as is a
  gap above a quarter of the budget, a bend radius at or below the disk
  radius, and a bend needing more than 4096 steps (`BudgetExceeded`). A
  fillet radius rounds each corner between two segments with a tangent
  arc (`IfcSweptDiskSolidPolygonal`); a fillet that overruns its
  segments, a disk at or above the fillet radius, a reversal and a corner
  beside an arc are refused by name. The derivation is in the module
  notes.

- `profile::profile_deviation` and `ProfileDeviation` (#232): how far a
  profile's exact boundary may lie from the rings `profile_rings` flattens
  it to. The chord budget for every segment family the flattener certifies,
  plus the largest merge of near-duplicate points, scaled by a derived
  profile's stretch; any other family (a clothoid, a non-positive spline
  weight) is unbounded by name.

### Changed

- Contour lowering splits a circular segment of half a turn or more into
  equal sub-arcs below half a turn instead of refusing it (#228, ADR 0053
  amended): an IFC arch is commonly one semicircle. The split vertices are
  the circle evaluated at their parameters, as a segment's ends are. Only a
  segment sweeping more than a whole turn is refused.

## [0.3.10] - 2026-10-02

### Added

- `sweep::swept_disk_within`, `sweep::fixed_reference_sweep_within` and
  `sweep::surface_curve_sweep_within`, with `sweep::SampledPath` (#231):
  sweeps that take a directrix sampler instead of fixed samples, and keep
  every point of the surface their section traces along a smooth
  directrix within a chord budget. Each span is bounded by the
  directrix's sagitta, plus the section's extra sagitta at its furthest
  point from the directrix, plus the wall quads' twist; the directrix is
  resampled at half the budget until every span fits. Along a circular
  arc this is exact (the stations are rotations about the arc's axis);
  for other smooth curves it is the second-order estimate. A path with
  end tangents places its end sections square to the curve. A path
  without them (a line, a polyline, a composite) is swept as given.

### Fixed

- `revolve::revolve` and `sweep::tapered_revolve` bound the distance from
  every point of the surface their rings sweep to the mesh by
  `tolerance.linear()` (#231). The step round the axis used to bound the
  sagitta alone, which misses the twist of the walls when the axis leaves
  the profile's plane or the section tapers: the step now also covers
  each wall quad's distance from its two triangles (a proved bound, zero
  for the planar trapezoids of an ordinary revolution), and a tapered
  turn's spiral is bounded by its own curvature. The module notes carry
  the derivation, including why a chorded profile and a chorded turn
  each need their own share of a surface budget. A budget beyond 4096
  steps is refused with `BudgetExceeded` instead of met by fewer steps.

## [0.3.9] - 2026-10-02

### Added

- Exact swept disks (#223). `swept_disk_exact::swept_disk_along_line_exact`
  sweeps a disk, optionally with a bore, along a straight segment into a
  capped cylinder; `swept_disk_along_arc_exact` sweeps it along a circular
  arc into a capped torus wedge, or a whole torus (with a toroidal bore for
  a hollow disk) for a full turn. Both are the exact extrusion and
  revolution of the disk placed with `ExactBRep::transformed`. Volumes
  match `pi r^2 L` and Pappus, every solid audits clean, and certified
  distances under tilted placements match their closed forms. Refused by
  name: a disk reaching the arc's axis, an arc beyond a full turn; invalid
  input: a non-positive radius, a bore not strictly inside the disk, a zero
  span or segment, a non-orthonormal arc frame.

### Fixed

- `revolve::revolve` (the mesh path) no longer builds an inside-out mesh
  for a negative sweep angle (#221). The station order it winds its walls
  from gave a positive signed volume only for a positive angle; since
  `angle` and `axis_direction` only matter through their product (the
  axis-angle rotation vector), a sign-of-`angle` check alone could not be
  made consistent with a sign flip in the axis direction too. The mesh's
  orientation is now settled from its own computed volume instead, so it
  matches `revolve_profile_exact`'s positive sign for either sign of angle
  and either axis direction. `tests/revolve_partial.rs` compared only
  `|volume|` between the two paths because of this; it now compares the
  signed volumes.
- `sweep::tapered_revolve` shared the same unconditional station reversal
  as `revolve::revolve` and so built the same inside-out mesh for a
  negative sweep angle (#221). Fixed the same way: orientation is settled
  from the built mesh's own signed volume rather than guessed from the
  station order.

## [0.3.8] - 2026-10-02

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
