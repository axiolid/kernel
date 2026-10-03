# Changelog

All notable changes to this crate are documented in this file.

The format follows [Keep a Changelog](https://keepachangelog.com/en/1.1.0/),
and this crate adheres to [Semantic Versioning](https://semver.org/).
Pre-1.0: the minor version is the breaking-change slot, per Cargo's own
caret rule for `0.x` versions.

## [Unreleased]

## [0.3.11] - 2026-10-03

### Added

- Sweep directrices resolve `TrimSelector::ArcLength` (#239) against an
  analytic basis, in the trim's sense, through
  `axiolid_reference::arc_parameter`; the mesh, piecewise and exact sweep
  paths read it like a parameter selector. An arc length on a
  curve-relation basis, or past the end of its basis, is refused by name.

## [0.3.10] - 2026-10-03

### Added

- `ReferenceExactCompiler` clips by half-spaces exactly (#234): a
  difference or intersection whose tool is a `HalfSpace`, or a
  `BoundedHalfSpace` (a polyline boundary in the plane), read through
  rigid placements, with the mesh compiler's semantics (`agreement` selects
  the normal side; the boundary is framed by its placement projected into
  the plane and swept along the normal). The half-space becomes a finite
  prism over a sound envelope of the subject plus a margin and goes through
  the general exact boolean under the #228 tolerance contract; the result
  does not depend on the margin (ADR 0080, amended). Clips compose with
  each other and with placed openings in either order: a wall minus its
  windows, clipped by two roof planes. Unions with a half-space, a
  half-space as the subject, a subject not built from placed extrusions, a
  boundary that is not a polyline, scaled placements and the general
  boolean's refusals are refused by name; a clip that removes everything
  is `Degenerate`.
- `ReferenceExactCompiler::compile_exact_with_report` and
  `compile_exact_batch_with_reports` return each body with a
  `BooleanReport` (re-exported from `axiolid-brep-boolean`, with
  `ToleranceDecision` and `ToleranceDecisionKind`) merged over every
  general boolean and clip beneath it (#236). An exact report means the
  body is the exact result of its operands as given and may be cited as
  exact; otherwise it bounds how far the operands were moved or turned,
  within the tolerance.

### Fixed

- Placed differences compile at `Tolerance::ZERO` (#236): openings placed
  with exact axis matrices (entries `0` and `+-1`), through, blind or
  flush, give the exact difference and an exact report.
- A boolean's deviation contribution says when its bound also holds for
  the exact boolean of the given operands (#236): measured against an
  exact compiler result whose report is exact, its detail reads "the exact
  boolean of the given operands" instead of "operands within tolerance".

- `ReferenceMeshCompiler::compile_mesh_with_deviation` reports a boolean
  `Certified` where `ReferenceExactCompiler` builds its exact result
  (#235): differences of placed extrusions, such as a wall with a round
  window, an I-beam with round holes through its web, or a slab with a
  round shaft, and walls clipped by roof half-spaces (#234). The boolean's
  mesh is measured against that exact B-rep by the certified branch and
  bound of #232, over each face's trimmed parameter domain (pcurves
  flattened to a certified chord bound; a cell is dropped only when it is
  certainly outside the face). The bound is relative to the exact
  compiler's result, which is the exact boolean of operands moved by at
  most the tolerance (#228), and is not claimed for the boolean of the
  unperturbed operands; the contribution's detail says so ("measured
  against the exact compiler's result, operands within tolerance").
  Nothing is derived from the operands' bounds, which are one-sided and
  say nothing about where the mesh boolean puts the cut. A boolean the
  exact compiler refuses stays `Unbounded`, named by the refusal (for
  instance "exact union or intersection of placed operands"), detail "no
  exact result". Only a deviation report pays for this, only for the
  booleans whose result is emitted (not the inner differences of a chain),
  and the search stops once the bound is within the requested budget.
  Release build, 1 mm budget: a 6 m wall with a 0.4 m round window in
  about 0.5 s, a 4 m I-beam with three round web holes in about 0.9 s.

### Changed

- The certified branch and bound covers a piece with a flat region of the
  mesh (edge-connected, consistently wound, coplanar triangles) as well as
  with single triangles, so edges inside a planar face or a cylinder facet
  no longer have to be resolved to the bound's own size. Bounds are as
  sound as before and settle with less work.

## [0.3.9] - 2026-10-03

### Added

- `ReferenceExactCompiler` compiles differences of placed extrusions
  exactly (#228): a wall or slab under any rigid placement minus openings
  under theirs, perpendicular to its extrusion (doors, windows) or
  parallel (shafts), with rectangle, circle and line-and-arc profiles,
  through, blind, flush or touching its edges, and several openings per
  body as nested differences. It runs the general exact boolean
  (`axiolid-brep-boolean`, ADR 0080); two unplaced sharp rectangles along
  `+z` keep the prism path. Unions and intersections of placed operands,
  operands that are not extrusions, a tool that is itself a boolean, scaled
  placements and configurations the general boolean refuses are refused by
  name; a difference that removes the whole subject is `Degenerate`.
  Faces that agree only up to rounding are read within the caller's
  tolerance; the result is then the exact boolean of operands moved by at
  most it, and a gap of ten tolerances is kept. One-segment semicircular
  arches compile (contour lowering splits them).

- `ReferenceMeshCompiler::compile_mesh_with_deviation` and
  `DeviationReport` (#232): next to the mesh, a certified upper bound on
  the distance from every point of the exact surface to the triangles, the
  paths that contributed (`DeviationContribution`, `DeviationPath`,
  `DeviationBound::{Proven, Certified, Unbounded}`) and whether the
  requested chord budget is met. #231's paths report the budget they are
  proven to; profiles report their flattening's bound (a derived profile
  its stretch); curved B-rep faces a per-triangle bound from the surface's
  second-derivative bounds plus each trim pcurve's lens; pipes along
  segments and arcs and primitive cylinders and cones the budget their
  constructions prove; disks swept along B-splines and ellipses a bound
  certified against the exact tube by branch and bound. Booleans, tapered
  extrusions, sectioned spines, bounded half-spaces, composites holding
  other curves and other frame laws are unbounded by name. An
  inherent method: the `MeshCompiler` contract is unchanged.

### Changed

- A swept disk's `fillet_radius` is honoured on polylines and composites
  of lines (`IfcSweptDiskSolidPolygonal`) instead of refused: each corner
  becomes a tangent arc of that radius (#232). A directrix with a corner
  and no fillet radius, which IFC leaves undefined and which used to be
  swept with sharp mitres, is now refused by name, as are a closed
  polyline directrix, a fillet that does not fit its segments, a disk
  radius at or above the fillet or bend radius, and a sweep range
  combined with a fillet radius (`UnsupportedInput`).

- A B-spline directrix with no corner knot reports its exact end
  tangents, so sweeps along it refine their stations and stand their end
  caps square to the curve, as along a conic (#232).
- Curved B-rep faces are also refined where a triangle's certified bound
  misses the chord budget, by its widest free edge, up to a vertex cap; a
  pass that runs out of depth or vertices is discarded for the
  measurement-only refinement. Trim edges are sampled until their
  certified chord bound fits too, where that converges (#232).

### Fixed

- A swept disk along a line, a polyline, or a trim or composite of lines,
  polylines and circular arcs keeps every point of the exact tube within
  the chord budget of the mesh (#232). The directrix is read as exact
  segments and arcs and swept by `axiolid_construct::pipe` (see its notes
  for the bound); it used to be sampled as a whole, each bend chorded for
  its centreline instead of the tube's outer side. At 1 mm: line + bend
  (R 0.1) + line, r 0.05, 1.184 -> 0.784 mm (502 -> 686 triangles);
  three bends out of plane 1.184 -> 0.803 mm (1652 -> 1882); at 0.1 mm
  R 0.04 / r 0.01 0.108 -> 0.077 mm, R 0.6 / r 0.3 0.115 -> 0.085 mm.
  A composite's sweep range is now cut at exact arc lengths rather than
  along its sampled chords. Lone circles and ellipses keep the #231
  sweep, and composites holding other curves (ellipse arcs, B-splines)
  are swept as sampled.

- Trim samples earcut skips as collinear are put back on a curved face
  (#232): a straight pcurve sampled into many points came back as one long
  triangle edge, a T-junction against the face across it. The fan of thin
  triangles that putting them back leaves is then flipped to a Delaunay
  triangulation of the trim polygon, in parameters scaled to the surface.

## [0.3.8] - 2026-10-02

### Added

- `exact_directrix` and `ExactDirectrix` are public (#230): they read a
  swept disk's directrix as one segment or one arc, exactly as
  `ReferenceExactCompiler` does, so a consumer building exact boundaries
  that must match the compiled solids no longer duplicates the reading.
  Corners, other curve families and unbounded lines are refused by name.

### Fixed

- Doubly curved meshes stay within the chord budget (#231). A revolution
  chorded its profile and its turn each to the whole budget, so the two
  deviations added inside a triangle: a torus (R 0.5, r 0.1) at 1 mm lay
  1.46 mm from its mesh. A revolution, a tapered one included, now chords
  its profile to half the chord budget and its turn to the other half
  (see `axiolid_construct::revolve` for the proof), and its turn follows
  `ExecutionOptions::with_chord_error` instead of the linear tolerance.
  Swept disks, fixed-reference and surface-curve sweeps along a circle or
  ellipse, plain or trimmed, give the section half the budget and refine
  the directrix until the section's far side fits the other half, and
  stand their end sections square to the curve rather than to its end
  chords (a disk r 0.1 along an arc R 0.2 lay 7 mm from its end caps at
  1 mm). Every point of the exact surface now lies within the budget of
  the triangles, checked by dense sampling of tori (R/r down to 1.2),
  spheres, partial and skew-axis revolutions, a revolved rounded
  rectangle, a tapered revolution and sweeps along arcs with r/R up to
  0.83, at 1 mm and 0.1 mm. A torus at 1 mm has 4928 triangles instead of
  3520. A budget beyond 4096 steps round an axis is refused with
  `BudgetExceeded` rather than met by a coarser mesh. Polyline, composite
  and B-spline directrices are swept as sampled, as before.

## [0.3.7] - 2026-10-02

### Added

- `ReferenceExactCompiler` compiles `Instance` nodes (#223): the source's
  exact B-rep is placed with `ExactBRep::transformed`, so extrusions,
  revolutions and swept disks under rotated, reflected and translated
  placements (nested instances compose) stay exact. A scaled or sheared
  instance is refused by name, never approximated.
- `ReferenceExactCompiler` compiles `SweptDisk` exactly along one segment
  or one arc (#223): a bounded line, a two-point polyline, a circle or a
  sub-range of one, a trim of a line or circle (across the seam, as the
  mesh path reads it) and a one-segment composite. A directrix with
  corners, any other curve, an unbounded line and a disk reaching its
  arc's axis are refused by name.

## [0.3.6] - 2026-09-28

### Changed

- The crates.io page is this crate's own `README.md`, with links to its
  API documentation, its reference page and the source (ADR 0078).

## [0.3.5] - 2026-09-28

### Fixed

- A boolean whose result touches itself is refused, not returned (#194):
  where operands meet tangentially -- a void tangent to its host's face --
  the solid has no material between two faces, and the mesh boolean keeps
  two copies of the vertices there, closed by index but pinched by
  position. Consumers welding by position saw an edge with four faces. The
  result is now checked, positions welded, for an edge with more than two
  faces or a vertex with separate fans, and refused with
  `GeomError::Degenerate` naming the edge or point of contact. Circular
  voids tangent along an axis direction no longer pinch at all (construct).
- Structural sections and rounded rectangles extrude to meshes (#193, via
  construct).

## [0.3.4] - 2026-09-27

### Added

- Curve-bounded planes (#192), as IFC `IfcCurveBoundedPlane` space-boundary
  connection surfaces carry them: a `SurfaceRelation::CurveBounded` over a
  planar basis compiles to a triangulated planar region with
  `MeshClosure::Surface`. Boundaries are read in the plane's parameters and
  mapped through its frame; the first is the outer loop, the rest holes.
  Straight boundaries (2D or 3D polylines, 3D ones on the parameter plane)
  are exact, curved ones chorded within the chord budget; the mesh faces
  along the plane's normal whatever the loop's winding. Crossing or
  degenerate loops, holes outside the outer loop and open boundaries are
  refused; a non-planar basis and `implicit_outer` get `UnsupportedInput`
  naming the capability.

### Changed

- Every solid of a B-rep is tessellated, not only the first (#111): a
  composite profile whose members do not touch is several solids in one
  B-rep.
- A solid's void shells are tessellated (#120). They were dropped as
  "boolean intent", which silently filled every cavity of an authored
  `IfcFacetedBrepWithVoids`-style B-rep and forced the exact booleans to
  refuse cavities. Each void is emitted facing into the cavity, so the mesh
  encloses the outer volume less every cavity; a void authored facing out
  of it (the STEP convention, reversed on use) is turned round, since a
  cavity can only remove material. A void shell that is not closed is
  refused rather than leaving the mesh open.

## [0.3.3] - 2026-09-25

### Added

- The reference compiler honours `ExecutionOptions::with_chord_error`
  (#165) everywhere it flattens a curve: profile arcs, circles and ellipses,
  sweep directrices, curved B-rep faces and edges, and CSG primitives. An
  instance scales the budget with its transform, like the tolerance, so it
  stays a world-space distance. Without a budget the chord error is the
  linear tolerance, exactly as before. Measured on a 5 mm disc extruded 1 m
  at `Tolerance::MILLIMETRE`: 10 % short by default, 2.6 % at a 0.1 mm
  budget, 0.16 % at 10 um and 0.01 % at 1 um.

### Fixed

- A surface model with a face whose outer bound encloses no area
  tessellates (#171): that face covers nothing, so it is skipped instead of
  refusing the whole model with `planar face bound has zero or non-finite
  area`. Real Nova MEP exports write pipe-fitting end caps as bowtie quads
  through the pipe axis (an annulus with a negative inner radius), each with
  signed area exactly 0; 42 fittings, pumps and valves in two models were
  refused over them. A declared solid still refuses such a face, and a
  zero-area hole or a non-finite bound is still refused everywhere.
- Closed authored meshes stay closed (#170). Planar faces of a
  `PolygonMesh` and of a B-rep were triangulated with earcut, which drops
  corners on a straight run and runs diagonals and hole bridges over corners
  of the same face; the neighbouring face still split that edge at the
  corner, so the mesh cracked (T-junctions). Every edge earcut invents is now
  split at each face corner on it, within a band of a thousandth of the
  linear tolerance (1 um at `Tolerance::MILLIMETRE`), so export noise on
  shared corners (1e-8 to 5e-8 m on real files) is judged the same on both
  sides. Authored ring edges are never split, and a thin triangle whose long
  side is authored is kept. Closure of an authored mesh is now read from its
  index connectivity instead of `audit_mesh`, which dropped real faces below
  its area threshold before counting edges. On two real ArchiCAD models this
  turns 505 authored-closed `IfcPolygonalFaceSet` products from `Surface`
  into `Solid`; no product that compiled before is refused.

## [0.3.2] - 2026-09-25

### Fixed

- A directrix trimmed from a circle or ellipse ACROSS its seam sweeps the arc
  the trim names (#168). The trimmed curve runs from `start` the way
  `sense_agreement` says, wrapping past the seam if it must; the directrix
  path used to sort the two trims and sample the complementary arc. Standalone
  that was silently wrong geometry (a `315 -> 45` degree bend swept the 270
  degree arc); inside a composite the ends no longer met and the sweep was
  refused as `composite directrix has a N unit gap`. In a real Revit rebar
  model that was 1,494 bent bars, all writing a bend as `270 -> 45`,
  `270 -> 15` or `270 -> 360` rounded just past the seam. A full turn rounded
  past the seam stays a full turn. A sweep `parameter_range` on such a trim is
  read in the trim's unwrapped interval, so `(330, 30)` and `(330, 390)`
  degrees name the same sub-arc of a `315 -> 45` trim; an end off the arc is
  refused as before. Profiles already honoured this.

## [0.3.1] - 2026-09-24

### Added

- `PolygonMesh` faces that are not plain triangles compile (#160): n-gons,
  concave faces and faces with holes (IFC4 `IfcIndexedPolygonalFaceWithVoids`)
  are triangulated in their own plane, keeping the authored positions and
  winding. Plain triangles keep their exact corner order as before. A face
  whose corners leave its plane by more than the linear tolerance, that has
  no area, or whose rings cross is refused with an error naming its index.
- B-reps with shells but no solid tessellate (#161): every shell is
  tessellated as authored and the result is reported as
  `MeshClosure::Surface` through `compile_mesh_reported`, even when the
  shell is closed. Collections are `Solid` only if every member is, and a
  boolean with a surface operand is refused. Authored meshes report
  `Solid` exactly when they are closed, consistently wound two-manifolds.

### Changed

- A `PolygonMesh` with non-triangular faces used to fail with
  `GeomError::Unsupported`; it now compiles. A B-rep with no solid and no
  shell is refused as "neither a solid nor a shell" instead of "has no
  solid".

## [0.3.0] - 2026-09-23

### Fixed

- `Instance` and `Collection` nodes no longer drop attribute channels and normals (#115). Both used to rebuild the mesh from positions and indices only, so a textured item lost its `uv` channel as soon as a product had a second item or was instanced — silently.
  - `Instance` carries channels through unchanged and now transforms normals by the inverse transpose instead of dropping them. Under a mirroring transform, corner-indexed channels and normals swap corners with the triangle.
  - `Collection` merges channels by name. A channel on only some members becomes corner-indexed, with the other members' triangles `UNMAPPED`; a channel on every member as per-vertex stays per-vertex. Members defining one name with a different width or blend drop it as `DropReason::IncompatibleChannels`.
  - Booleans keep the provider's channel fates, composed onto what each operand already went through, instead of discarding the evidence.

### Added

- `ReferenceMeshCompiler` implements `MeshCompiler::compile_mesh_reported`: the compiled mesh plus each channel's fate on its way to the root (worst over parallel members, sequential through booleans).
