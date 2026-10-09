# Changelog

All notable changes to this crate are documented in this file.

The format follows [Keep a Changelog](https://keepachangelog.com/en/1.1.0/),
and this crate adheres to [Semantic Versioning](https://semver.org/).
Pre-1.0: the minor version is the breaking-change slot, per Cargo's own
caret rule for `0.x` versions.

## [Unreleased]

### Changed

- `boundary_distance`, `boundary_clearance` and the `body_*` distance and
  clearance in space stop as soon as the two boundaries are shown to meet
  (#273): a pair the search pops with a lower bound of zero is first asked
  for a point of each element that the two share -- the closest points of
  two line edges, or a point of a line edge and the point of a planar face
  at the same place, the face's certified by its domain. Two such points
  within the rounding margin the lower bounds carry give `[0, d]` with
  them as the witnesses, the interval the refinement would only converge
  to. Touching or crossing elements no longer refine round the contact
  until the accuracy or the step budget: two boxes sharing a face or
  crossing stop within a few steps, and on the `boundary_distance`
  benchmark a wall standing on a slab under a general placement measures
  in about 35 us instead of about 0.3 s (where the old search also ran out
  of budget at `[0, 1e-3]`), two crossing walls in about 55 us instead of
  27 ms. A gap wider than the rounding margin is never taken for a touch,
  and boundaries that are apart measure bit for bit as before. Curved faces
  and circle or ellipse edges show no touch this way and are refined as
  before.

## [0.3.9] - 2026-10-03

### Added

- Plan distance, clearance and overlap between bodies of several exact
  solids (#237): `body_plan_boundary_distance`,
  `body_plan_boundary_clearance` and `body_plan_overlap` take a
  `PlacedBody` (or a slice of items) per side, as the `body_*` queries in
  space. The plan distance is the least over item pairs of the distance
  between their shadows on the XY plane, by one shared search, so far
  pairs are never refined; it is zero where an item's shadow overlaps or
  lies inside one of the other body's, and the witnesses name their items
  (`BodyDistance`). `body_plan_overlap` answers `BodyPlanOverlap`: an
  overlap of positive area with the two items that show it, a certified
  gap, or undecided. Items of one body may overlap or touch freely in
  plan, so no layout is checked. The placement is applied before
  projecting along the world's `z`: a body turned about `z` measures as
  its plan turned, and a tilted body casts its tilted items' shadows.

### Fixed

- `plan_overlap` no longer leaves undecided an overlap that two level
  planar faces show when their shadows overlap only in part (a block over
  part of a column's disc): when the distance search shows neither an
  overlap nor a gap -- its budget taken by walls and edges whose shadows
  merely cross, or its run ended by a plan distance of zero met between
  two boundary points -- a second search over patches of planar faces
  that are not vertical alone, coarsest pair first, looks for one.

## [0.3.8] - 2026-10-03

### Added

- Boundary distance and Hausdorff distance between bodies of several
  exact solids (#229), one `ExactBRep` per item in a slice:
  `body_boundary_distance` and `body_boundary_clearance` (the least
  distance over item pairs, by one shared branch and bound, so far pairs
  are never refined), `body_boundary_hausdorff_distance`,
  `one_sided_body_boundary_hausdorff` and its `_with_budget` form (the
  Hausdorff distance between the boundaries of the unions). Witnesses name
  their items (`BodyDistance`, `BodyHausdorffBounds`, `BodyHausdorff`).
  The Hausdorff queries first show every pair of items of a body apart (a
  plane with a certified positive gap between them, or boundaries
  certainly apart with neither inside the other), touching without a
  shared patch of face (a plane between them within rounding where one of
  them has no face), or in exact face contact. Items in face contact (a
  column on its footing, blocks sharing a wall, stacked steps) share a
  patch that is interior to their union: each face in contact is cut down
  to its free region by an exact arrangement of the faces' boundaries on
  their common plane (`axiolid-overlay`'s `ArcArrangement`, now an
  optional dependency under `exact`), and the union's boundary is
  measured. Contact is exact when both faces lie, in the B-rep's own
  numbers, on one plane normal to a coordinate axis, which a placement
  turning about that axis keeps. A gap is never glued: two items a
  sub-millimetre apart are measured as apart, both facing faces on the
  boundary, whatever the tolerance (gluing and widening the interval by
  the gap would not be sound: a lifted column's base disc is a whole
  radius from the glued boundary). Anything else is refused by name in
  the new `BodyMeasureError`: `ContactPlaneNotAxisNormal` (items touching
  on a plane no coordinate axis is normal to, such as walls turned in
  plan), `ItemsNearlyShareFace` with the gap (items interpenetrating by
  no more than the caller's tolerance, or faces on axis planes at
  coordinates differing below rounding), `ItemsShareFace` (a face in
  contact that cannot be cut: a B-spline face, an elliptical edge) and
  `ItemsOverlap`.
  Every `body_*` query takes each side as a `PlacedBody`: the items in the
  body's own frame and one rigid placement (a slice, vector or array of
  items converts with the identity). Contact is found and cut in the
  body's frame, where an IFC body's extruded items stand on axis-normal
  planes, and the placement then moves the cut boundary, so a body turned
  in any direction keeps its contact cut; witnesses are in the world.
  Items turned against each other within the body's frame are still
  refused as `ContactPlaneNotAxisNormal`, and a placement that is not
  rigid as `Placement`.
  Distance needs no such check. A translated multi-item body closes as
  fast as one solid: faces match across items, a free region is bounded
  through the face it was cut from (so a cut that falls out differently
  in the last bit for the moved copy still closes at once), and each
  item's own support point seeds the lower bound. Ten stacked steps
  turned in plan and moved close both ways at 1e-9 in about 0.14 s.

## [0.3.7] - 2026-10-02

### Added

- `one_sided_boundary_hausdorff_with_budget` (#227): the one-sided
  boundary Hausdorff distance with the caller's cap on the splits, so a
  consumer bounds the work per pair. The interval is sound at any budget,
  only wider when it runs out.

### Fixed

- `boundary_hausdorff_distance` took its matched fast path only for faces
  with identical pcurves (#227). Prisms from `boolean_arc_prisms_exact`
  are trimmed in world coordinates, so two of them a translation apart
  were treated as unmatched and closed at first order: [1.0, 1.59] mm at
  accuracy 1e-4 after the whole 200k-split budget. Faces are now also
  matched as translates, independently of how the B-rep was built: same
  family, axes and shape (radius, semi-axes, angle, a cone's apex), and
  trims equal after the parameter shift the translation induces (for a
  plane the shift across its axes, for a cylinder, elliptical cylinder or
  cone along its axis only), to within a gate of 1e-9 relative. The other
  face is re-charted by the shift and bounded as before, the measured trim
  residue folded in, so the bound is `|t|`; a turned or resized face, by
  however little, is never matched this way. The lower bound is seeded
  with the support point against each matched displacement, exactly `|t|`
  from a translate's boundary. Those prisms, square, round or with an arc
  in the section, moved by 1 mm or 0.2 m in any direction, now close to
  1e-6 without a split.

## [0.3.6] - 2026-10-02

### Added

- Certified Hausdorff distance between exact boundaries (#224), behind the
  `exact` feature. `boundary_hausdorff_distance` returns
  `BoundaryHausdorff`: intervals certain to contain the two-sided
  Hausdorff distance between the boundaries of two `ExactBRep`s and both
  one-sided ones, refined to a requested absolute accuracy;
  `one_sided_boundary_hausdorff` measures one direction. Each interval is a
  `HausdorffBounds`, as for meshes, with witnesses: a point certainly on
  the measured boundary realising the lower bound, and its nearest point
  found on the other. Lower bounds are such points' distances bounded below
  by the certified distance search; upper bounds come from branch and bound
  over face patches, each bounded by a matched face of the other boundary
  (same family, same trim: `|S_A - S_B|` over the patch, exact for a
  translation) or by the 1-Lipschitz distance about a point. Copies offset
  by a translation and identical copies close in a few dozen splits at any
  accuracy (two columns 0.1 mm apart to 1e-6, an arched opening moved
  0.2 m, a re-export within the accuracy); other pairs close at first order
  in the patch size. A fixed budget can end refinement early; the interval
  stays sound, only wider.

## [0.3.5] - 2026-09-30

### Added

- Certified mesh Hausdorff distance (#148). `hausdorff_distance` returns
  `MeshHausdorff`: intervals certain to contain the two-sided Hausdorff
  distance between two triangle-mesh surfaces and both one-sided ones,
  refined to a requested absolute accuracy; `one_sided_hausdorff` measures
  one direction. Each `HausdorffBounds` carries witnesses: the sample
  realising the lower bound and its nearest point on the other mesh. Lower
  bounds are sampled points' distances bounded below through each
  triangle's support function; upper bounds come from branch and bound over
  subdivided triangles (the distance to a triangle is convex, so a piece's
  maximum is at a corner) with a best-first BVH, and flat convex patches of
  the target (edge pairs, closed fans) bound pieces across seams with a
  certified hull excess. Rounding is accounted for with explicit margins,
  including the drift of rounded subdivision midpoints. Open, non-manifold
  and degenerate meshes are measured; empty meshes, bad indices, non-finite
  positions and invalid accuracies are refused with `HausdorffError`.
- Hausdorff distance between polylines (#147), one-sided and two-sided,
  in 2D and 3D: `polyline_hausdorff_distance`,
  `one_sided_polyline_hausdorff_distance` and their `_2d` forms. The
  supremum over a segment of the distance to the other polyline is found
  among its ends and the points where two features of the other polyline
  -- vertices, segment interiors -- are equally near, since the distance
  to one feature is convex along a line.
- `frechet_decide_certified` and `frechet_decide_certified_2d` answer
  whether the Fréchet distance is at most `eps` only when rounding cannot
  change the answer: `FrechetDecision::AtMost`, `MoreThan`, or
  `Undecided` within the error margin of the floating-point decision.

## [0.3.4] - 2026-09-30

### Added

- Plan relations between exact bodies (#217). `plan_boundary_distance`
  and `plan_boundary_clearance` certify the distance between the bodies'
  projections onto the XY plane, with the contract of `boundary_distance`:
  the interval contains the distance, and the witnesses lie on the
  boundaries with their projections `upper` apart. A shadow is its
  boundary's, so the same search runs with horizontal gaps, directions and
  enclosing discs; the result is zero where the shadows overlap, including
  a body standing inside another's footprint. `plan_overlap` returns
  `PlanOverlap::Overlapping` with a plan point inside both shadows when two
  planar faces, not vertical, share an open patch in plan (a column on a
  slab), `Disjoint` with a certified gap, or `Undecided`.

## [0.3.3] - 2026-09-28

### Changed

- The crates.io page is this crate's own `README.md`, with links to its
  API documentation, its reference page and the source (ADR 0078).

## [0.3.2] - 2026-09-27

### Added

- Fréchet distance between polylines (#147): `frechet_distance` (the
  continuous distance: the least critical value of the free space that the
  Alt-Godau decision accepts), `discrete_frechet_distance` (Eiter-Mannila,
  `O(nm)` time, `O(m)` memory) and the decision `frechet_at_most`, each
  with a `_2d` form. Empty polylines, non-finite points and an invalid
  leash are `FrechetError`s. Floating point, not certified.

## [0.3.1] - 2026-09-27

### Added

- Face domains trimmed by pair sections of two B-splines: the section's
  nodes are kept as breaks, where its parameter speed changes.
- Face domains for B-spline and lifted pcurves. Turning points come from a
  dense scan of the derivative's signs with bisection, and a B-spline's
  knots are kept as breaks.

- Turning points of `Curve2::Implicit` pcurves, isolated with interval
  bounds, so faces trimmed by them get a certified domain.

### Fixed

- `FaceDomain::contains` counts winding -1 as inside (a face whose loops
  run clockwise in its parameters). It also moves a point by whole periods
  of a torus's second angle into the face's box.
- Face integrals split an implicit pcurve at its cell boundaries, where
  its parameter's speed jumps, so adaptive quadrature converges there.

- `FaceDomain::contains` decides a point whose ray passes through a
  boundary vertex or along a boundary piece, by simulation of simplicity
  (the ray runs just right of the point), instead of answering `None`. The
  centre of a diamond, whose every axis ray meets a corner, used to be
  undecidable.

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

- `ExactMeasureError` is unchanged in shape (exhaustive, four variants), so
  this is a patch release. A face the module cannot integrate -- a surface
  family it has no integral for, pcurves that bound no domain, a surface it
  cannot evaluate there, an integral short of its error bound -- is
  `NonPlanarFace`, whose name predates curved faces, with the reason. The
  display text reads "cannot integrate a face (reason)".
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
