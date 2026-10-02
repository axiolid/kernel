# Changelog

All notable changes to this crate are documented in this file.

The format follows [Keep a Changelog](https://keepachangelog.com/en/1.1.0/),
and this crate adheres to [Semantic Versioning](https://semver.org/).
Pre-1.0: the minor version is the breaking-change slot, per Cargo's own
caret rule for `0.x` versions.

## [Unreleased]

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
