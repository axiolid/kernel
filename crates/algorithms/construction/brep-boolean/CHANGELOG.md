# Changelog

All notable changes to this crate are documented in this file.

The format follows [Keep a Changelog](https://keepachangelog.com/en/1.1.0/),
and this crate adheres to [Semantic Versioning](https://semver.org/).
Pre-1.0: the minor version is the breaking-change slot, per Cargo's own
caret rule for `0.x` versions.

## [Unreleased]

## [0.1.5] - 2026-10-03

### Added

- `BooleanError::ToleranceExceeded` (#251): a decision within tolerance
  was asked for beyond the caller's tolerance, any at all at
  `Tolerance::ZERO`. The report's session refuses it by name (and asserts
  in debug builds) rather than return a result its report misdescribes;
  no reading asks for one.

- `BooleanError::UnsupportedContact` (#243): a plane read as touching a
  cylinder within tolerance is crossed by a curve whose meeting with their
  contact ruling this stage cannot place (a curve of the plane that is no
  line, a curve of the cylinder that is no ruling or conic). Refused by
  name rather than sewn against exact roots that disagree with the
  reading.

### Fixed

- At `Tolerance::ZERO` the report is empty again (#251, the #236
  contract). A support reading the exact predicate rejects was taken when
  its `f64` measure was exactly `0`, and reported as a decision with
  `linear` 0: a round hole across a wall under a general rotation, whose
  axis and the wall's normals round to the same `f64` vector though the
  numbers given (a frame orthonormal only to rounding) are not
  perpendicular, was read `PlanePerpendicularToAxis` at zero. A tolerance
  with a zero part now never takes a reading the exact predicate rejects
  (coincident supports; a plane parallel or perpendicular to a cylinder's
  axis or touching it): the exact answer stands, and the general closed
  form cuts the hole, exactly, with an empty report. Point readings
  already recorded only residues above the rounding floor and within the
  tolerance, never at zero. Both positive parts: unchanged.
- A round web hole touching an I-beam's flange is no longer refused when
  the beam has root fillets (#249). The hole is then tangent to each
  fillet cylinder where the fillet meets the flange, and the hole/fillet
  section (perpendicular cylinders, axes skew by the difference of their
  radii) is a quartic with a double point there: two loops round the
  fillet crossing at the fillet/flange edge. It was refused ("a section
  curve the general boolean does not build"): the trace over the fillet
  face's parameter box met the crossing on the box's edge. Now (new
  `tangency` module):
  - a trace that cannot be decided in a face's box of a periodic carrier
    is taken over the carrier's whole turn, where the crossing is inside
    and both loops end at it (ADR 0077); a loop closing through such a
    vertex is cut there like at any other cut;
  - a section on a face is cut at a smooth edge of the face (a fillet
    running tangent into the flange or web) where the edge itself crosses
    the section's other surface, not where the section touches the
    adjacent surface (a double root that rounding splits or loses);
  - a traced section on the face's own surface keeps the traced curve as
    its pcurve, instead of a second trace over the face's box;
  - with the flange read as touching the hole within tolerance (#243),
    the double point is placed on the contact, and a traced section on a
    cylinder in a contact is cut where it crosses the contact ruling.
  A hole that is exactly tangent (exact axes, dyadic sizes) gives an exact
  result with an empty report at `Tolerance::ZERO`; under a general
  placement the report carries the `PlaneTouchesCylinder` reading. A hole
  a fraction of the tolerance into or short of the flange, which meets
  each fillet in two arcs the contact would have to join (no single move
  of one operand makes both readings hold), is refused by name
  (`BooleanError::UnsupportedContact`). Holes cutting the fillets
  transversally, or ten tolerances into the flange, are decided exactly.

- A plane touching a cylinder within tolerance is read the same way by
  every face pair (#243). A round hole tangent to a planar face (an
  I-beam's web hole touching the flange) placed under a general rotation,
  or reaching a fraction of the tolerance into or short of the face (the
  #234 roof plane in a column), was refused ("split face pieces do not
  close", "point too close to a face boundary", "do not sew"): the
  plane/cylinder pair read one contact ruling while the cap disks' chords,
  the circles of faces across the axis and the plane's own edges kept the
  exact roots, `2 sqrt(2 r d)` apart. The reading is now taken once per
  pair of supports, over the largest common box of their faces, and every
  crossing of a curve on one of the two by the other is placed on the
  contact ruling, as the moved plane gives it (`TangentCrossing`); points
  where the given and the moved plane disagree on a side never decide a
  piece or region. The result is the exact boolean of operands whose plane
  moved by at most the reported `PlaneTouchesCylinder` distance. Exactly
  tangent pairs are still decided by the exact predicates, and succeed at
  `Tolerance::ZERO` with an empty report.
- Cuts merged within tolerance on a circle or ellipse on both sides of its
  parameter origin (one point named at `0` and at `2 pi`) averaged to the
  opposite side of the curve, so a flush hole's circle was split in the
  wrong place and the result did not sew (#243). They are unwrapped to one
  turn first.
- A section's window reaches `2^-20` of its diagonal past the two faces'
  common box (#243). Its planes ran through the faces' extremes, so a
  section touching a face there (a tool resting on a face, at
  `Tolerance::ZERO`) left the window next to the touching point, and the
  piece between could not be classified.

## [0.1.4] - 2026-10-03

### Added

- `BooleanReport::rounding_floor` (#244): the absolute rounding floor that
  applied to a result, in the operands' length unit. Constructed points
  closer than it count as one point and are not reported, so a consumer
  widening distances on a result, exact or not, widens them by it.
  `BooleanReport::extent` is the extent it was scaled from, with its
  definition (the largest `max(|x|, |y|, |z|)` over the operands' vertices
  and their analytic surfaces' origins widened by their radii), and the
  public `ROUNDING_FACTOR` (`2^-40`) is the relative factor. `merged` keeps
  the larger floor.

## [0.1.3] - 2026-10-03

### Added

- `boolean_with_report` returns, with the result, a `BooleanReport` of
  the within-tolerance decisions that fired (#236): one
  `ToleranceDecision` per `ToleranceDecisionKind` (coincident supports,
  contact, a point on an edge, a section along an edge, a tangent crossing,
  merged points, an iso-curve, a plane parallel or perpendicular to or
  touching a cylinder), with the furthest it moved (`linear`) or turned
  (`angular`) the operands. An empty report means the result is the exact
  boolean of the operands as given. `boolean` is unchanged.

### Fixed

- Placed differences succeed at `Tolerance::ZERO` (#236). Naming a point
  on the curve or surface it was evaluated from compared its `f64`
  round-trip residue with the linear tolerance, so at zero every boolean
  failed with `BooleanError::Evaluation`; that bookkeeping, and welding
  two evaluations of one vertex, now allow the rounding of the
  operands' extent (`2^-40` of it).
- Readings about the operands' own surfaces (coincident supports, a plane
  parallel or perpendicular to or touching a cylinder) ask an exact
  dyadic predicate first: faces exactly coplanar, parallel or
  perpendicular, as operands placed by matrices with entries `0` and
  `+-1` meet, are decided exactly and read nothing within tolerance.

## [0.1.2] - 2026-10-03

### Added

- `BooleanError::NearCoincidence`: features chained within tolerance of
  each other over more than the tolerance, which no single perturbation
  within tolerance reconciles, are refused by name (#228).
- The crate docs state what a result guarantees: which decisions are
  exact, which are taken within the caller's tolerance, and that the
  result is then the exact boolean of operands moved by at most it.

### Fixed

- Operands placed by independent rigid motions cut each other (#228). Their
  coincident, parallel or tangent faces agree only up to rounding, which
  the exact predicates saw: an opening flush with a wall's faces was
  undecided, an arched opening's tangent jambs left slivers or no cut, and
  a plane `1e-17` off parallel or perpendicular to a cylinder's axis gave
  an ellipse nothing could evaluate. Within tolerance, a section now runs
  along an edge it lies on, a tangent crossing is cut once where it meets
  the edge, cuts and boundary splits closer than tolerance are one point,
  and such a plane cuts the cylinder in rulings or a circle.

### Changed

- Every within-tolerance decision uses the caller's `Tolerance` only:
  distances its linear part, with no built-in `1e-9` floor, and the
  coincidence and contact tests on normals and axes its angular part
  instead of a fixed `1e-9`. With `Tolerance::METRE` results are
  unchanged.
- Work between faces that cannot meet is skipped (#228). Each face and
  edge gets a sound box from its surface or curve in closed form, enlarged
  by the linear tolerance; face pairs whose boxes are apart are not
  sectioned, a section line or conic is cut only against edges whose
  boxes meet both faces' common box (and at that box's boundary),
  classification rays skip faces whose boxes they miss, and split faces
  ignore cut points outside their boxes. A wall losing ten placed windows
  one at a time went from 1.67 s to 0.28 s in release
  (`cargo bench -p axiolid-benchmark --bench exact_openings`).

## [0.1.1] - 2026-09-28

### Changed

- The crates.io page is this crate's own `README.md`, with links to its
  API documentation, its reference page and the source (ADR 0078).

## [0.1.0] - 2026-09-27

### Added

- Surfaces tangent along a whole curve and crossing there: the line of
  contact is a section edge (the tangent-contact rule now applies to
  closed forms only), and pieces leaving a vertex with the same direction
  and bend are ordered by their chords a small way out.
- A traced pair that only touches adds no section, as the closed forms'
  touching does. A pcurve whose trace on its face cannot be decided falls
  back to the section's own space curve read on the face, and failed
  traces are not repeated. Each surface pair's closed form is computed
  once per boolean.
- Two B-spline faces meeting each other: their section is traced in both
  faces' parameter boxes and carried on both (`Curve3::PairSection`), its
  pcurve on each read from the solve. The last refusal by face type is
  gone; sections whose branches cross where the surfaces touch stay
  refused by name.
- Faces that wind round their surface without a seam edge (a dome bounded
  by its rim alone, a can by its two rims), as files may deliver them. The
  boolean first gives each a seam edge along the iso-curve where its loops
  wrap, joining the lower loop, the seam, the upper loop or a pole, and the
  seam back into one loop.
- Sections through a sphere's pole or a cone's apex are cut there, and the
  collapsed pole piece is split where they end. Meridians and latitudes,
  cone rulings and circles, and a torus's tube and ring circles get their
  straight pcurves, affine in the curve's own parameter.

- B-spline faces (#167, ADR 0075 stage 3) against planes, quadrics and
  tori:
  - The section is traced on the spline (ADR 0077).
  - On the analytic face its pcurve is the same space curve read in the
    face's parameters (`Curve2::Lifted`), so it shares the edge's
    parameter.
  - An edge next to a B-spline face is cut where it meets the section's
    other surface.
  - Classification rays meet B-spline faces through the spline trace.
  - Two B-spline faces meeting each other are refused by name
    (`UnsupportedSection`).

- Faces on spheres, cones, tori and elliptical cylinders (#167, ADR 0075
  stage 2), meeting in any section #119 builds:
  - A section with no line or conic is traced inside one face's parameter
    box (ADR 0077).
  - Every section on every analytic face gets an exact implicit pcurve, cut
    out of the other surface's traced equation between the section's ends.
  - A sphere's pole or a cone's apex closes loops as a collapsed piece that
    is no edge.
  - Seam circles are cut by the cone of normals along them.
  - Section branches that cross (a Steinmetz pair) are split where they
    meet.
  - Frame components that are only rounding residue are cleared before
    intersecting.

- Operands that touch (#167, ADR 0075 stage 2): faces on one surface share
  their overlap (each face's edges are imprinted on the other, and a region
  on the other solid's boundary is kept once by normal agreement); sections
  along an existing edge split only the other face; tangent contact adds no
  section; pieces leaving a vertex in one direction are ordered by
  curvature; solids meeting along an edge are paired radially around it so
  each stays manifold. Cavities go to the smallest solid around them, in
  results of several solids too.

- `section_edges` (#167, ADR 0075 stage 1): the exact intersection curves of
  two exact B-reps' faces, each trimmed to where it lies inside both faces.
  Crossings with a boundary edge are found against the adjacent face's
  surface, or across a seam against the plane through the ruling.
- `boolean(a, b, operator, tolerance)` (#167, ADR 0075 stage 1): the exact
  union, intersection and difference of two exact solids whose faces lie on
  planes and cylinders and meet in lines, circles and ellipses -- not only
  vertical columns. Regions are classified by exact ray parity with
  certified face membership and sewn into shells; cavities become voids.
  Every result audits clean and measures exactly.
- `split_face` (#167): a plane or cylinder face cut along its section edges
  into regions, traced in the face's parameters with exact pcurves (lines,
  conics, rulings, circles about the axis, `Sinusoid2` for oblique cuts).
