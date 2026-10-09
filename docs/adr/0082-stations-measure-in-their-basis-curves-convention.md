# 0082 — Stations measure in their basis curve's convention

- **Status:** Accepted
- **Date:** 2026-10-03
- **Deciders:** Friedrich Schrödter
- **Supersedes:** —

## Context

Issue #241 (consumer: openbimrs/ifc#307). Alignment geometry names
places by a distance along a curve with offsets in the curve's frame:
a point or placement at a station, an offset curve whose lateral and
vertical offsets are given at several stations, a solid whose cross
sections stand at stations and are interpolated between them, and a
surface of open sections joined point by point. IFC 4.3 carries these as
`IfcPointByDistanceExpression`, `IfcAxis2PlacementLinear`,
`IfcOffsetCurveByDistances`, `IfcSectionedSolidHorizontal` and
`IfcSectionedSurface`. A consumer must not evaluate curves itself, so the
relation has to be stored exactly in the graph and resolved by the
kernel.

Three things are ambiguous unless named:

1. **Which distance.** On an elevated or banked centreline the authored
   measure is plan distance; its 3D arc length runs ahead by the grade
   factor (0.125 m per 100 m at 5%). On every other curve the only
   intrinsic measure is arc length, and on a B-spline or an ellipse that
   needs quadrature and a root find (#239).
2. **Which frame.** The offsets are lateral, vertical and longitudinal,
   so they need three axes. A banked curve carries its own roll (ADR
   0081); a plain 3D curve does not, and the curve-evaluation contract's
   `frame_at` documentation disagrees with its reference provider about
   which axis is up (#242).
3. **What lies between stations.** A station run defines values at
   stations only.

## Decision

We will add `GeometryNode::CurveStation`, `CurveRelation::OffsetByStations`,
`SolidOperation::StationedSpine` and `SurfaceRelation::SectionedSurface`
to `axiolid-model`, all sharing one `Station { distance, offsets }`, and
resolve them in `axiolid-evaluate` (`station`) and `axiolid-mesh-compile`
(`station`).

- **Distance follows the basis curve.** Plan distance on `Elevated3` and
  `Banked3` (their parameter, the contract's `PlanDistance`), arc length
  on every other curve, from the curve's start (parameter 0, or a
  polyline's or B-spline's first domain parameter). A distance must lie
  in `[0, L]`; a negative or non-finite one is refused when the node is
  pushed, one beyond `L` (by more than `ARC_LENGTH_TOLERANCE * max(1, L)`)
  when it is resolved. Numerical arc length is resolved by #239's
  `arc_parameter`, to the same tolerance.
- **The frame is the basis curve's section frame:** tangent, lateral to
  the LEFT, `up = tangent x lateral`. A 2D curve lies in `z = 0` with
  `up = +Z`; a `Banked3` uses its rolled section; every other 3D curve,
  `Elevated3` included, uses the reference-up frame against `+Z` that
  `ReferenceCurveEvaluator::frame_at` returns (lateral horizontal, up
  leaning with the grade). `StationFrame::Plan` replaces it with the
  upright frame (horizontal tangent, horizontal left normal, `+Z`). A
  resolved frame is presented in the provider's layout, `x` tangent,
  `y` up, `z` right, which is the provider's actual behaviour, not the
  contract's current wording (#242).
- **Offsets** `(lateral, vertical, longitudinal)` move the point along
  `(lateral, up, tangent)`. A profile at a station maps its `x` to
  lateral and `y` to up, so its normal is the tangent.
- **Between stations everything is linear in distance:** the offsets,
  and each section point with its counterpart (a closed profile's by
  ring and vertex index, so sections must share a ring structure; an open
  section's by tag, every section carrying the same tag sequence).
- **Accuracy.** A station point is the exact curve's at a measure within
  `ARC_LENGTH_TOLERANCE * max(1, s)` where the measure is numerical, and
  at `s` itself where the parameter is the measure. Meshes between
  stations are refined by a midpoint and flatness test, a sample rather
  than a proof, so their deviation is reported `Unbounded` by name; the
  exact compiler refuses a station-placed spine by name.

## Alternatives considered

| Option | Why not |
| --- | --- |
| Arc length on every curve, also elevated ones | Contradicts how alignments are authored and the existing `PlanDistance` convention; misplaces stations by the grade factor. |
| Store an explicit convention field on each station | A second source of truth: a plan distance is meaningless on a curve with no plan, and an arc length on an elevated curve would need its own 3D quadrature for no consumer. |
| Frenet frame on plain 3D curves | Undefined on straight pieces and flips at inflections; the reference-up frame is what `frame_at` already returns. |
| Resolve stations in the source adapter | Curve evaluation is the kernel's job (the consumer's ADR 0004), and every adapter would repeat it. |
| Match closed-profile vertices by tags | `Profile` carries no vertex tags; index matching with a refused structural mismatch is exact for the parameterised and polyline profiles alignments use. |

## Consequences

**Positive**

- Distance-along-curve placements, offset curves, sectioned solids and
  sectioned surfaces lower exactly into the graph and mesh through the
  reference compiler.
- One frame rule for every consumer, consistent with `frame_at`.

**Negative / costs**

- Station-placed meshes carry no certified deviation bound.
- A basis must be an atomic curve node; a station along a curve relation
  or an instance is refused by name.
- An explicit section orientation (an `IfcAxis2PlacementLinear`'s own
  `Axis`/`RefDirection`) is not carried; only the section and plan frames
  are.

**Follow-ups / risks to watch**

- #242: once the contract text and provider agree, re-check that the
  layout named here still matches.
- A certified bound for station-placed sweeps, and exact compilation of
  the constant-section-along-a-line case.

## Relation to existing code

- `crates/representations/modeling/graph/src/station.rs`
- `crates/algorithms/parametric/evaluate/src/station.rs`
- `crates/execution/compile/src/station.rs`
- `crates/algorithms/parametric/evaluate/src/station/seam.rs` and
  `crates/representations/analytic/curve/src/seam.rs` (#263)
- ADR 0081 (banked section frame), #239 (`arc_parameter`).
- `crates/execution/compile/tests/station_placement.rs` (#264).

## Amendment 2026-10-03: explicit orientation and tag matching (#246)

The consumer (openbimrs/ifc#307) needs two things the decision above left
out: a station placement may state its own `Axis` and `RefDirection`,
which must not be parallel, and sections must be matched by tag wherever
tags are given, closed profiles included.

- **What the vectors mean.** buildingSMART IFC 4.3, `IfcAxis2PlacementLinear`:
  "Relative placement axes (Axis and RefDirection) are relative to the
  curve used for linear referencing provided in IfcPlacement Location
  (IfcPointByDistanceExpression BasisCurve), maintaining the relationship
  to the tangent of the curve"; `Axis` is "the exact direction of the
  local Z Axis", `RefDirection` "the direction used to determine the
  direction of the local X Axis ... if RefDirection is omitted, the
  direction is taken from the curve tangent", and rule WR2 forbids them
  parallel or anti-parallel. We therefore store them as COMPONENTS IN THE
  STATION'S BASE FRAME, never in world coordinates: `(a, b, c)` is
  `a * tangent + b * lateral + c * up` of the `StationFrame` the station
  names (section or plan), whose `(tangent, lateral, up)` is the
  placement's local `(X, Y, Z)` (tangent, left, up). Which base frame the
  consumer's curve frame is (grade and bank kept or dropped) is the
  consumer's `StationFrame` choice, stated per relation as before.
- **Orthonormalisation.** Gram-Schmidt with the axis primary: `up' =
  axis / |axis|` exactly, `tangent' = normalise(r - (r . up') up')` for
  the unit reference direction, `lateral' = up' x tangent'`. Defaults:
  axis `(0, 0, 1)`, reference `(1, 0, 0)`. A zero or non-finite vector,
  or a pair (given or defaulted) whose sine is at most
  `ORIENTATION_TOLERANCE = 1e-9`, is refused by name when the node is
  pushed (`axiolid-model`) and again when a frame is turned
  (`axiolid-evaluate`, `SectionFrame::oriented`).
- **What it turns.** The frame a resolved station presents and the plane
  a section's profile is placed in (profile `x` along `lateral'`, `y`
  along `up'`). The offsets stay in the base frame: they locate the
  origin, as `IfcPointByDistanceExpression` does independently of the
  placement's axes. Between two sections the unit axis and unit reference
  direction are interpolated linearly, component-wise in the base frame,
  and orthonormalised again; an interpolated pair that degenerates is
  refused by name.
- **Tags.** A run of sections is tagged throughout or not at all (mixed is
  refused by name); every section carries the same SET of tags, none
  repeated. Open sections: the tags run along the polyline in the first
  section's order or in reverse (joined reversed); anything else would
  cross the sheet and is refused. Closed sections: the tags name the
  contour vertices, outer ring then holes, each ring from its first
  segment's start in its authored sense; the profile must be a polygonal
  contour, optionally under a `Derived` transform (a mirrored guardrail).
  After winding every ring (outer counter-clockwise, holes clockwise,
  tags turning with their vertices), each ring of the first section must
  map onto one ring of every other, outer onto outer, in the same cyclic
  order; so a section may start a ring elsewhere, list its holes in
  another order or be authored with the other winding, nothing else.
  Index matching remains for untagged sections.
- **Representation, additively.** `CurveStation`, `StationedSection` and
  `StationedOpenSection` have public fields and the #241 variants public
  struct fields, so adding a field to any of them breaks construction. We
  add instead `StationOrientation`, `OrientedCurveStation` (a new
  `GeometryNode` variant) and `SectionAtStation` (`#[non_exhaustive]`,
  built by `new`/`with_tags`/`with_orientation`, so later fields stay
  additive) with two appended variants, `SolidOperation::SectionsAtStations`
  and `SurfaceRelation::OpenSectionsAtStations`, the general forms of
  `StationedSpine` and `SectionedSurface`, which compile through the same
  code. `StationFrame` could not carry the vectors: it derives `Eq` and
  `Hash`, which floats do not have.

| Option | Why not |
| --- | --- |
| Axis and reference direction in world coordinates | Contradicts the schema's reading, and a section's orientation would not follow the curve between stations. |
| Tags as a sequence, matched by index | What the #241 surface did; it cannot express a section authored from another start or in the other direction, which is exactly what tags are for. |
| Tags on flattened vertices of any profile | A chord-dependent vertex has no authored name; tags on curved or parametric profiles are refused by name instead. |
| Any permutation of tags | Lofting a ring onto a permuted ring crosses the walls; only rotations (and, for open sections, reversal) keep the section's edges. |

Still open: a station along an instance or a curve relation, and offset
curves extended past the first and last station (non-blocking for the
consumer).

## Amendment 2026-10-09: seams (#263)

The consumer (openbimrs/ifc#346) needs a stated frame where a station
lands on a tangent discontinuity of its basis, and lets the PREVIOUS
segment's tangent govern there (IFC4X3 ADD2 8.9.3.48.3); its section and
offset runs across such a discontinuity are mitred (8.8.3.35.1). The
evaluators took the next piece (`ElevationLaw::piece_at`: "seams belong
to the piece that starts there"; a polyline's spans) without saying so,
and on a polyline or a B-spline the numerical arc-length inverse could
land a hair either side of a vertex, so the frame there was not even
deterministic.

- **What a seam is.** An interior place where two pieces of the basis
  meet: a polyline's vertex, a B-spline's knot of multiplicity at least
  its degree, an arc-length chain's join or an intrinsic law's seam
  (smooth: the heading is continuous by construction), an elevated
  curve's profile seam (the grade may jump), a banked curve's cant seam
  (the roll may jump) and pivot seam (the grade of its point path may).
- **The rule.** A station within `ARC_LENGTH_TOLERANCE * max(1, s)` of a
  seam that is not smooth is ON it, and is read at the seam's own
  distance from the piece its `SeamSide` names: `Outgoing`, the piece
  that starts there, which is the default and what every evaluator
  already read, or `Incoming`, the piece that ends there, read as the
  curve truncated at the seam at its end (a polyline's previous segment,
  a B-spline's previous span through the reversed spline, a profile, cant
  and pivot law cut at the seam). At the curve's start and end there is
  one piece, and both sides read it. A consumer with a coarser precision
  snaps its distance to the seam position the kernel reports.
- **Representation, additively.** `SeamSide` (`#[non_exhaustive]`,
  default `Outgoing`) lives in `axiolid-curve`, which both `axiolid-model`
  and `axiolid-evaluate` depend on; neither may depend on the other.
  `OrientedCurveStation` is `#[non_exhaustive]`, so it gains a field,
  `seam`, set by `with_seam_side`; `CurveStation` has public fields and
  cannot, so `CurveStation::with_seam_side` returns the oriented station
  in its base frame, and a plain `CurveStation` keeps reading the
  outgoing piece. `Station`, shared by every run, is unchanged: a run
  needs no side (below). The evaluator reads a side through
  `station_section2_on` / `station_section3_on`; the side-less functions
  are the outgoing reading. A frame placement at a station (#264) reuses
  the oriented station and so its side.
- **Runs.** `OffsetByStations`, `StationedSpine`, `SectionsAtStations`,
  `SectionedSurface` and `OpenSectionsAtStations` read the outgoing piece
  at their first station and the incoming one at their last, the pieces
  they lie on. A seam strictly inside the run whose tangents (in the
  run's `StationFrame`) differ by more than `SEAM_TANGENT_TOLERANCE =
  1e-9` rad gets a section in its MITRE plane, through the seam's point
  normal to `n = (t_in + t_out) / |t_in + t_out|`, whether it falls
  between two stations or on one: each section point is placed in both
  sides' frames (offsets and orientation as usual, the longitudinal
  offset aside), each projected along its own side's tangent onto the
  plane -- where that side's extrusion of the section meets it -- and the
  midpoint taken; the two coincide when the outgoing frame is the
  incoming one turned about `t_in x t_out` (a plan turn, a grade break on
  a straight plan). The longitudinal offset moves the section along `n`.
  `|t_in + t_out| / 2`, the cosine of half the turn, at most
  `MITRE_TOLERANCE = 1e-6` is a near reversal and refused by name, as is
  a seam whose two sides do not share their point. A seam whose tangents
  agree -- a chain join, a jump in a banked curve's roll alone -- is
  sampled as before.
- **Sampling beside a mitre.** The interval on either side of a mitre is
  refined as its piece alone would be (the test reads that piece's plain
  section at the seam), then joined to the mitred section. A sampled
  midpoint section the mitre plane would cut is not inserted; an authored
  section it cuts is refused by name. The deviation of these meshes stays
  `Unbounded` by name.
- **Where the seams are.** `station_seams2` / `station_seams3` list a
  curve's seams from its stored data: distance in the station measure,
  native parameter, `smooth` (the data guarantees one frame on both
  sides) and `exact` (no quadrature: false only for a B-spline's corner
  knot). `exact_station_seams2` / `exact_station_seams3` refuse a curve
  with an inexact seam by a typed `UnsupportedInput`, as they do a chain
  whose parametric piece turns a corner inside itself.
  `axiolid_mesh_compile::station::seams` reads a graph curve: an atomic
  curve exactly, and a 3D composite or trim of lines, polylines and
  circles -- what a plain composite curve lowers to -- at the running sum
  of its pieces' lengths in the direction of travel; any other relation
  and a 2D relation are refused by name.

| Option | Why not |
| --- | --- |
| A side field on `Station` or `CurveStation` | Both have public fields: a new field breaks every struct literal; and a run's side is decided by the run. |
| Exact equality with the seam distance instead of a tolerance | The polyline's numerical arc-length inverse lands within the tolerance either side of a vertex; equality would leave the frame at a seam to rounding. |
| Mitre by rotating the incoming frame half way about `t_in x t_out` | Equal to the projection when the frames are compatible, and undefined in what to do with the roll between them when they are not (a turn on a grade); the projection lies in the plane either way. |
| Refuse a mitre whose two projections differ | Refuses every turn on a grade under the reference-up frame, the case alignments have most. |
| `SeamSide` in `axiolid-model` | `axiolid-evaluate` cannot depend on the graph; a second enum would be two sources of truth. |

Still open: seam positions of 2D relations, a roll that jumps at a seam
inside a run (sampled across as before), stations along a curve relation
(still refused), a chain's parametric piece turning a corner inside
itself (read as before), and a certified bound for a mitred run.

## Amendment 2026-10-09: nodes placed at stations (#264)

The consumer (openbimrs/ifc#311) places a curve in the frame of a station
on another curve: an `IfcSegmentedReferenceCurve`'s segments are curve
segments placed by an `IfcAxis2PlacementLinear`. `Instance` takes a
resolved `Transform3` only, so the consumer would evaluate the station
itself and bake the frame, losing the relation to the base and any claim
about its accuracy; and a station along an instance was refused (#246).

- **The relation.** `InstanceAtStation { source, station:
  OrientedCurveStation }`, a new `GeometryNode` variant appended last:
  the station-framed form of `Instance`. The station stays symbolic and
  is resolved at evaluation with its `StationFrame`, its orientation and
  its `SeamSide`, the oriented station's own rules, so a placement on a
  seam reads the piece the station names. One node places a curve, a
  solid or a surface; a placed curve is a 3D curve whatever its source's
  dimension, a placed solid or surface keeps its family, and a placed 2D
  profile is no longer a profile.
- **The axes.** The source's local `x`, `y`, `z` map onto the oriented
  frame's `tangent'`, `lateral'` (to the LEFT) and `up' = tangent' x
  lateral'`, and its local origin onto the station's point (offsets read
  in the base frame, as for every station). That is the linear
  placement's own reading quoted in the #246 amendment (local `X` the
  tangent, `Y` the left lateral, `Z` up), applied to the turned frame. It
  is a rigid, right-handed motion (`SectionFrame::placement`). It is
  deliberately NOT the section mapping (profile `x` along lateral, `y`
  along up, normal along the tangent) and not the provider layout a
  resolved station presents (`x` tangent, `y` up, `z` right): a placed
  curve is laid along the basis, a section across it. A 2D source curve
  lies in its local `z = 0`, the plane of `tangent'` and `lateral'`.
- **Exactness.** The frame is exact, rounding aside, only on a line basis
  (`station_frame_is_exact2` / `station_frame_is_exact3`), itself placed,
  if at all, in exact frames. On every other basis the distance is read
  by the arc-length inverse (an estimate) or the point by quadrature, so
  `ResolvedPlacement::exact` is false, the mesh compiler reports the
  placement `Unbounded` under `DeviationPath::StationPlacement`, and the
  exact compiler refuses it by name; on a line it places the source's
  exact B-rep rigidly.
- **Refusals.** As the oriented station's, by name: a distance past the
  basis's length, a vertical tangent where the section or plan frame
  needs a horizontal direction, a basis with no tangent (a zero-direction
  line), a degenerate orientation (refused when pushed), a basis that is
  not an atomic or placed curve.
- **Stations along a placed curve.** A station whose basis is a placed
  curve is its source's station carried by the placement: a rigid motion
  keeps arc length, and one that keeps `+Z` keeps plan distance and
  carries the source's section frame (2D, reference-up, banked) onto the
  placed curve's own, so the two readings agree. The source must be an
  atomic curve or itself a placed curve (placements compose); its seams
  are the source's. A placement whose frame tilts `+Z` (beyond `1e-12`) is
  refused by name: there the carried frame is not the placed curve's own
  reference-up frame against `+Z`, and an elevated source's plan distance
  is not the placed curve's.
- **Compilation.** A placed curve sweeps as a directrix, sampled (a 2D
  source in its own `z = 0`, then moved); its exact pieces are not read,
  so an exact swept disk along it is refused by name.

| Option | Why not |
| --- | --- |
| `Instance` with a baked transform | The consumer would evaluate curves itself, the station's relation to its base and its seam side would be lost, and an inexact frame would read as exact. |
| A `CurveRelation` variant for curves only | A second relation would be needed for solids; one node places both, as `Instance` does. |
| A reference to a station node instead of an inline station | The station would have to be pushed first for every placement; an inline `OrientedCurveStation` is the same data and keeps one node per placement. |
| Map local `x`, `y` onto lateral and up, as sections do | Contradicts the linear placement's reading (`X` along the tangent) the consumer needs; sections stand across the curve, placed curves run along it. |
| Stations along a tilted placement in the carried frame | Would silently differ from the reference-up frame every other 3D curve gets, and from the plan distance on an elevated source. |

Still open: stations along a curve relation (a composite of placed
segments included) and along a tilted placement (refused by name), exact
pieces of a placed directrix, and a certified bound for a placement on a
curved basis.
