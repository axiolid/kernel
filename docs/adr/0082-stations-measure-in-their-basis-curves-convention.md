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
- `crates/algorithms/parametric/evaluate/src/station/composite.rs` and
  `crates/execution/compile/src/station/basis.rs` (#285).
- `crates/contracts/operations/curve-evaluate/src/contract.rs` and
  `crates/algorithms/parametric/evaluate/src/provider.rs` (#286).
- `crates/representations/analytic/curve/src/path.rs`,
  `crates/execution/compile/src/station/basis.rs` (`curve_path`) and
  `crates/algorithms/parametric/evaluate/tests/curve_evaluate_path.rs`
  (#290).
- `crates/algorithms/parametric/evaluate/src/station/offset.rs`,
  `crates/execution/compile/src/station/basis.rs` (`offset_by_stations`)
  and `crates/execution/compile/tests/station_offset.rs` (#289).
- `crates/algorithms/parametric/evaluate/src/station/certify.rs`,
  `crates/algorithms/parametric/evaluate/src/station/window.rs`,
  `crates/contracts/operations/curve-evaluate/src/interval.rs`,
  `crates/algorithms/parametric/evaluate/tests/station_certify.rs`,
  `crates/execution/compile/tests/station_window.rs` and
  `scripts/probe_seam_window_mutants.py` (#294).
- `crates/representations/analytic/curve/src/path.rs` (`PathPiece::open_end`),
  `crates/algorithms/parametric/evaluate/src/station/composite.rs`
  (`StationPiece::from_path_piece`),
  `crates/algorithms/parametric/evaluate/tests/station_certify.rs`,
  `crates/algorithms/parametric/evaluate/tests/curve_evaluate_path.rs` and
  `scripts/probe_curve_path_open_end_mutants.py` (#298).
- `crates/algorithms/parametric/evaluate/src/station/composite.rs`
  (`CompositeBasis::joint_intervals_within`),
  `crates/contracts/operations/curve-evaluate/src/contract.rs`
  (`path_joint_distances_within`),
  `crates/algorithms/parametric/evaluate/tests/station_certify.rs`,
  `crates/algorithms/parametric/evaluate/tests/curve_evaluate_path.rs` and
  `scripts/probe_joint_window_mutants.py` (#305).

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

## Amendment 2026-10-09: curve relations as station bases (#285)

The consumer (openbimrs/ifc#311, and ifc#346 part 2) measures stations
along curve relations: a plain composite curve, and a segmented reference
curve whose segments are placed at stations of a base that is itself a
composite of placed segments over a composite. Stations along any
relation were refused (#246, #264).

- **The representation is unchanged.** A basis may be a `Composite`, a
  `Trimmed`, a `SurfaceCurve` whose 3D curve governs, or an
  `InstanceAtStation`, nested and placed in any combination, 2D or 3D.
  The graph already accepts them where a curve is expected; only the
  resolvers change. An `Offset`, a `ParameterCurve`, an
  `OffsetByStations` and an `Instance` stay refused by name.
- **Flattening.** A resolver reads the relation as a composite directrix
  is read: segments in order, one whose `same_sense` is false traversed
  backwards; an untrimmed line inside a composite is its parameter domain
  `[0, 1]`; a trim of an atomic curve spans the parameters its selectors
  name (a point inverted, an arc length measured, a closed conic's trim
  possibly across its parameter seam, at most one turn); a trim of a
  relation spans its parameter selectors read as arc length along it, so
  a relation measured in plan distance cannot be trimmed (refused by
  name); a placed curve is its source carried by the placement. The
  result is a list of spans of atomic curves, each between two distances
  in its curve's station measure, reversed or not, under a composed rigid
  placement (`axiolid-evaluate`'s `StationPiece`), measured as one
  `CompositeBasis`. The evaluator holds this neutral form; the compiler
  owns the graph walk.
- **Which distance.** The distance runs end to end: the lengths of the
  pieces before a point plus its distance into its own piece, each in its
  own curve's convention. Those must agree. Every piece plan-measured (an
  elevated or banked curve, trimmed, reversed, or placed by a motion that
  keeps `+Z`): plan distance. No piece plan-measured (lines, conics,
  polylines, B-splines, intrinsic curves and chains, 2D or 3D, placed or
  not): arc length. Both: refused by name, since no one distance runs
  through it. A composite of 3D segments placed at stations of an
  elevated curve is therefore measured by its OWN segments' arc length:
  a placed segment is a rigidly moved copy of its source, and its source
  (a line, a clothoid) is arc-length-measured whatever base it was placed
  on. Only a composite whose pieces are themselves elevated or banked
  curves is measured in plan distance.
- **Joints.** Consecutive pieces must meet: the end of one within
  `JOINT_TOLERANCE * max(1, |p|)` (`1e-9`, relative to the joint's
  largest coordinate) of the next one's start. Otherwise the joint is
  refused by name, as a reversed piece where the pieces meet end to end
  or start to start (a sense that was not declared), as a gap with its
  size elsewhere. Every interior joint is a seam under the #263 rule,
  never smooth (a declared `Transition` is not a guarantee in the data):
  a station within the arc-length tolerance of it reads the piece its
  `SeamSide` names at the joint's distance, the incoming piece at its
  end, the outgoing one at its start, each at its own point (pieces that
  meet within the joint tolerance are not snapped together). At the
  composite's ends, and at each piece's own ends, the piece is read from
  inside, so a trim that starts on a seam of its curve never reads what
  it trims away. Runs mitre at a joint like at any seam; a joint whose
  pieces are farther apart than the mitre's own position tolerance is
  refused there by name. `station::seams` lists the joints and each
  piece's seams inside it; a joint after an ellipse or a B-spline piece
  is not exact (its length is a quadrature) and is refused typed by the
  exact reading.
- **Reversed pieces.** A piece traversed backwards is read at `end - u`,
  its seam sides swapped, its tangent and lateral negated and its up
  kept: the reference-up and planar frames of the reversed curve exactly;
  for a banked curve the same rolled section seen from the other end.
- **Tilted placements.** #264 refused a station along a curve placed in a
  frame that tilts `+Z`. A rigid motion keeps arc length, so for an
  arc-length-measured source the station is well defined: its point and
  tangent are carried, and its frame is the placed curve's own
  reference-up frame against `+Z` (`SectionFrame::carried`), which is the
  frame every other 3D curve gets. That is what a segment placed in a
  grade-leaning section frame needs. An elevated or banked source so
  placed stays refused by name: its plan distance and its frame are not
  the placed curve's.
- **Exactness.** A frame on a composite is exact, rounding aside, only
  where every piece up to and including the one it is read on is a line
  placed, if at all, exactly: those pieces' lengths make the distance and
  the read piece makes the frame. A placement there reports `exact`; any
  other is `Unbounded` by name under `DeviationPath::StationPlacement` and
  refused by the exact compiler, as for #264.

| Option | Why not |
| --- | --- |
| Measure every composite by 3D arc length | A composite of elevated pieces is authored in plan distance; it would misplace stations by the grade factor, as the original decision says. |
| Plan distance for a composite of segments placed on an elevated base | The segments are rigid copies of arc-length-measured curves; their plan length is not stored anywhere, and a tilted segment's plan distance would need a quadrature the data does not ask for. |
| Snap the two sides of a joint onto one point | A substitute position; reading each side's own point is honest, and the mitre refuses a joint too wide to mitre. |
| Refuse a declared reversed piece (`same_sense = false`) | Legitimate data; the reversed reading is exact for every frame the kernel builds. |
| Keep refusing tilted placements | Leaves the segments of a segmented reference curve on a grade unmeasurable, although arc length and the reference-up frame are both well defined there. |

Still open: a trim of a plan-measured relation, a closed conic's trim of
more than one turn, an exact compilation of a placement on a composite
whose read piece is not a line, and the seam side on the
curve-evaluation contract's point, tangent and frame queries (#286, see
below).

## Amendment 2026-10-09: the seam side on the curve-evaluation contract (#286)

A consumer that frames a placement through `CurveEvaluator` (openbimrs/ifc#409,
IFC4X3 ADD2 8.9.3.48.3: the previous segment's tangent governs) could not
name a side: `point_at`, `tangent_at` and `frame_at` read the piece the
measure falls in, the outgoing one on a seam. The same station lowered as
geometry could read the incoming one (#263), so the two disagreed exactly
where the consumer cares.

- **The contract reads a side additively.** `point_at_on`,
  `tangent_at_on` and `frame_at_on(curve, at, SeamSide)` are new trait
  methods WITH defaults, so no provider breaks (a defaulted trait method
  is a minor change, and `cargo semver-checks` accepts it). The default
  answers `Outgoing` with the side-less method, which reads what the
  outgoing side reads exactly on a seam. Every other side is refused by a
  typed `UnsupportedInput` naming `SEAM_SIDE_UNSUPPORTED`, on and off a
  seam: the contract does not know where a provider's seams are, and an
  `Incoming` answered with the outgoing frame is the silent substitute
  this issue is about. `SeamSide` is `axiolid-curve`'s, re-exported, so
  there is one enum.
- **What a side means is the station rule.** A measure within
  `ARC_LENGTH_TOLERANCE * max(1, s)` of a seam whose frame may jump is ON
  it and is read at the seam from the named side's piece; off a seam, and
  at the curve's ends, both sides agree with the side-less answer. The
  reference provider does not restate the rule: it asks the station module
  whether the measure is on a seam (the predicate `station_section3_on`
  itself uses) and, if so, answers with `station_section3_on`, for both
  sides, so `Outgoing` a hair before a seam reads the outgoing piece at
  the seam, as a station does. Elsewhere it answers side-lessly, so its
  refusals and its domain (a line's negative distance, a circle past one
  turn) are unchanged. The side-less queries keep reading the piece the
  measure falls in, without the tolerance: changing them is a behaviour
  change no consumer asked for, and they agree with `Outgoing` on a seam
  itself, rounding aside.
- **Which measure is located.** A distance in the provider's convention,
  which on every family it measures is the station measure (plan distance
  on an elevated or banked curve, arc length on a polyline), and the
  native parameter of an elevated or banked curve, which is its plan
  distance. A polyline's or a B-spline's native parameter is not a
  station measure: `Outgoing` there is the side-less reading, `Incoming`
  is refused by name.
- **Frame and up.** On a seam the frame is the station's section in the
  provider layout (`x` tangent, `y` up, `z` right) for `+Z`; for another
  reference up, the reference-up frame of the side's point and tangent
  against it; a banked curve's rolled section only against `+Z`, refused
  otherwise, as `frame_at` refuses it.
- **#242.** The contract text now names the axes the provider has always
  returned, `x` tangent, `y` up, `z` right (ADR 0063's own wording); the
  layout this ADR resolves frames in is the documented one. No behaviour
  changed, so consumers of `frame_at` need no change.

| Option | Why not |
| --- | --- |
| A sided `CurveMeasure` (`Distance { at, side }`) | `CurveMeasure` says how a number locates a place; a side is a second question, and a provider unaware of a new variant would refuse every sided query, `Outgoing` included. |
| Default `Incoming` to the outgoing answer off a seam | The default cannot tell a seam; answering anywhere would answer wrongly on one. |
| Snap the side-less queries to seams too | A behaviour change of the existing queries for every consumer; the sided ones carry the rule. |
| Locate a polyline's native parameter as a seam | A vertex index is not a station measure; a consumer holding one maps it to a distance first, or asks side-lessly. |

Still open: a point or tangent on a seam whose side's piece is vertical
(the station frame refuses it, so the sided point and tangent do too), and
sides on 2D curves, which the contract does not evaluate.

## Amendment 2026-10-10: curve paths on the curve-evaluation contract (#290)

A consumer that frames a placement through `CurveEvaluator` (openbimrs/ifc#418:
a linear placement on a plain composite curve, a trim, or a segmented
reference curve) could only ask about one `Curve3`. Since #285 the same
station lowered as geometry is measured along the curve relation, so the
consumer had to refuse every relation basis the kernel could already lower.
The contract cannot see the graph (`axiolid-model`) or the evaluator
(`axiolid-evaluate`), which is what reads a relation.

- **A neutral value one layer down.** `axiolid-curve` gains `path`:
  `CurvePath` (pieces in order), `PathPiece` (the span `[start, end]` of a
  curve in its station measure, `reversed`, an optional rigid `placement`
  and whether it is exact) and `PathCurve` (`Two(Curve2)`,
  `Three(Curve3)`). It is #285's `StationPiece` as an owned value: it holds
  and composes (reverse a path, place it, join paths), it measures nothing.
  It owns its curves because `axiolid-curve` is a geometry data-plane
  crate, whose values carry no borrowed references so that a native backend
  can copy them across FFI (`axiolid-model`'s `native_backend_readiness`
  test enforces it).
  `axiolid-curve` is already a dependency of the contract, of
  `axiolid-evaluate` and of `axiolid-mesh-compile`, so no crate edge is
  added and no allowlist changes. `PathCurve` and `PathPiece` are
  `#[non_exhaustive]`: an offset piece (#289) is a new variant, refused by
  name by an evaluator that does not know it.
- **One reading.** `CompositeBasis::from_path` borrows the path's curves,
  checks each piece against its curve as `StationPiece::between` does, and
  measures the pieces as `CompositeBasis::new` does; `CompositeBasis::path`
  hands them back. `StationPiece` keeps its public shape (it is published)
  and converts to and from `PathPiece`; it keeps borrowing its curve, so a
  station on a graph relation, resolved again for every station of a run,
  copies no curve. The compiler's flattening is handed out as
  `axiolid_mesh_compile::station::curve_path`, its pieces converted, so a
  consumer's path is the one a station on the relation is resolved along
  (checked bitwise).
- **The contract reads a path additively.** `path_point_at_on`,
  `path_tangent_at_on` and `path_frame_at_on(path, at, SeamSide)`, the
  plain `path_point_at`, `path_tangent_at` and `path_frame_at` (the
  outgoing side), `path_distance_convention` and `path_frame_is_exact_at`
  are trait methods with defaults: the sided ones refuse by a typed
  `UnsupportedInput` naming `CURVE_PATH_UNSUPPORTED`, the plain ones
  delegate to them with `Outgoing`, the convention is `Unsupported` and no
  frame is claimed exact. A provider that does not implement paths
  therefore refuses them all by name and never reads a path as its first
  piece or a joint from the wrong side. The conformance suite reads a line
  and a quarter arc meeting at a right angle, forwards and reversed, on the
  joint and within the seam tolerance either side: each side must read its
  own piece, the plain queries the outgoing one, a distance off the path
  and a path with a gap must be refused, the convention must be arc length
  and nothing on or after the arc exact; or every path query refused as
  unsupported, with no convention reported.
- **What a path means is the composite station rule.** The distance runs
  end to end in the pieces' common convention, every interior joint is a
  seam read from the piece `SeamSide` names at its own point, a reversed
  piece is read from its end with tangent and lateral negated, a placed
  piece's frame is carried, and a frame is exact only where every piece up
  to and including the one read is an exactly placed line: all as the
  #285 amendment states, because the reference provider does not restate
  it. It builds `CompositeBasis::from_path` and answers with
  `CompositeBasis::section_on` (point, tangent, and against `+Z` the
  station section in the provider layout), so the contract and a resolved
  station agree bitwise; against another reference up the frame is the
  reference-up frame of the section's point and tangent, and a banked piece
  is refused there, as `frame_at` refuses a banked curve. A path has no
  native parameter (its pieces' parameters do not run end to end): a
  `CurveMeasure::Parameter` is refused by name.
- **Which distance.** A path is measured in the station measure, numerical
  arc length on an ellipse or a B-spline piece included, to the accuracy
  this ADR states; the per-curve `distance_convention`, which refuses a
  distance on those families, is unchanged.

| Option | Why not |
| --- | --- |
| Keep the contract per curve; a reference-side helper resolves a relation for the caller | The consumer would depend on `axiolid-evaluate` (or the compiler) to place along a relation, the dependency the contract exists to avoid, and a second provider could not take part. |
| Move `CompositeBasis` itself down to `axiolid-curve` | Its reading needs arc-length inversion, the evaluators and the seam rule, all in `axiolid-evaluate`; the representation crate would become an evaluator. |
| Turn `StationPiece` into a re-export of `PathPiece` | `StationPiece`, `StationCurve` and their inherent constructors (which measure) are published; a foreign type cannot carry them, so the change breaks semver. They convert instead. |
| Borrow the curves in the path (`&Curve3`), as `StationPiece` does | Saves a copy per path, but `axiolid-curve` is a data-plane crate whose values must cross FFI without borrowed references. |
| Resolve stations internally through an owned `CurvePath` | Each station of a run would copy the relation's curves, or rebuild and recheck the composite per query; the internal pieces keep borrowing the graph. |
| A `CurveMeasure` or `Curve3` variant for a path | `Curve3` is an atomic value and a relation is not one; a provider unaware of a new measure would refuse it anyway, and the side is a separate argument as #286 decided. |
| Plain path queries that default to refusing | A provider that implements the sided queries would also have to implement three plain ones that can only mean `Outgoing`; delegating cannot be wrong where the sided query is right. |

Still open: an offset piece (#289), a native parameter along a path, a
point or tangent on a joint whose side's piece is vertical (refused, as
the station frame refuses it), and 2D paths as such (a 2D piece lies in
`z = 0` and is answered in 3D).

## Amendment 2026-10-10: offset curves as station bases (#289)

The consumer (openbimrs/ifc#414) lowers `IfcOffsetCurve2D`,
`IfcOffsetCurve3D` and `IfcOffsetCurveByDistances` as station bases: a
kerb or rail line beside an alignment, a placement on it. Offset relations
were refused as a basis (#241, #285), and a curve path could not carry one
(#290).

- **The graph is unchanged.** `CurveRelation::Offset { basis, distance,
  reference_direction }` and `CurveRelation::OffsetByStations` are station
  bases wherever a basis may be, nested in composites and placements. The
  neutral path gains `PathCurve::Offset(Box<PathOffset>)`, appended last:
  `PathOffset { base, law }` is the offset of ONE base piece (a span of an
  atomic curve, reversed or placed, not itself an offset, with no seam of
  its curve inside it) by an `OffsetLaw`: `Planar { distance }`,
  `Directed { distance, reference_direction }` or `Linear { start, end,
  frame }`, all `#[non_exhaustive]`. `axiolid-evaluate` reads it as
  `StationCurve::Offset(StationOffset)`, also appended; `offset_pieces`
  builds them and the compiler flattens a graph offset with it.
- **What an offset is.** `C(v) = B(v) + D(v)`, `v` the base's station
  measure. A planar offset is `distance` along the base's left lateral:
  the offset curve 2D's "anti-clockwise rotation through 90 degrees from
  the tangent" (IFC4 ADD2, `IfcOffsetCurve2D.Distance`); a 3D one is
  `distance` along `normalise(V x T)` ("in the direction V x T where V is
  the fixed reference direction and T is the unit tangent",
  `IfcOffsetCurve3D`), normalised so the distance is constant; a tangent
  parallel to `V` is refused by name, as is a planar offset of a 3D curve
  (no reference direction). An offset by stations is ADR 0082's own:
  offsets linear in distance between consecutive stations, in the run's
  `StationFrame`, from its first station to its last.
- **Which distance: the offset's own length, from its start** (the base
  piece's start, the first station). That is its arc length, or, where the
  base is plan-measured (an elevated or banked curve), its own PLAN length,
  the arc length of its plan projection. An alignment is authored in plan
  distance and so is everything staked out beside it; an offset keeps its
  base's convention (`StationCurve::convention`), so an offset of an
  alignment joins plan-measured pieces and its own neighbours under #285's
  rule, and its 3D arc length (which depends on the grade) is not used. It
  is NOT the base's distance: beside a curve of curvature `k` the offset
  runs `1 - d k` times as fast, and the issue asks for its own measure.
  The frame is the offset's own section frame: the planar frame of its own
  tangent where it runs level at one height, else the reference-up frame
  of its own point and tangent; beside a banked curve, the base's rolled
  lateral made perpendicular to its own tangent.
- **Seams.** The flattening splits the base at every seam of every piece
  (a joint, a polyline vertex, a corner knot, a grade or cant seam, a chain
  join) and at every station of a by-distances law, and offsets each span
  as one piece; the pieces must meet as composite pieces must (#285). So
  every seam of the base and every break of the distance law is a joint of
  the offset, a seam under #263, never smooth. Where the base turns a
  corner or its roll jumps and the two sides' offsets do not meet, the
  offset is refused by name ("an offset across a corner"): an offset curve's
  basis must have "a well-defined tangent direction at every point"
  (`IfcOffsetCurve2D/3D`), and the mitre a run of sections uses there is a
  plane for sections, not a curve a distance runs along. Where they meet
  (the offset vanishes at the corner, or a lateral offset across a grade
  break, whose lateral is level on both sides) the corner is a seam of the
  offset. An untrimmed atomic line inside an offset is its parameter
  domain `[0, 1]`, as inside a composite.
- **Exactness.** A constant offset of a line is a line beside it: its
  frame is exact where the base is exactly placed (`frame_is_exact_at`,
  `PathCurve::is_line`), and a placement there is exact. A constant offset
  without a longitudinal part of a circle whose offset stays in planes
  normal to its axis (a 2D or level circle under a planar, section or plan
  law; any circle with `V` along its axis) is a circle of radius `r - a`,
  `a` the displacement towards the centre: its length is `(r - a) / r`
  times the base's, a closed form, so seams after it are exact; its frame
  is read through the base's arc-length inverse, as any circle's, and not
  claimed exact. Every other offset is read numerically: its velocity is
  the base's tangent (times `1 / |T_xy|` in plan distance) plus a
  five-point difference of the displacement (step `min(L / 8, 0.01)`,
  one-sided at the piece's ends), its length an adaptive 8-point
  Gauss-Legendre quadrature to `OFFSET_TOLERANCE = 1e-9` relative, inverted
  by a safeguarded Newton step inside a panel. That is an estimate, the
  difference quotient keeping the base's `1e-12` out of reach: a frame on
  it is never exact, a placement on it is `Unbounded` by name and refused
  by the exact compiler, and a seam after it is refused typed by the exact
  reading (`station::seams`).
- **Degenerate offsets are refused by name**: a circle whose offset radius
  `r - a` falls to `CUSP_TOLERANCE * r` or below ("collapses"); an offset
  whose speed in its base's direction falls to `CUSP_TOLERANCE` of the
  base's or below at a quadrature node ("cusp or turning back"; the
  adaptive quadrature concentrates nodes at the kink `|1 - d k|` has
  there); an offset lying in one horizontal plane whose polyline through
  its panel quarter points crosses itself within the piece ("crosses
  itself"). A trim of an offset relation is refused by name too: an
  offset curve takes its basis's parameter, which is not its own length,
  and a trim of a relation reads its parameter as that length (#285).
- **Curve paths.** `station::curve_path` hands an offset's pieces out as
  `PathCurve::Offset`, and the contract's `path_*` queries answer them
  through `CompositeBasis::from_path`, bitwise as a resolved station.

| Option | Why not |
| --- | --- |
| Measure an offset in its base's distance | Not the offset's own measure: beside a curve it runs at `1 - d k` of the base's rate, so a station would not lie at its own distance along the offset it names. |
| 3D arc length beside an elevated or banked curve | Contradicts the plan convention the alignment and everything beside it are authored in, depends on the grade, and would make an offset unjoinable with plan-measured pieces. |
| Refuse offsets of plan-measured bases | Leaves the consumer's main case (a line beside an alignment) unmeasurable, although its plan length is well defined. |
| Mitre an offset across a corner, as runs of sections do | The mitre is a plane sections stand in; on a curved piece the trimmed side has no closed form, and offset curves require a tangent-continuous basis. |
| Hand out a closed form as a plain line or circle piece | The path would lose the offset the consumer lowered, and a piece borrowing its curve cannot hold a curve computed on the fly; the closed form is read inside the offset piece instead. |
| Differentiate the displacement analytically | Needs every family's curvature and frame derivative (banked roll, elevated grade, chains); a difference of the small displacement is family-blind and accurate to about `1e-10`. |

Still open: an offset across a corner (mitred), an offset of an offset, a
trim of an offset by its basis's parameter, a self-crossing out of a
horizontal plane or across pieces, exact frames on an offset circle, and a
certified bound for anything placed on a numerical offset.

## Amendment 2026-10-10: certified joint distances and seam-snapping windows (#294)

The consumer (openbimrs/ifc#423) lowers stations without linking a
provider (its ADR 0004), and IFC 4.3 ADD2 8.9.3.48.3 lets the previous
segment's tangent govern a station "within precision limits" of a tangent
discontinuity: the model's `Precision`, up to a millimetre, far wider than
`ARC_LENGTH_TOLERANCE * max(1, s)`. So it snapped stations onto joints
itself, which needs the joints' distances; since #289 those follow offset
pieces whose length is a quadrature (beside an elevated or gradient curve,
a spiral, a B-spline), an estimate the consumer cannot reproduce. Derived
placements (openbimrs/ifc#353, #409, #418) go through a caller's
`CurveEvaluator` and need the joints there. Option (b) of the issue, a path
built from piece descriptions, was not taken: (c) covers the lowering
without it.

- **(a) Certified joint distances.** `CurveEvaluator::path_joint_distances(
  path, tolerance)`, a defaulted trait method (a minor change; the default
  refuses by a typed `UnsupportedInput` naming
  `PATH_JOINT_DISTANCES_UNSUPPORTED`, never handing out an estimate), answers
  each interior joint with a `DistanceInterval { lower, upper }` that
  CONTAINS its true distance and is no wider than `tolerance`: the sum of
  the true lengths of the pieces before it, each in the path's convention
  (an offset's own plan length beside a plan-measured base), between the
  points the provider reads as each piece's ends. A length read from stored
  data is exact (a line, a circle, a polyline, an intrinsic curve, a chain,
  an elevated or banked curve, a constant offset of a line or of a circle in
  planes normal to its axis), and a path of them has point intervals. The
  reference provider answers with `CompositeBasis::joint_intervals`, so a
  joint certified through the contract is the one a station is resolved
  against. The conformance suite checks a line and quarter arc (the joint
  at 10 exactly) and a widening offset of a line (`sqrt(101)`, a quadrature
  for a provider), the tolerance refusals and a path with a gap; a provider
  may refuse the query as unsupported.
- **How a numerical length is bounded** (`axiolid-evaluate`'s
  `station::certify`). Where the measured curve `P` is evaluated in closed
  form at its native parameter -- a line, a circle, an ellipse, a B-spline,
  and any offset of one, measured by arc length -- by chords and turning:
  on a panel of width `h` the chord `c` is a lower bound of the arc; with
  `M` a certified bound on `|P''|` there, every tangent lies within `beta =
  h M / (2 m)` of the midpoint's, `m = |P'(mid)| - h M / 2`, so the chord's
  projection on the midpoint tangent, and so `c`, is at least `cos(beta)` of
  the arc. `M` is the base's certified `|B''|` (`bound::curve_derivative_
  bounds*`, #232) plus the displacement's, from the unit-vector bounds
  `|u'| <= |w'| / r`, `|u''| <= |w''| / r + 3 |w'|^2 / r^2` (`r` a lower
  bound of `|w|`) for the tangent, the reference-up lateral, its up and a 3D
  offset's `V x T`. Points carry `POINT_ROUNDING = 1e-14` of the largest
  coordinate of the data; the rounding of a whole chain of chords is bounded
  once, since the first-order parts telescope along it (a few `e`, not one
  per chord). Where the base is parameterised by its plan arc length with a
  curvature law -- an intrinsic curve or a chain (2D), the plan of an
  elevated curve (measured by plan length) -- by the midpoint rule on the
  offset's speed in closed form, `f = sqrt((1 + b' - a k)^2 + (a' + b k)^2 +
  c'^2)` (`a` lateral, `b` along the plan tangent, `c` vertical, `c'` only
  beside a 2D base measured in 3D; in a section frame beside an elevated
  curve `b = longitudinal cos g - vertical sin g`, `tan g` the grade), within
  `h^3 / 24 sup |f''|`, the supremum assembled from suprema of `k, k', k''`
  and of `z'` to `z''''` by the Leibniz and chain rules, with `inf (1 + b' -
  a k) > 0` (a cusp is refused). Panels are halved until each is within its
  share of the tolerance; a budget of `MAX_CERTIFIED_PANELS` is
  `BudgetExceeded`. Uncertifiable lengths are refused by a typed
  `UnsupportedInput` naming `UNCERTIFIED_LENGTH`: an offset of a banked or a
  3D intrinsic curve, of a chain's parametric piece, beside an intrinsic
  profile in a leaning frame; an offset by distances along an ellipse or a
  B-spline (its law runs in the base's numerical arc length, so its points
  carry the inverse's estimate); a 3D offset direction that is not vertical
  beside a plan-parameterised base; a base placed in a frame that tilts
  `+Z`; a skewed frame; a tolerance the points' rounding already exceeds.
- **(c) The seam-snapping window.** `OrientedCurveStation` gains
  `seam_window` (`#[non_exhaustive]`, so additively), a distance in the
  basis's station measure set by `with_seam_window`, default `0`; the graph
  refuses a negative or non-finite one when the node is pushed. Wherever an
  oriented station is resolved -- `station::resolve`, `station::placement`
  and so every `InstanceAtStation`, compiled or used as a basis -- a station
  whose window `[d - w, d + w]` contains exactly one seam whose frame may
  jump (a composite's joint, a polyline's vertex, a B-spline's corner knot,
  a grade, cant or pivot seam; not a smooth one) is read at that seam's own
  distance, so the #263 rule takes over and its `SeamSide` names the piece;
  its point is the seam's. The seam's distance is exact where it is read
  from stored data, else a certified interval no wider than
  `WINDOW_CERTIFICATION = 1e-3` of the window (`CompositeBasis::
  snapped_distance`, `snapped_distance2`/`3`). An interval inside the window
  is within it, one wholly outside is not; two seams within one window are
  refused naming `SEAM_WINDOW_AMBIGUOUS`, a seam whose interval straddles the
  window's edge naming `SEAM_WINDOW_STRADDLED`. A window of `0` is the
  arc-length tolerance alone, every existing station's reading. The window
  is a 1.1 field of the graph's wire format (ADR 0085 note).

| Option | Why not |
| --- | --- |
| (b) A path built from piece descriptions, the kernel computing every distance | The consumer's lowering needs no joint distance once the station carries its window; the contract query covers derived placements. |
| Snap with the kernel's estimated joint distances | The estimate is what the consumer could not trust; a station a hair outside the window would read the wrong side silently. |
| A certified quadrature from derivative bounds of the speed (Gauss-Legendre with a remainder) | Needs bounds on the speed's higher derivatives, which for an offset involve the base's fourth and fifth; the chord bound needs only `|B''|`, `|B'''|`, already certified. |
| Snap to the nearest seam when two lie in the window | The side would not say which seam the station means: a guessed piece. |
| A window on `Station` or `CurveStation` | Both have public fields (#263's reasoning); a run decides its own sides. |

Still open: certified lengths of offsets beside banked and 3D intrinsic
curves, of offsets by distances along ellipses and B-splines, and of bases
placed in tilted frames (refused by name); tolerances below about `1e-7` of
the length on B-spline offsets (the chord bound's panel budget); a window on
the curve-evaluation contract's `path_*` queries themselves (a caller with
a provider snaps with `path_joint_distances`).

## Amendment 2026-10-12: open curve-path ends (#298)

A consumer with no provider (openbimrs/ifc#427) builds a `CurvePath` over an
offset beside a gradient curve, a clothoid or a B-spline. A piece is a span
`[start, end]` in its curve's station measure, for an offset its own
length, which has no closed form there; so it could not state the end, and
`path_joint_distances` (#294 (a)), which takes the path, could not tell it.
#294's option (b), a path built from piece descriptions, is not taken: an
open end is the smaller, additive change.

- **The representation.** `PathPiece` is `#[non_exhaustive]` with public
  fields, so it gains one, `open_end: bool`, set by `PathPiece::to_end(curve,
  start)` and `PathPiece::whole(curve)` (`to_end` from `0`). Its `end` then
  holds `Scalar::INFINITY` and is not read: infinite rather than NaN so that
  an open piece still equals itself, and `PathPiece::length` and
  `CurvePath::length` are infinite, not a wrong finite number, until an
  evaluator resolves it. No open start is needed: the station measure starts
  at `0` on every curve, so "from the curve's start" is the stated `start =
  0` (`PathPiece::from_start(curve, end)`). An offset's base piece may itself
  be open (the whole of a B-spline); `PathPiece::is_stated` says whether a
  piece, its base included, states every length. `CurvePath` is not in the
  graph's wire format (ADR 0085), so no format version changes.
- **Resolution: the evaluator's own length.** An open end is resolved, when
  the path is read (`StationPiece::from_path_piece`, so
  `CompositeBasis::from_path` and every `path_*` query), to where the curve
  ends by the evaluator's own length of it -- the end `StationPiece::whole`
  reads and a graph compiler's flattening states: an unbounded line's
  parameter domain end, `|d|`, as inside a composite; a length read from
  stored data exactly (a circle's turn, a polyline, an intrinsic curve or a
  chain, an elevated or banked curve's plan length, a constant offset of a
  line or of a circle); otherwise its numerical estimate, an ellipse's or a
  B-spline's arc-length quadrature to the arc-length tolerance, an offset's
  own table to `OFFSET_TOLERANCE`. The resolved piece is then checked as a
  stated one, so a start past the end or not finite is refused by name; a
  stated end keeps every check it had.
- **Distance arithmetic runs through those numbers.** Which joint a
  distance is on (within the arc-length tolerance of the nominal joint, the
  #285 rule) and how far into its piece are read against the resolved
  lengths. So an open path reads BITWISE as the stated path holding them
  (`CompositeBasis::path` hands that path back), which is the path
  `station::curve_path` hands out for the same relation and the composite a
  resolved station is read on: the consumer's provider-free path and the
  lowered graph agree exactly. Where a joint TRULY lies is
  `path_joint_distances`' certified interval: an open end's piece ends at the
  curve's own end (a B-spline's read at its last domain parameter exactly,
  an offset's at its base's end up to the rounding of its inverse), so the interval contains the true length
  of the whole curve, and the nominal joint lies within the estimate's
  tolerance of it, not necessarily inside the interval. No query is refused
  for falling within that uncertainty: a caller holding a distance authored
  against the true geometry, and needing a side, certifies the joint with
  `path_joint_distances` or, lowering, snaps by the station's seam window
  (#294 (c)).
- **The contract** says so in the path queries' documentation, and names
  `OPEN_PATH_END_UNSUPPORTED` for a provider that reads stated spans only;
  it refuses an open end by a typed `UnsupportedInput`, never reading the
  unset `end`. The conformance suite reads the whole of a widening offset of
  the whole of a line (true length `sqrt(101)`, no length stated) and a line
  to its end, forwards and reversed, against closed forms, refuses a distance
  past its end, and certifies its joint; a provider may refuse it as
  unsupported. `station::curve_path` keeps handing out stated pieces: its
  consumers may read `PathPiece::length`, which an open end would make
  infinite.

| Option | Why not |
| --- | --- |
| (b) A path built from piece descriptions, every length left to the kernel | A second, parallel description of what `CurvePath` already carries; an open end on the existing piece covers the consumer additively. |
| Fix an open end at the certified interval's midpoint | Needs a certification (and a tolerance, and a panel budget) on every path read; refuses the offsets whose length cannot be certified (an offset by distances along a B-spline) that read today; and the path would no longer read as the composite a station on the lowered relation is resolved on. |
| Refuse a distance within the certified uncertainty of a joint | The same costs, and a refusal near every numerical joint for every caller; the uncertainty (about `1e-9` of the length) is the accuracy this ADR states for numerical offsets, and a caller who needs a side has the certified interval and the window. |
| A `SpanEnd { At(Scalar), CurveEnd }` field replacing `end` | Changing the type of a public field breaks semver; a second field holding the stated end would be two sources of truth. |
| An open start as well | The measure starts at `0` on every curve; a flag for it would state nothing a number does not. |
| A NaN `end` as the marker | A stated non-finite end is refused today and must stay refused; and an open piece would not equal itself. |
| Hand out open pieces from `station::curve_path` | Changes its output for consumers that read piece lengths, for no reading gained: its stated ends are exactly what an open end resolves to. |

Still open: a path's total length on the contract (with an open last piece
a caller without the evaluator learns it from no query; `CompositeBasis::
length` has it), and an open end on a closed conic, which runs from `start`
to the end of the turn, never across its parameter seam.

## Amendment 2026-10-12: seam-snapping windows and sides on run stations (#301)

The consumer (openbimrs/ifc#429) lowers sectioned solids and surfaces
(`IfcSectionedSolidHorizontal`, `IfcSectionedSurface`) and offsets by
distances along bases whose joints need quadrature (offsets beside gradient
curves, spirals, B-splines). Without a provider it cannot compute where
those joints fall, and #294 gave the snapping to the kernel only for single
oriented stations; run stations had no window and no side, so such runs had
to be refused.

- **Representation, additively.** `Station`, `StationedSection` and
  `StationedOpenSection` have public fields and `OffsetByStations` public
  struct fields, so none of them can gain a field (#263's reasoning). The
  #246 pattern is used instead: `SectionAtStation` is `#[non_exhaustive]`
  and gains `seam: SeamSide` and `seam_window`, set by `with_seam_side` and
  `with_seam_window`; a new `#[non_exhaustive]` `RunStation { station,
  seam, seam_window }` (`new`, the same two `with_` methods, `From<Station>`)
  is the station of a new variant appended to `CurveRelation`,
  `OffsetsAtStations { basis, stations: Vec<RunStation>, frame }`, the
  general form of `OffsetByStations`. Defaults are `Outgoing` and `0`, so
  `SectionsAtStations` and `OpenSectionsAtStations` without them read as
  before, and `StationedSpine`, `SectionedSurface` and `OffsetByStations`
  stay window-free. The graph refuses a negative or non-finite window on any
  run station when the node is pushed, by the #294 detail.
- **Snapping.** Before a run is read, each station is snapped by the #294
  rule against the run's directrix, through the same
  `axiolid_reference::station::window` the oriented station uses (no second
  copy of the rule): within its window `[d - w, d + w]` of exactly one seam
  whose frame may jump, it moves to that seam's own distance, exact or
  located by its certified interval; two such seams refuse naming
  `SEAM_WINDOW_AMBIGUOUS`, a seam straddling the edge `SEAM_WINDOW_STRADDLED`;
  a window of `0` keeps the distance. The run then uses the snapped distance
  everywhere: where the section stands, the linear interpolation of
  offsets, profiles and orientations between stations, the run's span, and,
  for an offset by stations used as a basis (#289), where its law breaks.
  Snapped distances must still increase strictly: two stations snapped onto
  one seam, or one snapped past its neighbour, are refused by name ("no
  longer increase strictly").
- **Mitre interaction (#263).** A station snapped onto a seam strictly
  inside the run whose tangents differ stands on the mitre's knot: its
  section (an offset's one point) IS the mitre section, each point placed
  in both pieces' frames with its offsets and orientation and projected
  onto the bisector plane. Both frames are used, so the station's side is
  NOT read there; a section authored 0.6 mm before such a joint without a
  window is cut by the mitre plane and refused, with a window it is the
  mitre section, whatever its side. The run's first station reads the
  outgoing piece and its last the incoming one, the pieces the run lies on,
  whatever side they name, as a curve's own start and end have one piece.
  An interior station on a seam whose tangents agree but whose frame may
  jump (a banked curve's roll) is not mitred: its section or offset point
  is placed in the frame of the piece its side names (`Outgoing`, the
  default, being today's reading). So the side matters exactly where the
  run itself does not decide the frame. An offset by stations used as a
  basis does not read the side: each of its stations is a joint of the
  offset, read by whatever station is measured along it.
- **Wire format 1.2** (ADR 0085 note): the new fields are omitted at their
  defaults, and `OffsetsAtStations` is a new variant; a graph that uses
  neither keeps writing 1.0 or 1.1.

| Option | Why not |
| --- | --- |
| Fields on `Station` | Public fields: every struct literal breaks, and `StationedSection`, `CurveStation` and `OffsetByStations` would carry a window they cannot honour consistently. |
| One window per run instead of per station | A run's stations stand beside different joints, and the consumer's precision applies per station; a side per station is needed where the roll jumps. |
| Read an interior station's side at a mitre (one piece's plain section) | A section standing on a joint where the tangent turns would cross the next piece's sections; the mitre is the only section both neighbours meet. |
| Honour a side on the run's first or last station | The run lies on one piece there; the other piece's frame would place the end section on a piece the run does not cover. |
| Snap with the kernel's estimate of a numerical joint | As in #294: the estimate is what the consumer could not trust. |

Still open: windows on `StationedSpine`, `SectionedSurface` and
`OffsetByStations` themselves (their general forms carry them); a side at a
run's ends that disagrees with the run's reading is ignored, not refused.

## Amendment 2026-10-13: joints near a window on the curve-evaluation contract (#305)

`path_joint_distances` (#294 (a)) certifies every joint of a path, so one
piece whose length cannot be certified anywhere on it (an offset by
distances along a B-spline, `UNCERTIFIED_LENGTH`) refused every derived
placement on that path, however far from it (openbimrs/ifc#436). A station
with a seam window (#294 (c), #301) already certifies only the seams it
needs; a caller deriving through the contract had no equivalent.

- **The query.** `CurveEvaluator::path_joint_distances_within(path, around,
  window, tolerance)`, a defaulted trait method (a minor change), answers
  `(joint, interval)` pairs in order: `joint` the index among the path's
  joints (its position in `path_joint_distances`' answer), `interval`
  certified as that query certifies it, no wider than `tolerance`.
  Returned: every joint whose interval meets `[around - window, around +
  window]`, edges included (`DistanceInterval::meets`, the same test the
  seam-snapping window uses to tell a seam outside it). Left out: only
  joints certified to lie outside it. A centre that is not finite or a
  window negative or not finite is `InvalidInput` (`window_span`).
- **What is certified.** Every piece BEFORE the window, since its length
  moves every joint after it; the certified intervals are summed as
  `joint_intervals` sums them and widen with each numerical piece, within
  the tolerance. Each joint that can fall WITHIN the window. Nothing AFTER
  it: a joint is known to lie past the window once the certified lower end
  of the joint before it plus a certified lower bound of its piece's
  length with no quadrature passes the window's end. That bound is the
  piece's chord between the points read as its ends (the plan chord on a
  plan-measured piece, since a plan length is never shorter than it), less
  `100 * POINT_ROUNDING` of the coordinates and the length. So the reading
  stops at the first such joint, and an uncertifiable piece after the
  window, or across its end where its chord alone carries its own end past
  it (a window inside the piece, far from its end), is never certified.
- **An uncertifiable piece before the window, or within reach of it,
  refuses the query by name** (`UNCERTIFIED_LENGTH`): every joint after it
  moves with its length, so where the window lies against them is unknown.
  No coarser bound is substituted. One that reaches every such piece (the
  chord below; above, the base's length scaled by the offset's largest
  speed factor) would be wider, by a share of the piece's length, than any
  tolerance a placement asks for, so the joints after it could be returned
  only by breaking the width the contract states.
- **Agreement.** The reference provider (`CompositeBasis::
  joint_intervals_within`) shares the tolerance among ALL pieces by their
  estimated lengths, as `joint_intervals` does (no certification needed
  for that), and sums in the same order, so each joint it returns carries
  bitwise the interval `path_joint_distances` gives it wherever that
  succeeds; a joint whose interval lies inside the window is always
  returned, and one left out truly lies outside. The seam-snapping window
  keeps its own sharing among the pieces before it (#294 `perf`): a
  station's snapped distance needs no agreement with a whole-path answer.
- **The default** delegates to `path_joint_distances` and keeps the joints
  that meet the window, as the plain path queries delegate to the sided
  ones: it cannot be wrong where that query is right. A provider that does
  not certify joints therefore refuses it naming
  `PATH_JOINT_DISTANCES_UNSUPPORTED`, and one that certifies whole paths
  answers windows without change, if no more often. The conformance suite
  asks four lines turning at right angles (joints at 10, 20, 30) and the
  widening offset then line of #294 for windows on, between, over and past
  their joints: each returned joint certified and meeting its window, every
  joint truly within a window returned, each agreeing with the provider's
  whole-path answer, a bad window or tolerance refused; or the query
  refused as unsupported.

| Option | Why not |
| --- | --- |
| A per-joint query, `path_joint_distance(path, joint, tolerance)` | The caller holds a distance, not a joint index: it would have to know which joints lie near it, which is what it cannot tell without the certificates; and every joint before it would be certified again per call. |
| A coarse certified bound for an uncertifiable piece before the window | Wider than any tolerance a placement asks for, so the joints after it could not be returned anyway; deciding only that none falls in the window would need an upper bound on an offset whose law runs in a numerical arc length, a new certificate for little gain. |
| Share the tolerance among the pieces up to the window only (as the seam window does) | Fewer panels when long numerical pieces follow, but the windowed interval would differ from `path_joint_distances`' for the same joint: two certificates of one joint a caller would have to reconcile. |
| Return only the joints whose interval lies inside the window | A joint straddling the window's edge can fall within it; leaving it out would tell the caller there is none. The caller sees the interval and decides, as the seam window refuses a straddling seam by name. |
| Default refusing by name, not delegating | A provider that certifies whole paths would refuse windows it can answer correctly. |

Still open: joints after an uncertifiable piece (certified lengths of
offsets by distances along ellipses and B-splines, #294); a piece's own
seams (a polyline's vertex, a B-spline's corner knot) on the contract,
which only the seam-snapping window certifies.
