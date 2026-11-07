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
- `crates/representations/analytic/curve/src/path.rs` (`PathBound`),
  `crates/algorithms/parametric/evaluate/src/station/path_bound.rs`,
  `crates/algorithms/parametric/evaluate/src/station/composite.rs`
  (`StationPiece::from_path_piece`),
  `crates/algorithms/parametric/evaluate/tests/curve_path_bounds.rs`,
  `crates/execution/compile/tests/station_curve_path.rs` and
  `scripts/probe_curve_path_bound_mutants.py` (#304).
- `crates/execution/compile/src/directrix.rs` (`point_parameter`),
  `crates/execution/compile/tests/point_trimmed_spline.rs` and
  `scripts/probe_point_trim_mutants.py` (#306).

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
  name; a trim of a trim is read in its basis's parameter instead, see
  the amendment of 2026-11-06); a placed curve is its source carried by the placement. The
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

## Amendment 2026-10-13: run stations at the basis's ends (#303)

An offset by distances continues its first and last offsets unchanged out
to its basis's ends (IFC 4.3 ADD2 8.9.3.42.3; consumer openbimrs/ifc#433).
A consumer whose last stated station stops short of the end must add one at
the basis's length, and for a basis whose length needs quadrature (an offset
beside a gradient curve, a spiral, a B-spline, an ellipse) a caller without
a provider cannot know it. This is #298's open path end applied to run
stations.

- **Representation, additively.** `RunStation` is `#[non_exhaustive]`, so
  it gains a field, `basis_end: Option<BasisEnd>`, with `BasisEnd { Start,
  End }` a new `#[non_exhaustive]` enum, set by the constructors
  `RunStation::at_basis_start(offsets)` and `RunStation::at_basis_end(
  offsets)`. `None`, the default, is the stated `station.distance`, as
  before. Such a station states the distance `0`, which is not read. Only
  `OffsetsAtStations` carries it; `Station`, `OffsetByStations` and the
  sections are unchanged.
- **Resolution: the evaluator's own length.** The start resolves to `0`,
  where the station measure starts on every curve. The end resolves to
  the basis's length by the kernel's own measure, exactly as #298 resolves
  an open path end: the length `station::curve_path` states for the basis,
  as `CompositeBasis` measures that path (an atomic curve's
  `StationPiece::whole`, an unbounded line's parameter domain `|d|`, a
  composite's running sum of its pieces' lengths), bitwise. A run with an
  end station therefore reads BITWISE as the same run stated at that
  number: swept as a directrix, flattened as a station basis, and handed
  out by `curve_path` (where its basis is not itself an offset, which is
  refused as an offset of an offset as before, #289). The length is read
  once, and only for a run that has an end station.
- **Order and windows after resolution.** The graph refuses, when the node
  is pushed, a start station that is not the run's first, an end station
  that is not its last, and either stating a distance other than `0`; the
  stated distances still increase strictly from a start station's `0`. The
  end's place needs the length, so the resolver refuses by name resolved
  distances that do not increase strictly ("at or past the basis's end").
  The #301 snapping then runs on the resolved distances, so an end
  station's window snaps from the resolved end onto the one seam within it,
  and the snapped distances are checked again.
- **Wire format 1.3** (ADR 0085 note): the field is omitted when `None`.

| Option | Why not |
| --- | --- |
| A flag on `OffsetsAtStations` continuing its first and last offsets to the basis's ends | A second mechanism for what an explicit end station says; it would change a variant whose fields are public, or need yet another variant. |
| A `RunPosition { At(Scalar), BasisStart, BasisEnd }` field | Two sources of truth beside `station.distance`, a public field; `Option<BasisEnd>` keeps the distance where it is and its default writes nothing. |
| A NaN or infinite `station.distance` as the marker | A non-finite distance is refused today and must stay refused, and the wire format carries no non-finite number. |
| Read an end station's stated distance as a hint, or ignore it | A number that is not read must not be stated: a distance other than `0` is refused rather than silently dropped. |
| End stations on sections (`SectionAtStation`) as well | The consumer rule is the offset's; a section run states where each section stands (`IfcSectionedSolidHorizontal` places every section explicitly). |

Still open: a run station at the basis's ends on sectioned runs
(`SectionAtStation`), and on an `OffsetsAtStations` along an offset basis
flattened as a path or a station basis (an offset of an offset, refused as
before; swept as a directrix it reads).

## Amendment 2026-10-13: curve-path piece bounds at a parameter or a point (#304)

After #298 a caller with no provider (openbimrs/ifc#436) could leave a
piece's end open, but a trim of an ellipse or a B-spline at an interior
native parameter, or at a point, starts (or ends) at a station-measure
distance that is itself a quadrature: it still could not state the piece,
nor any span between stations along a relation containing one.

- **The representation, additively.** `PathPiece` is `#[non_exhaustive]`, so
  it gains two fields, `start_bound` and `end_bound: Option<PathBound>`, and
  a new `#[non_exhaustive]` `PathBound { Parameter(Scalar), Point(Point3) }`
  (a point in the curve's own coordinates, before the piece's placement,
  `z = 0` on a 2D curve). `PathPiece::between_parameters`,
  `PathPiece::bounded`, `with_start` and `with_end` set them; the stated
  `start` or `end` beside a bound holds `-INFINITY` or `INFINITY`, not read,
  as an open end's does, so `PathPiece::length` is infinite and `is_stated`
  false until an evaluator resolves it. A stated distance keeps its meaning
  and its checks; an end both open and bounded is refused by name. Both
  bounds lie in the curve's direction, so a reversed piece is bounded as a
  forward one. `CurvePath` is not in the wire format, so no version changes.
- **Points.** The evaluator resolves a point to the one parameter naming it
  (`StationCurve::parameter_at_bound`), within the joint tolerance `1e-9 *
  max(1, |p|)`, the tolerance within which pieces meet and a compiler
  accepts a trim point: in closed form on a line, a circle or an ellipse
  (`invert2`/`invert3`, a closed conic's angle in `[0, 2 pi)`, so its seam
  is one parameter), by projection onto each segment of a polyline, and on a
  B-spline by sampling (64 per control point, and every knot) and Newton
  from each local minimum within its neighbours. Places are merged where the
  curve stays within tolerance between two parameters. None within
  tolerance is refused naming `PATH_POINT_OFF_CURVE` (never projected);
  several naming `PATH_POINT_AMBIGUOUS` (a self-crossing, a closed polyline's
  or B-spline's seam); a point on any other family (intrinsic curves, chains,
  elevated and banked curves, offsets) naming `POINT_BOUND_UNSUPPORTED`.
- **Distances: the measure `between_parameters` reads.** A parameter is
  resolved to the curve's measure there: exact where it is read from stored
  data (a line's `t |d|`, a circle's `r t`, a polyline's, a curve
  parameterised by its measure), otherwise the arc-length quadrature of an
  ellipse or a B-spline. A piece bounded at both ends IS
  `StationPiece::between_parameters` of its parameters; a closed conic's
  span runs round as a graph compiler's trim in its sense does, with its
  arithmetic (from `t0` by `t1 - t0` within one turn, else wrapped, whole
  turns apart a whole turn), so the bounded path reads BITWISE as the path
  `station::curve_path` flattens the trim into, at parameters and at points
  (checked against the compiler). One bound beside a stated or open other
  end is folded into a closed conic's first turn, an end not after the start
  lying a turn further. An offset takes no bound of its own (its parameter
  is its base's; its base piece may carry one).
- **Certification is unchanged.** `joint_intervals` and the contract's
  `path_joint_distances` bound each numerical piece's length between the
  points the kernel reads as its ends, which for a bounded piece are the
  curve's at the resolved measures (the parameters up to the rounding of the
  arc-length inverse); the tests check that they contain an adaptive
  quadrature of the curves written out independently.
- **The contract** documents bounds in the path queries, names the three
  refusals and `PATH_PIECE_BOUND_UNSUPPORTED` for a provider that reads
  stated spans only; the conformance suite reads a line from a parameter to
  a point and a circle from it round its seam to a parameter, forwards and
  reversed, certifies its joint, and requires an off-curve point and a closed
  polyline's seam to be refused.

| Option | Why not |
| --- | --- |
| A `SpanEnd { Distance, Parameter, Point, CurveEnd }` replacing `start`/`end` | Changing the type of public fields breaks semver (#298's reasoning). |
| Fold `open_end` into the bound | `open_end` is released; a `CurveEnd` bound would be a second way to say it. |
| Project an off-curve point onto the curve | A guessed bound; the rules forbid substitute geometry, and a compiler refuses such a trim point too. |
| Take the first of several parameters | Which one is a guess; the caller states the parameter it means. |
| A point's tolerance carried by the bound | A second tolerance policy beside the joint tolerance the compiler already applies to trim points; addable later as a variant. |
| `between_parameters` semantics (refuse `t1 <= t0`) on a closed conic | A path would not read as the compiler's trim of the same conic, at points across the seam in particular. |

Still open: a point bound on a curve with no inversion here (it is refused,
not iterated), a tolerance per point, and a point on an offset.

## Amendment 2026-10-13: graph trims by points on a B-spline (#306)

After #304 a `CurvePath` piece bounded by a point on a B-spline resolves,
but the reference compiler refused the same trim stated in the graph: a
trim's point selector was inverted by `invert2`/`invert3`, which have no
B-spline arm, so lowering and derivation disagreed (openbimrs/ifc#436).

- **One inversion.** A graph trim's point selector on a 2D or 3D B-spline
  (`TrimSelector::Point2`, lifted to `z = 0`, or `Point3`) is resolved by
  the inversion a point bound takes, `StationCurve::parameter_at_bound`:
  the same point-bound tolerance `1e-9 * max(1, |p|)` (the point's own,
  not the trim's), the same sampling and Newton refinement, the same
  merging of places. The compiler's `directrix::point_parameter` is the
  one place that dispatches it, for every trim it reads: a station basis
  (2D and 3D), a sweep directrix and its pieces. So a trim of a B-spline by
  points flattens to exactly the path of the piece bounded by those
  points, bitwise, with and against its sense.
- **Refusals by the #304 names.** A point off the B-spline or one it
  passes through more than once is refused naming `PATH_POINT_OFF_CURVE`
  or `PATH_POINT_AMBIGUOUS`, when the trim has no other selector to read;
  a trim preferring its points falls back to a parameter selector beside
  a refused point, as it does on every other curve, and keeps the
  generic refusal for those curves' failed inversions.
- **`invert2`/`invert3` keep their contract.** They stay closed-form only
  (exact or refused, no iteration, the tolerance their caller passes);
  the B-spline arm is the compiler's, delegating to the evaluator's
  station inversion, rather than a second iterative path inside them that
  every caller of the closed-form inversion (surface sections, NURBS
  operations) would inherit.

| Option | Why not |
| --- | --- |
| A B-spline arm in `invert2`/`invert3` | Their callers rely on closed-form exactness to their own tolerance; the point-bound tolerance is fixed by the joint rule, so the arm would ignore the tolerance argument it is given. |
| `locate2`/`locate3` (the evaluator's other B-spline point search) | A second search: no ambiguity refusal, no merging, its own sampling, so a trim and a bound could disagree. The probe's mutant replacing the shared call with it is killed. |
| The compiler's scale-of-all-selectors tolerance for B-splines | A tolerance the bound does not use; the two paths would refuse different points. |

Still open: a polyline's point selector, which the shared inversion also
resolves, is still inverted by `invert2`/`invert3` and so refused in a
graph trim; and a point selector on a trim of a relation.

Evidence: `crates/execution/compile/src/directrix.rs`
(`point_parameter`, `first_point`),
`crates/execution/compile/src/station/basis.rs` (`parameter2`),
`crates/execution/compile/tests/point_trimmed_spline.rs` and
`scripts/probe_point_trim_mutants.py`.

## Amendment 2026-10-13: certified bounds and exact bodies for runs of sections (#247)

Runs of sections at stations (`StationedSpine`, `SectionsAtStations`,
`SectionedSurface`, `OpenSectionsAtStations`) were meshed by bisection
until a midpoint and a flatness test passed, and reported
`DeviationBound::Unbounded("stationed sections: ... sampled, not
bounded")`; the exact compiler refused them. A consumer (through
ifc-geometry's IFC4X3 sectioned solids) could not cite even a deck along a
straight polyline whose rectangle widens from 2 m to 4 m.

- **The exact surface.** Between two stations vertex `k` stands at `A_k(s)
  = O(s) + lat L + vert U + long T + x_k L_t + y_k U_t`, the offsets in the
  base frame and the profile point in the frame turned by the interpolated
  orientation, all linear in distance (the meshing rule of #241 and #246).
  At a fixed `s` the placement is rigid, so the exact section's edge is the
  segment between its vertices: the surface is `P(s, v) = (1 - v) A_k(s) +
  v A_l(s)` per edge. A curved profile is the exact profile only where every
  section is that one profile (the interpolation between a profile and
  itself is the profile); different curved profiles have no vertex
  correspondence and stay unbounded by name.
- **The bound per cell** (`axiolid-mesh-compile` `station::bound`). With
  `l_k` the chord of `A_k` between the mesh's vertices and `B` the bilinear
  patch over a wall quad, `|P - B| <= max(e_k, e_l)` (`P` and `B` are the
  same convex combination of `A` and `l`) and `B` is within
  `loft::quad_deviation` of the two triangles. `e_k` is `h^2 / 8 sup |A_k''|`
  per piece, the cell split at every seam strictly inside it (a joint whose
  tangents agree, a change of curvature or grade, a jump of a banked
  curve's roll), plus the largest distance of a breakpoint value `A_k(c+-)`
  (and of an end read from another side than the cell's, #301) from the
  chord. `sup |A''|` is assembled from the frame law's derivative bounds
  (`axiolid-evaluate` `station::piece_frame_bounds`): `T'` and `T''` from the
  curve (`1/r` and `1/r^2` on a circle; `D_2 / sigma` and `2 D_3 / sigma + 6
  D_2^2 / sigma^2` on an elevated or banked curve, #252's), the axes through
  the unit-vector rule (`|u'| <= |w'| / rho`, `|u''| <= 2 |w''| / rho + 6
  |w'|^2 / rho^2`, `rho` a certified lower bound of the horizontal tangent)
  for the reference-up, plan and rolled frames, the roll `asin(g)` bounded
  from the cant law's certified suprema under both conventions. A turned
  axis adds its components' derivatives (orientations interpolated
  linearly and orthonormalised, bounded by the same rule). Where the
  distance is read by the arc-length inverse (a circle, a polyline: a
  speed constant per panel, so the rule is exact up to rounding), each
  vertex stands within `delta = 2 ARC_LENGTH_TOLERANCE max(1, s)` of the
  distance asked, and `delta (S + 2 Q)` per vertex covers the shifted
  vertices, the slivers at the ends and the longer cells.
- **Refinement.** The stations are bisected until the sampling test passes
  AND the cell's certified bound fits the budget (half of it beside a
  profile's flattening in a spine), so the report meets the request
  wherever the run is bounded; the reported value is the worst cell's,
  recomputed on the mesh's own vertices.
- **Mitres.** A cell ending on a mitre's plane is exact where the piece is
  straight and the section constant across the cell (every vertex runs on
  a line through the plane, and the wall is a planar strip cut by it):
  only its `quad_deviation` counts. Any other cell beside a mitre is
  unbounded by name: the definition of the cut surface at a mitre is the
  projection of each side's section, which a changing section or a curved
  piece does not make a ruled surface.
- **Unbounded by name.** A basis framed numerically (an ellipse, a
  B-spline, an intrinsic curve, a chain, an offset; a skewed circle; an
  elevated or banked curve with another plan or an intrinsic profile), a
  tangent too steep to bound its horizontal part, a roll reaching a
  quarter turn, the mitre cells above, and different curved profiles.
- **Exact bodies** (`ReferenceExactCompiler`). A constant section along a
  line or an open polyline whose vertices lie on one line is the profile
  extruded along the tangent in the section's plane (oblique under a
  tilted orientation); along a horizontal circle (2D, or 3D in a
  horizontal plane), unturned and without a longitudinal offset, it is the
  profile revolved about the circle's axis, which lies in the section's
  plane along its up, by `(s1 - s0) / r`; a linear taper of straight-edged
  sections in one orientation along a straight basis is the planar
  polyhedron `axiolid-construct` `loft_exact` builds, each wall a plane, a
  wall whose corners are not coplanar refused by name. The frames are read
  in closed form (the line itself, the circle's angle `s / r`). Every
  other run is refused by name; the stations are read first, so a station
  beyond the basis is refused as the mesh compiler refuses it.

| Option | Why not |
| --- | --- |
| Branch and bound against an exact parameterisation (#232's tube search) | The mesh's vertices sit at known surface parameters, so the closed per-cell argument is direct, cheaper, and bounds the refinement as it runs. |
| Keep the sampling test as the only refinement and report the bound | The report would often exceed the budget it was asked for; bounding while refining makes `meets_requested` hold wherever the run is bounded. |
| Measure the run against a fitted surface | A fitted surface is not the exact one; the bound would be relative to an approximation. |
| Bound a mitre's cut on curved pieces from the curvature | The cut surface is the projection of each side's section along its tangent, not a sweep; bounding it needs its own definition, left open. |
| An exact taper whose walls are warped, as bilinear B-spline patches | A degree (1, 1) patch is exact, but the B-rep builders here make planar walls; refused by name until a consumer needs it. |

Still open: cells beside a mitre on curved pieces or with a changing
section; numerically framed bases (their frame needs a certified
arc-length inverse); exact bodies along a circle with a changing section
(not ruled by lines) and with a turned section, and warped tapers.

## Amendment 2026-10-14: varying offset laws along ellipses and B-splines, and adaptive certified panels (#307)

The #294 certificate left two gaps a real deck found (openbimrs/ifc#443:
a sectioned solid along an offset by distances beside an ellipse, its
offsets 0.5 to 0.75 between stations 2 and 20, its sections windowed).
An offset by distances along an ellipse or a B-spline was refused naming
`UNCERTIFIED_LENGTH`, so every windowed station after its start was
refused; and uniform refinement ran out of panels at a window's
certification (`1e-3` of a `1e-5` window, `1e-8`) on ordinary curves.

- **The law read through a certified length.** The law runs linearly in
  the base's measure, on an ellipse or a B-spline its arc length, a
  quadrature. The base is measured by chords first; its arc length
  `lambda(t)` from the piece's first end to a native parameter `t` is then
  known to an interval (the panels before `t` and a sub-panel to it). The
  measured curve is the offset whose law runs linearly in the base's true
  arc length between the points the kernel reads as the piece's ends,
  `v(t) = va + (vb - va) lambda(t) / Lambda`, so it passes exactly
  through both ends; a point read at the interval's middle is off it by
  at most the law's rate times `(vb - va)` times the ratio's half-width,
  an error the chain of chords carries as it carries rounding (its
  first-order parts telescope). The table is made fine enough that this
  takes a small share of the tolerance. `P''` and `P'''` take the law's
  derivatives along the base, `v' = (vb - va) |B'| / Lambda` and its next
  two, and the frame's unit vectors' third derivative, `|u'''| <= |w'''| /
  r + 9 |w'| |w''| / r^2 + 15 |w'|^3 / r^3`, which needs `|B''''|`:
  `bound::curve_fourth_bound*` (a conic's is its `|B''|` bound, a
  B-spline's the fourth hodograph through the quotient rule).
- **A higher-order deficit.** A chord panel's deficit `L - c` is the
  integral of `|g|^2 / (|P'| (1 + cos theta))`, `g` the velocity across
  the chord, which has zero mean over the panel. So `int |g|^2 <= M2^2 h^3
  / 12` (three times tighter than the turning bound `c / cos(beta)`, which
  is kept as a ceiling), and with `M3 >= |P'''|` and `G = (P'(b) - P'(a))
  / h` (within `M3 h / 4` of `P''(mid)`), `int |g|^2` is known from both
  sides to within `h^4`, with no further derivative evaluated. The panel's
  width falls as `h^4`; panels are cut where `P''` may jump.
- **Adaptive panels.** The panels sit in a queue by width and the widest
  is halved first until the widths sum to an eighth of the tolerance
  (`PANEL_SHARE`), the rest left to the points' rounding; the midpoint
  rule beside plan-parameterised bases is refined the same way. The
  budget `MAX_CERTIFIED_PANELS` is unchanged: the deck's offset certifies
  to `1e-8` at a station 30 along in about 4,400 panels, a whole turn of
  the ellipse's offset by distances (47 long) to `1e-9` in 21,000. The
  committed benchmark `tools/benchmark/benches/certified_length.rs`
  measures the cost of the deck's joints and of a windowed station 30
  along it: about 6 ms at `1e-6`, 17 ms at `1e-8`, 34 ms at `1e-9`.
- **No interval got wider.** Every interval the #294 tests certify was
  logged before and after: each is 0.26 to 0.79 times its old width, or
  unchanged where it is exact or rounding alone.

| Option | Why not |
| --- | --- |
| Refuse the varying law and certify only constant offsets | The deck's windowed stations stay refused, and the law's dependence on the base's length is boundable. |
| Raise the panel budget for the uniform refinement | A panel's width falls only as `h^3`, so the panels grow as the tolerance to the `-1/2` rather than the `-1/3`, a cost paid in every windowed station; and a uniform share wastes the budget where the bound is already tight. |
| Richardson extrapolation of chord sums | Its remainder needs a derivative bound one order higher than the deficit bound here, for the same order. |
| Read the law at the kernel's own arc-length estimate | An estimate, not a bound: a point would be off the measured curve by an unknown amount. |

Evidence: `crates/algorithms/parametric/evaluate/src/station/certify/`
(`panels.rs`, `beside.rs`), `crates/algorithms/parametric/evaluate/src/bound.rs`
(`curve_fourth_bound*`),
`crates/algorithms/parametric/evaluate/tests/station_certify_varying.rs`
(references by a Runge-Kutta integration in the native parameter),
`crates/execution/compile/tests/station_window_varying_offset.rs` and
`scripts/probe_varying_offset_certify_mutants.py`.

## Amendment 2026-10-14: offsets beside planar 3D intrinsic curves (#302)

An offset beside a `Curve3::Intrinsic` was refused naming
`UNCERTIFIED_LENGTH`, every one of them, though a stand-alone composite of
segments with a clothoid lowers to one that is a 2D intrinsic curve lifted
to a level plane. Such a curve is now certified as the 2D case is, by the
midpoint rule on the offset's speed in closed form (the #294 bullet above),
measured in 3D so that a vertical rate lengthens it.

- **Recognised from the stored data, exactly.** The torsion law is
  identically zero by `Intrinsic3::is_planar` (a structural test on the
  stored coefficients; a law zero only once evaluated, such as a frozen
  sine cancelling its mean, is not read as zero), and the start tangent
  and normal have no `z` component, compared with `0.0`, not to a
  tolerance. Then every tangent is horizontal, the reference-up frame is
  the 2D frame lifted (`lateral` the left normal, `up` `+Z`), and the
  curvature to the left is the law times the start frame's handedness in
  plan.
- **Still refused by name**: a torsion law not identically zero, and a
  planar curve in a tilted plane, whose reference-up lateral leaves the
  plane, so the offset is not the 2D case in that plane. The stations and
  joints of the planar curve itself were already exact (its stored length).

| Option | Why not |
| --- | --- |
| Test planarity numerically (sampled torsion, or the points' distance from a plane) | An estimate decides which bound applies: a nearly planar curve would be certified with a model it does not satisfy. |
| Accept a start frame tilted within a tolerance | The reference-up lateral then leaves the plane by that tilt, and the closed-form speed no longer bounds the offset's. |
| Certify planar curves in any plane | Beside a tilted plane the lateral turns out of it with the tangent, a different closed form for a case no consumer has. |

Evidence: `crates/algorithms/parametric/evaluate/src/station/certify.rs`
(`level_curvature`),
`crates/algorithms/parametric/evaluate/tests/station_certify_planar_intrinsic.rs`
(references by RK4 on the heading and chords Richardson-extrapolated) and
`scripts/probe_planar_intrinsic_certify_mutants.py`.

## Amendment 2026-10-18: certified frame laws for frames read by quadrature (#316)

After #307 the openbimrs/ifc#443 deck compiled, but its deviation stayed
`Unbounded(UNCERTIFIED_FRAMES)`: a frame on an ellipse, a B-spline or an
offset is read through estimates (the arc-length inverse, an offset's own
length table, an offset's tangent by a difference quotient), so
`piece_frame_bounds` refused it. The #247 amendment left it open.

- **The exact frame law.** The surface a run is bounded against is the
  one on the curve the length certificates measure (#294, #307: for an
  offset by distances, the law linear in the base's true arc length
  between the points the kernel reads as the piece's ends), framed by its
  own tangent at its true arc length from the piece's start. The kernel's
  estimates are measured against it, never trusted.
- **Each read, compared.** A frame read at `u` stands at a native
  parameter `t`, where the exact point and velocity are known in closed
  form to a certified error (`panels::Measured`, the reading the length
  certificate rests on). The read point's and each read axis's distances
  from the exact frame at `t` are computed and widened by that error
  (`FrameLawBounds::point_error`, `axis_error`); the exact frame's true
  arc length to `t`, bounded by chord panels over the piece measured once
  (`beside::BaseTable`), gives the measure. At a B-spline's corner knot
  the kernel reads the named side in closed form at the knot, so there the
  read axes are the exact ones. On a composite, a piece after joints
  measured by quadrature adds its start's distance from the certified
  joint (`CompositeBasis::joint_intervals`) to the measure, which also
  closes that gap for circle and polyline pieces after such joints.
- **The law over a span.** From the certified suprema `M2 >= |P''|`, `M3
  >= |P'''|` between two reads (the base's `|B''|` to `|B''''|` and the
  frame's unit vectors beside an offset, #307's) and the speed's floor
  `sigma`: `|O'| = 1`, `|T'| = |O''| <= M2 / sigma^2`, `|T''| <= (M3 /
  sigma + 3 M2^2 / sigma^2) / sigma^2 + M2^2 / sigma^4`; the axes follow
  from the tangent by the unit-vector rule as before.
- **The cell.** `axiolid-mesh-compile`'s `station::bound` adds `e_p + e_a
  (|lat| + |vert| + |long| + c (|x| + |y|))` per vertex, `c` one for the
  base frame and `sqrt 3` for a turned one, beside `delta (S + 2 Q)`. A
  run certifies its frames' quadratures to a thousandth of its budget; on
  the deck the measure is about `1e-8`, the point error `1e-10` and the
  axis error `2e-9`, against a `5e-4` budget.
- **Still unbounded by name**: an intrinsic curve, a chain, offsets beside
  them or beside elevated and banked curves, a skewed circle, a steep
  tangent, a curve read by quadrature placed in a frame that tilts `+Z`,
  a placed offset piece, and pieces after joints whose distances cannot
  be certified. The deck's own run crosses the corner its widening law
  makes at station 2 (its slope jumps from 0 to `0.25 / 18`), where it is
  mitred on a curved piece: that cell stays `MITRED_SECTIONS`, the mitre
  rule of the #247 amendment, which this amendment does not change. Runs
  along the deck's offset that do not cross it are certified.

| Option | Why not |
| --- | --- |
| Bound the estimates themselves (the arc-length inverse's and the difference quotient's errors from higher derivatives) | The difference quotient's truncation needs the displacement's fifth derivative, the base's sixth, and every estimate's error would need its own certificate; comparing each read with the exact frame needs none of them and is tight to rounding. |
| Certify the arc length to every read from the curve's start | The cost grows with the distance and the number of cells; one table of chord panels per piece answers every read with one sub-panel. |
| Reuse the kernel's own frame as the exact one | An estimate, not the exact law: the run would be bounded against a surface no one can state. |

Evidence: `crates/algorithms/parametric/evaluate/src/station/certify/reading.rs`,
`crates/algorithms/parametric/evaluate/src/station/frame_bound.rs`
(`FrameLaw`), `crates/execution/compile/src/station/bound.rs`,
`crates/algorithms/parametric/evaluate/tests/station_frame_law.rs`,
`crates/execution/compile/tests/station_bound_quadrature.rs` (the exact
solid framed independently, by a fine quadrature) and
`scripts/probe_quadrature_frames_mutants.py`.

## Amendment 2026-11-03: clothoid and other intrinsic plans (#361)

Context: an alignment family (rheonda/cad#139) needs a gradient curve whose
plan is one intrinsic curve: lines, arcs and transitions as one piecewise
law, which openbimrs/ifc lowers so. Its runs were `UNCERTIFIED_FRAMES`.

- **The frame law is closed form.** The heading is the law's integral, so
  with `|k|` and `|k'|` bounded per seam-free piece, `|T'| <= sup|k|` and
  `|T''| <= sup|k'| + sup k^2` in arc length. On an elevated curve they are
  composed with the profile in plan distance (#252's `D_2`, `D_3`). The axes
  follow by the unit-vector rule. There is no separate twist term: the
  reference-up frame's turn about the tangent is fixed by the tangent and
  `+Z`, and the plan frame turns with the heading.
- **Only the point is a quadrature.** Each read point is compared with a
  certified `IntrinsicPositionTable`: the composite two-point
  Gauss-Legendre rule on the closed-form tangent, remainder `w^5 / 4320
  |f''''|` with `|f''''| <= |k'''| + 4|k||k''| + 3k'^2 + 6k^2|k'| + k^4`
  (the Peano kernel of the rule keeps one sign, so the vector form needs no
  `sqrt 2`), panels cut at the law's seams, checkpointed every 64 panels so
  that any read costs at most 64 panels and does not depend on the order of
  reads. Each read axis is compared with the closed-form axes. The distances
  plus the table's error and the heading's rounding are `point_error` and
  `axis_error`. The measure is exact: plan distance is the parameter.
- **Families.** A 2D intrinsic curve; a level planar `Intrinsic3`; an
  elevated curve over an intrinsic plan with a closed-form profile; every
  curvature law with closed-form bounds of `k` to `k'''` (clothoid, Bloss,
  sine and cosine transitions), with an orthonormal start frame.
- **The soundness rule.** `piece_frame_bounds` reports `point_error = 0`, so
  it keeps returning `None` for every intrinsic curve; only the reader, which
  compares what it reads, certifies them.
- **Still unbounded by name.** Chains (cubic pieces), banked curves over
  intrinsic or chain plans (cant, out of scope), intrinsic space curves with
  torsion or tilted, vertical clothoids (`ElevationLaw::Intrinsic`), and
  offsets beside intrinsic curves.
- **Measured.** On a 5 km alignment at a position error of `1e-6` (a
  thousandth of a 1 mm budget), 512 points read in bisecting order cost 198 ms
  with the midpoint positions and 3.9 ms with the table, table build
  included; a run of sections along the gradient curve compiles and
  certifies in 367 ms at a bound of 0.49 mm
  (`tools/benchmark/benches/station_clothoid.rs`).

| Option | Why not |
| --- | --- |
| Certify the kernel's 8-point Gauss read itself | Its remainder needs the 16th derivative per panel, and the read is used everywhere; comparing needs no certificate of the estimate. |
| Reuse #317's midpoint positions | `O(h^2)` makes a 10 km alignment at a thousandth of a millimetre budget need about a million panels, its answers depend on the order of reads (the refinement bisects, so reads go backwards and restart the march), and it measured fifty times slower here. |
| Fresnel functions per clothoid piece | Degree-one laws only, and a certified Fresnel evaluation is its own project. |

Evidence: `crates/algorithms/parametric/evaluate/src/bound/position_table.rs`
(`IntrinsicPositionTable`),
`crates/algorithms/parametric/evaluate/src/station/certify/intrinsic.rs`,
`crates/algorithms/parametric/evaluate/src/station/frame_bound.rs`
(`FrameLaw`),
`crates/algorithms/parametric/evaluate/tests/intrinsic_position_table.rs`
and `station_frame_law_intrinsic.rs` (points from the Fresnel series, closed
forms and a Simpson reference, never the kernel's quadrature),
`crates/execution/compile/tests/station_bound_clothoid.rs` and
`scripts/probe_clothoid_frames_mutants.py`.

## Amendment 2026-11-03: sections that cannot fold (#361)

Context: a station-placed section stands in the plane normal to its frame's
tangent. A vertex on the inner side of a bend at or past its centre of
curvature has the planes stop advancing past it: the sections cross, and the
loft is a self-intersecting shell that was reported certified. A swept disk
as wide as its bend is refused for the same reason (`horn_torus`).

- **The condition.** A vertex at `(a, b)` in the section frame with the
  longitudinal offset `c` moves along the plane normal `N` at `F = O'.N + c'
  - a L.N' - b U.N'`, in station distance. `F` is affine in `(a, b)`, so the
  smallest `F` of a section stands at a vertex, and each vertex is affine in
  `s` between two stations. On a plan arc `F = 1 - a kappa`; on a gradient
  curve in the section frame `F = sigma (1 - b kappa_v)`, `kappa_v` the
  vertical curvature; the plan frame's `U = +Z` never turns toward the
  horizontal tangent, so vertical curvature never folds it; `c' <= -sigma`
  folds a section on a line.
- **The coefficients are closed form.** `FrameLaw::bend` bounds `O'.N`,
  `L.N'` and `U.N'` over a seam-free span from the jet `(O', O'')` of the
  basis (a line, a circle, a polyline's segment, an elevated or banked curve
  over those, an intrinsic curve; a reversed piece negates `O'`, a placement
  maps both) at the span's middle and the frames' Lipschitz constants from
  `FrameLawBounds`. `FrameLawBounds::forward` bounds `O'.N` from below for
  every frame with bounds, quadrature frames included.
- **Exact decision.** The compilers bisect each seam-free span of an
  unturned run with a bend: a span is clear once the lower bound of `F` at
  every vertex, less the profile's flattening bound times `|beta|` and the
  rounding, is positive; otherwise the span's middle is a witness whose
  coefficients are known to rounding, and a vertex not certified positive
  there is refused. A touch (a vertex within rounding of the centre) is a
  non-manifold solid, so it is refused like a fold. A span still undecided
  after 48 halvings lies within rounding or the profile's flattening bound of
  a fold and is refused too, by
  invalid input like the mesh compiler's other refusals of solids that do not
  exist; the bisection stops with a budget error by name after 65536 spans.
- **Conservative bound.** A turned section, or a frame read by quadrature
  (an ellipse, a B-spline, an offset), has no closed-form bend. The rate of a
  vertex along the turned plane's normal is bounded below from the forward
  floor, the orientation's components widened by what it turns within the
  span, the frames' derivative bounds and the profile's flattening bound,
  to 8 halvings. It never refuses: where it cannot clear a span the run is
  reported `Unbounded` by name (`FOLD_NOT_EXCLUDED`), like the other runs
  whose frames are not certified.
- **Local only.** The planes advance at every station distance; a run that
  meets itself far away (two legs of a loop of the basis) is not detected. A
  sectioned surface, and an open section run, are not checked here, and #407
  tracks both. A turned section or a frame read by quadrature is never
  refused; #406 asks for an exact decision there. The exact compiler, which
  cannot build what the mesh compiler cannot vouch for, refuses such a run by
  name as unsupported.

| Option | Why not |
| --- | --- |
| Refuse every run whose bound is not positive, conservative cases too | A turned section or an ellipse's frame past the centre of a tight bend is a fold only by the bound, not by the geometry; refusing on an inequality that can be loose would refuse valid solids. Unbounded by name is the honest report. |
| Check meshed cells for self-intersection | The mesh can step over a thin fold at the sampled distances, and the exact compiler has no mesh. |
| Check a swept shape's volume sign | Needs the meshed solid, and tells nothing for a touch. |
| Refuse only where the bound is certainly negative | A touch within rounding is a non-manifold solid; the refusal must include it. |

Evidence: `crates/algorithms/parametric/evaluate/src/station/frame_bound/bend.rs`
(`BendBounds`, `FrameLaw::bend`),
`crates/algorithms/parametric/evaluate/tests/station_bend.rs` (coefficients
against finite differences of the frames' tangent and closed forms),
`crates/execution/compile/src/station/bound/fold.rs`,
`crates/execution/compile/tests/station_fold.rs` (plan arc touch and past it,
a reversed arc, a vertical sag and a crest in the section and plan frames, a
fold inside a clothoid, a receding longitudinal offset, a turned section and
an ellipse reported unbounded) and `scripts/probe_station_fold_mutants.py`.

## Amendment 2026-11-04: changing curved profiles (#361)

Context: `profile_rings` flattens each section on its own, so two different
curved profiles (a disk growing, a rounded rectangle whose corners grow, an
ellipse becoming a disk) get unrelated vertex sets. The loft between them
interpolates vertices that have no exact meaning, and the run was reported
`Unbounded` ("different curved profiles").

- **Pairing.** `paired_profile_rings` lowers each profile to oriented
  segments (a rectangle, section or contour through its exact contour, a
  circle and its bore as full-turn circles, an ellipse as a full-turn
  ellipse, a derived profile as its basis with the transform applied to the
  flattened points after the joints are snapped) and matches rings and
  segments by index after orienting each ring by the exact signed area of
  its segments (a closed form for lines, circles and ellipses; a zero or
  non-finite area does not pair). A ring reversed in one profile and not in
  another, a mirrored profile against an unmirrored one, does not pair:
  reversing keeps the mirrored start point, so each vertex would match a
  point of the opposite side. Every pair of
  segments is cut into one common number of equal steps of its own
  parameter, `n = max over the profiles of ceil(|dt| sqrt(m2 / (8 e)))`,
  with `m2` the supremum of `|C''|` over the parameter (a line 0, a circle
  `r hypot(|x|, |y|)`, an ellipse `hypot(a |x|, b |y|)`, times the largest
  stretch of a derived profile).
- **The bound holds between the stations.** The interpolated exact section
  at fraction `u` is `(1 - u) A(f) + u B(f)`, a curve with second
  derivative `(1 - u) A'' + u B''`, so `h^2 m2 / 8` bounds its distance from
  the chords at `f = k / n`, and those chords are the interpolation of the
  two rings' edges. The distance a joint moved when snapped, times the
  stretch, is added, and the whole inflated by `1 + 1e-12`. This is the
  profile's share of the run's bound and the margin of the fold check.
- **Half step.** A run whose every profile is a circle starts each ring half
  a step off the axes, as the circle flattener does (#194); a mixed run does
  not.
- **What does not pair** returns `None`, and the run keeps today's
  `Unbounded` result by name: a polyline or B-spline segment (no closed-form
  bound on `|C''|` that pairs by fraction), a centre line or composite, a
  different ring or segment count, consecutive vertices within the linear
  tolerance and more than `1 << 20` vertices. The exact compiler still
  refuses a curved section that changes between stations.

| Option | Why not |
| --- | --- |
| Match the rings by arc length or nearest vertex | The blend of two curves at equal arc length is not a curve with a second-derivative bound from the two; the bound would need a model of the matching. |
| Flatten each profile alone and keep the stations' flattening bound | Between the stations the vertices pair unrelated points, so nothing bounds the interpolated sections. |
| Pair polyline and B-spline segments by fraction too | A polyline's second derivative is a delta at its vertices and a B-spline's bound needs its knots; the first pairs only by vertex count, the second is not claimed. |
| Match segments by family | A line against an arc is a valid blend (the line has no curvature); the family adds a refusal and no certainty. |

Evidence: `crates/algorithms/construction/construct/src/profile/paired.rs`,
`crates/algorithms/construction/construct/tests/paired_profile_rings.rs`
(vertices on the exact curves, a dense sampling of the interpolated section
within the deviation, a mirrored derived profile refused, a scaled, a
non-uniformly scaled and a thin half-annulus band oriented and paired, the half step,
the profiles that do not pair),
`crates/execution/compile/tests/station_bound_curved_taper.rs` (a disk along
a line and an arc, a rounded rectangle, an ellipse becoming a disk, against
a dense sampling of the exact surface) and
`scripts/probe_curved_taper_mutants.py`.

## Amendment 2026-11-03: certified placements at stations on curved bases (#405)

Context: an alignment family (rheonda/cad#139) places elements along a
gradient curve over an intrinsic plan, or a plan curve alone, and refuses an
unbounded placement. A node placed at a station (#264) was exact on a line
and `Unbounded` on every other basis; #361 certified runs of sections, not
one placement.

- **The exact placement.** The source's local `x`, `y`, `z` on the oriented
  axes of the exact frame law (#316, #361) at the station's distance,
  snapped by its window (#294) and read from its seam side (#263), its
  origin at the offsets read in the exact base frame. The kernel's frame is
  compared with it, never trusted.
- **One read, compared** (`axiolid-evaluate` `FrameLaw::error_at`). On an
  intrinsic plan, alone, level or under a closed-form profile, the read
  point is compared with a certified `IntrinsicPositionTable` and the axes
  with the closed-form heading (#361); on a circle, a polyline, and an
  elevated or banked curve over a line or a circle, the frame is read in
  closed form at its parameter and charged its rounding (`POINT_ROUNDING`
  of the largest magnitude the read works with, its point, the curve's data
  and the placement's translation, since `r - r cos a` rounds at the scale
  of `r`; `AXIS_ROUNDING` per axis): `e_p`, `e_a`.
- **The distance's resolution.** Where the distance is resolved by the
  arc-length inverse (a circle, a polyline) the read stands within `delta =
  2 ARC_LENGTH_TOLERANCE max(1, s)` of `s` (#247); on a composite, a piece
  after joints measured by quadrature adds its start's certified distance
  from the true joint (#316). The exact frame at `s` lies within `delta
  |O'|` of the read point and `delta max(|T'|, |L'|, |U'|)` of its axes,
  the bounds taken over `[s - 2 delta, s + 2 delta]` on the side read. Plan
  distance and an intrinsic curve's arc length are the parameter: `delta =
  0`. A station within `delta` of a seam it is not on could stand on the
  other side's law: unbounded by name.
- **The placement** (`axiolid-mesh-compile` `station::placement_accuracy`,
  `PlacementAccuracy`, `PlacementBound`). With offsets in the base frame and
  an orientation turning the axes by unit components (`sqrt 3 e` for a
  turned frame, `e` for the base one): `point_error = e_p + e_a (|lat| +
  |vert| + |long|)` plus the rounding of placing the offsets, `axis_error =
  c e_a` plus, for a turned frame, `1e-14 / |axis x reference|`, the
  conditioning of Gram-Schmidt on nearly parallel inputs, `e_p` and `e_a` with the resolution added. The mesh compiler
  certifies the quadratures to a thousandth of the chord budget, as for runs.
- **The body.** A local point `q` lands within `point_error + axis_error
  |q|_1` of its exact place; `|q|_1` is convex, so `R = max |v|_1` over the
  source mesh's vertices covers every point of the mesh, and every point of
  the exact source lies within the source's bound of one of them. `F =
  point_error + axis_error R` grows every contribution of the instance, as
  the placement's rounding does (#295), and is named `Certified(F)` under
  `DeviationPath::StationPlacement`; the source is re-meshed with room for
  rounding and `F` when they would carry it over the budget. The exact
  compiler still refuses every basis but a line.
- **Still unbounded by name**: frames read by quadrature (an ellipse, a
  B-spline, an offset curve), chains, banked curves over intrinsic or chain
  plans, intrinsic space curves with torsion or tilted, a stretched or
  skewed start frame, elevated curves with an intrinsic profile or another
  plan, a skewed circle, a steep tangent, a roll reaching a quarter turn, a
  basis placed at a station whose frame is not exact (a run of sections
  along it too, #408: `FrameLaw::bounds` reads the placement's exactness),
  and a station within
  its distance's resolution of a seam it is not read on.

| Option | Why not |
| --- | --- |
| A `bound: Option<PlacementBound>` field on `ResolvedPlacement` | `placement` takes no tolerance, `None` names no reason, and every internal reader of a placement would build a position table it never reads. |
| `certified_placement` refusing by name | A refusal is an error: the mesh compiler would tell "cannot certify" from "cannot resolve" by matching errors, and a line's exact frame would read as merely certified. |
| Name the frame's error only under `StationPlacement` | A node's bound is its worst contribution, not their sum: the source's bound and the frame's would be under-reported. |
| An L2 radius with `sqrt 3` | The axes' errors are bounded per axis; `|q|_1` is tight and never larger. |
| Read a circle's angle in closed form to make it exact | Changes the reading every station rests on, for a term of `2e-12 s`. |

Evidence: `crates/algorithms/parametric/evaluate/src/station/frame_bound/error.rs`
(`FrameLaw::error_at`), `crates/execution/compile/src/station/accuracy.rs`,
`crates/algorithms/parametric/evaluate/tests/station_frame_error.rs`,
`crates/execution/compile/tests/station_placement_bound.rs` (the exact placed
body from closed forms and Fresnel series) and
`scripts/probe_station_placement_bound_mutants.py`.
Measured: 1000 blocks of 4 x 2 x 3 m, one every 5 m along the gradient
curve over the 5 km alignment, compile with their report at a 1 mm chord
error in 97 ms (median of 10 samples, 93 to 106 ms, on a machine under a
load of 70), every placement certified with a frame error of 0.32 um
(`tools/benchmark/benches/station_clothoid.rs`); a measurement, not a
performance claim.

## Amendment 2026-11-06: a trim of a trim reads its basis's parameter (#424)

A trim reads its selectors in the parameter of the curve it trims
(ISO 10303-42 `trimmed_curve` inherits its basis's parameter). The
directrix reads a trim of a trim so (`crate::directrix`: a sweep range on a
trimmed circle is angles of the circle), and since #424 the polyline
compiler does. A station basis read it differently: a trim of a relation,
and a trim is one, spanned its parameter selectors as arc length along
the relation, so `trim(trim(circle, 0, pi), pi / 4, pi / 2)` was the arc between
the distances `pi / 4` and `pi / 2` along the inner trim, the angles `pi /
(4 r)` to `pi / (2 r)` of a circle of radius `r`, and not the arc a
directrix or a polyline of the same relation runs along. One relation read two ways is
silent wrong geometry in whichever consumer reads it the other way (an IFC
sweep range on an `IfcTrimmedCurve` circle is angles of the circle).

- **Decision.** A trim whose basis is, through the curves placed at
  stations, a trim is read in the inner trim's basis's parameter, as the
  directrix and the polyline read it. The outer selectors are validated
  against the inner span and narrow it: where the inner basis is a closed
  conic they are read modulo its period into the inner trim's arc, so
  selectors on either side of the seam name the same part of the arc; else
  they must lie within the inner span, to a slack of
  `1e-12 * max(1, |x|)` in parameter units, the evaluator's own rule for a
  distance, so a nested and an unnested station basis admit the same window
  (a station basis, read from the graph alone, has no linear tolerance; the
  directrix and the polyline admit the tolerance along the curve, ADR 0101). A trim of a composite stays measured as
  arc length along the composite, which is the composite's parameter, and a
  trim of a trim of a composite reads the same arc length.
- **What did not change.** Only a trim of a trim (or of a placed trim)
  reads differently; a trim of an atomic curve, of a composite and of an
  offset (refused) read as before.
- **Selectors.** Only parameter selectors are read in a relation's
  parameter. A point or arc-length selector, or a Cartesian preference, on
  a trim of a trim is refused by name; this is the open question of
  axiolid/kernel#432, which covers the directrix, the polyline and the
  stations alike.
- **Consumers.** This changes the span a station reads on a nested trim,
  and so the point of every station on one. A consumer that wrote a
  nested trim's outer selectors as distances along the inner trim names the
  basis's parameters now; a distance along it is an `ArcLength` selector,
  which #432 would read on a relation.

Evidence: `crates/execution/compile/src/station/trim.rs` (the span a trim
keeps, shared with the polyline compiler),
`crates/execution/compile/src/station/basis.rs` (`flatten_in`) and
`crates/execution/compile/tests/station_trim_of_trim.rs` (a station on a
trim of a trim, bitwise against the polyline's ends and the swept
directrix's volume).

## Amendment 2026-11-07: sheets of open sections cannot fold either (#407)

Context: the amendment on sections that cannot fold (#361) left a sectioned
surface and an open-section run "not checked here, and #407 tracks both".
Their section planes fold exactly as a spine's do, and the sheet was reported
certified (or, where a vertex touched the centre, refused only as a collapsed
triangle).

- **The same condition, margin 0.** A sheet is decided by the spine's
  decision (`SectionLaw::fold`) at the vertices of its polylines. An open
  section is its own exact polyline, so there is no flattening margin; `F`
  is affine in the section coordinates, so the minimum stands at a polyline
  vertex. Tags join the vertices before the check, so a section authored
  reversed is judged on its joined vertices. The exact decision for unturned
  sections on closed-form bend frames (a touch within rounding included)
  and the conservative, never-refusing bound for turned sections and
  frames read by quadrature are unchanged. The refusal reads "station-placed
  sections fold" for spines and sheets alike.
- **Stricter than regular.** The rule is that the planes advance past every
  section point. Refused although the sheet is regular: a sheet wholly past
  the centre of a plan arc, and a helicoid crossing its axis.
- **Local only.** A sheet that meets itself far away is not detected (#407
  stays open for it).

| Option | Why not |
| --- | --- |
| Refuse only singular sheets (`F = 0` and the in-plane velocity parallel to the section tangent) | Needs new certified in-plane frame coefficients from axiolid-evaluate, and moving to the strict rule later would be breaking, while relaxing it later is additive. |

Evidence: `crates/execution/compile/src/station.rs` (`sectioned_surface`),
`crates/execution/compile/src/station/bound/fold.rs`,
`crates/execution/compile/tests/station_fold_sheet.rs` (both kinds in both
frames: a plan arc's centre reached, touched and cleared, a receding
longitudinal offset, a sheet wholly past the centre, a ramp's axis, reversed
tags, a turned sheet reported unbounded) and
`scripts/probe_station_fold_mutants.py`.

## Amendment 2026-11-07: certified points of elevated and banked curves over intrinsic plans (#393)

Context: an elevated or banked curve over an intrinsic plan is written
against the plan's arc length, and its point is the plan's point, which the
evaluator estimates (`intrinsic_point`, budgeted from the plan's turning
variation, not from how fast the heading oscillates). Two certificates read
that estimate as exact. On a plan of curvature `0.01 sin(50 s)` the
estimate is 1.4 cm off at 96.5 m; over `[96.4375, 96.5625]` the exact curve
lies 2.8348e-3 from the chord of the evaluated ends, and `chord_bound3`
reported 1.9531e-5 (the exact-chord bound `h^2 / 8 k`), 145 times too little
(2.5e-6 against 5.7e-4 at amplitude 0.002). A 0.3 m disk along amplitude
0.002 over `[95, 98]` at 1 mm was reported `Certified(1.1957e-3)`, and the
exact tube lay 2.8435e-3 from its mesh. A chain plan with a cubic piece and
a line plan under a vertical clothoid were reported certified on points
with no error at all (`Certified(9.999e-4)`).

- **Certificates read certified points.** `ElevatedPositionTable` takes the
  plan's point from an `IntrinsicPositionTable` (ADR 0082, #361) and adds
  the closed-form height and, for a banked curve, the pivot. Its error is
  the plan table's; the height's and the pivot's own rounding is not
  charged, as for every closed-form point in the tube certificate.
  `IntrinsicPositionTable::max_error` is one uniform bound for the table.
- **`chord_bound3` charges the evaluated ends.** For an elevated or banked
  curve over an intrinsic plan the exact-chord bound plus the larger end's
  distance of `evaluate3` from the table's point and that point's error, as
  `chord_bound2` does for the plan itself (#414). The bound of a point span
  is that distance, read from the next representable span. Over a line or a
  circle it is unchanged bit for bit.
- **The tube reads the table.** The tube wall and its end disks take their
  centre from a table built at a thousandth of the budget; the derivatives
  are closed form. The exact tube lies within `e` of the table's tube
  (`e` the table's `max_error`), so the certificate runs at `target - e` and
  adds `e`. A table that cannot be built, or `e` at or above the budget,
  leaves the tube unbounded.
- **Chain plans and intrinsic profiles are unbounded by name.** A tube
  along an elevated or banked curve over a chain plan, or under an
  intrinsic profile, reports `Unbounded(UNCERTIFIED_DIRECTRIX_POINTS)`.
  Chain plans need a position table over pieces (`TODO(#438)`, which
  also owns `chord_bound3` over them); intrinsic profiles need certified
  heights (`TODO(#481)`). Both were certified before, which #361 (item 4)
  had meant to keep refused or unbounded by name.
- **Resolution of the estimate itself** is `TODO(#480)`.

| Option | Why not |
| --- | --- |
| Certify the evaluator's eight-point Gauss read | It needs the sixteenth derivative of the heading's exponential; the #361 amendment rejected it. |
| Make `evaluate3` return table points | It changes every consumer's bits and builds a table per call. |
| A per-cell error in the `certify::Patch` trait | Not needed: the uniform error is a thousandth of the budget. |
| Charge only the comparison in the tube | Reading the certified point is tighter and needs no trait change (the deviation from the wording of #393 is commented on the issue). |

Evidence: `crates/algorithms/parametric/evaluate/src/bound/elevated_positions.rs`,
`crates/algorithms/parametric/evaluate/tests/elevated_positions.rs` (points
against the Fresnel series and the closed-form profile and pivot, the fast
harmonic against compensated Simpson), `tests/intrinsic_position_table.rs`
(`max_error`), `crates/execution/compile/src/deviation/paths.rs`,
`crates/execution/compile/tests/elevated_directrix.rs` (exact centres, the
chain plan and the intrinsic profile unbounded by name) and
`scripts/probe_elevated_positions_mutants.py`.

## Amendment 2026-11-07: canted curves over intrinsic plans (#396)

Context: the amendment on intrinsic plans (#361) read a gradient curve over
a clothoid but left a banked curve over one unbounded by name: the frame law
had no tangent law for it, its stations were not charged, its placements
were unbounded and its section folds were only checked conservatively
(never refused). The roll is closed form in the cant law; only the point
needed the position table.

- **One reader.** `IntrinsicReader` gains a banked arm rather than a second
  reader. The point is the plan table's point raised by the profile's height
  and the pivot, read from the laws the kernel reads (the pieces ending at
  the station on an incoming read, `banked_ending_at`, the truncation the
  seam code already used). The tangent is the plan heading lifted by the
  grade and the pivot's rate, and the axes are the rolled section the kernel
  frames, factored out of `banked_section` as `rolled_axes` so reader and
  kernel share one formula (bitwise unchanged).
- **The tangent law** is the elevated one with `Kind::Banked`: the roll's
  derivative bounds from the cant law (`roll_bounds`), as over a line or a
  circle. The jet is `banked_derivative` and `banked_second_derivative`, so
  a fold is an exact decision, as for banked curves over a line or circle.
- **The charge.** The horizontal normal carries the heading's error
  (`2 r / flat`), the square up that and the tangent's, and each rolled
  axis the sum of the two and the roll's rounding through `asin`
  (`16 eps / cos(roll)`). A held-rail cant (derived pivot) is not read.
- **The soundness rule stays.** `piece_frame_bounds` returns `None` for
  these curves: its point error is zero, and the position is a quadrature.

| Option | Why not |
| --- | --- |
| Certify banked frames inside `piece_frame_bounds` | Its `point_error = 0` would be unsound for a point read by quadrature. |
| A separate banked reader | It would duplicate the table, the seams and the charges. |

A canted run over a chain plan stays unbounded by name (`TODO(#394)`).

Evidence: `crates/algorithms/parametric/evaluate/src/station/certify/intrinsic.rs`,
`crates/algorithms/parametric/evaluate/src/banked.rs` (`rolled_axes`),
`crates/algorithms/parametric/evaluate/tests/station_frame_law_intrinsic.rs`
(both conventions and a pivot), `tests/station_frame_error.rs` (on and off
the seams, the pivot's rate jumping),
`crates/execution/compile/tests/station_bound_clothoid.rs`,
`tests/station_placement_bound.rs`, `tests/station_fold.rs` (a fold inside a
canted transition found and refused) and
`scripts/probe_clothoid_frames_mutants.py`.

## Amendment 2026-11-07: cells beside a mitre (#335)

Context: the amendment on certified bounds (#247) bounded a cell beside a
mitre only on a straight piece with a constant section, on the argument that
the mesh quad is then the exact wall, and left every other cell unbounded by
name (`MITRED_SECTIONS`): a mitre on a curved piece, a changing section, a
seam inside the cell. That left openbimrs/ifc#443's deck, a run along an
offset that starts to widen at a corner of the offset, without a bound
although its volume is right. The straight-leg argument was also wrong where
the mesh's mitre row is not on each side's own cut.

- **One exact solid at a mitre.** Each side's surface is continued past its
  end by extruding its end section along its end tangent (the projection of
  the amendment on seams, #263) and cut by the mitre plane through the seam
  point moved along the plane's normal by the longitudinal offset there. The
  incoming side keeps the points on or behind the plane, the outgoing side
  those on or ahead of it. Where the two cuts differ (the inner side of a
  curved piece, a turn on a grade in the section frame, a longitudinal
  offset) the boundary also holds the planar faces of the plane inside one
  cut and not the other. The mesh keeps its single mitre row, the midpoint
  of the two sides' projections moved by the offset; no topology changes.
- **The bound, closed form per cell.** For a cell with the mitre at one end
  and `c` the exact section of the cell's own side at the knot, `P` its
  projection along the end tangent onto the plane, `M` the mesh's mitre row:
  the per-vertex bound of the amendment on certified bounds holds against
  the quad `(far row, c)`; the extrusion `(c, P)` is exact; a point of the
  quad at equal height differs from the quad `(far row, P)` only in the
  plane, by the kink `|h(c)| |oblique(w)| / (w . n)` (`w = c - far row`),
  second order on a curved piece and first order in a changing section's
  rate; a point up to the plane's read error past the plane is clamped back
  onto it; `(far row, P)` and the mesh quad `(far row, M)` differ by
  `D = |M - P|`; and a face of the plane lies, by the winding number of the
  two cuts, within the crossing error of the mesh's mitre row, whose edges
  are mesh edges. The frame laws' certified derivative bounds give every
  read error (the plane's point, the tangents at the knot and across it).
  Taking the cut points at equal height, not at equal distance, is what
  keeps the kink second order: the two vertices of an edge cross the plane
  at different distances, and pairing at equal distance charges an in-plane
  slide of about the section's width times the half turn.
- **Straight legs, corrected.** A cell on a straight piece with a constant
  section keeps its closed form but is no longer the mesh quad alone: the
  mesh's mitre row is not on each side's cut when the two sides' cuts
  differ, so `D`, the row's distance from the plane and the plane's read
  error enter. A turn on a grade (reported 0.0987, the exact walls stray
  0.1007) and a longitudinal offset of 0.5 (reported 3.6e-11, the exact
  walls stray 0.5) were reported certified and wrong; they now report the
  gap, 0.2236 and 0.5.
- **Transversality, refusals, refinement.** The argument needs each vertex
  to cross the plane once over the cell and the far row to lie strictly on
  the cell's side. When it does not (a corner too sharp for the piece's
  turning over the cell, or a section at the far end on the plane) the cell
  is unbounded by name, `MITRE_NOT_TRANSVERSAL`, and so is a curved or
  changing cell between two mitres, `MITRE_ON_BOTH_ENDS`; refinement splits
  both, since a shorter cell crosses the plane more transversally, and the
  halves of a cell between two mitres have one mitre each. A cell is also refined
  while its bound exceeds the budget unless its floor, the largest `D`
  (the mesh's mitre row's distance from the projection), already exceeds it:
  no splitting removes that gap. A run with a longitudinal offset whose far
  row lies between the mesher's unshifted mitre plane and the bound's shifted
  one stays unbounded (`MITRE_NOT_TRANSVERSAL`), because the two cuts
  disagree (#484). Refinement reads a mitre's knot as the
  mitred row the mesh stands, not the plain section of the piece next to it;
  the chord of a cell beside a mitre is measured from the exact section of
  its own side there.
- **How it composes.** The global fold check of #407 reads each side only up
  to its mitre plane, since the two sides overlap before the cut by design.
  The transversality test is its own bound along the fixed mitre normal and
  does not use the fold check (#406). A frame law with certified bounds that
  lands later (#393, #396) gets its mitred cells bounded with no further
  work; a roll mismatch across the seam shows up as `D`.

| Option | Why not |
| --- | --- |
| Branch and bound against the clipped exact surface (`crate::certify`, as the swept-disk mitre beside an arc, #248) | It needs the built mesh, so it cannot steer refinement, it has no inside test for general sections on the planar faces, and it runs under a work cap. |
| Pair the two vertices of an edge at equal distance | They cross the plane at different distances; the in-plane slide charged is about the section's width times the half turn, 20 to 100 times too loose in simulation. |
| Treat the mesh's mitre section as the exact one | The exact walls would then depend on how the mesh mitres. |
| Mesh each side's own cut plus the faces between them | Tight for turns on a grade, but a change of the mesh's topology; not done until a consumer needs it. |

Evidence: `crates/execution/compile/src/station/bound/mitre.rs` (the
definition, the terms and the argument in its module documentation; unit
tests of `kink` and `reach_past` against a dense supremum),
`crates/execution/compile/tests/station_bound_mitre.rs` (an oracle sampling
the exact mitred solid piece by piece, faces on the shifted plane: a turn on
a grade, a longitudinal offset of 0.5 and 0.2, a flat turn, two arcs plain
and turned, a line beside a tight arc, a widening section, a grade break over
a circular plan, a short arc between two mitres, a corner too sharp, and a
sheet), `crates/execution/compile/tests/station_bound_quadrature.rs` (the
deck, certified at 4.88e-4 within its 1e-3 request, and a B-spline's corner
knot) and `scripts/probe_station_mitre_mutants.py`. Every test asserts the
sampled exact solid stays within the reported bound.

## Amendment 2026-11-07: certified positions of arc-length chains, and tubes over chain plans (#438)

Context: the amendment on elevated and banked curves over intrinsic plans
(#393) refused a chain plan by name, because the chord bound of an elevated
curve over one under-reported the evaluated ends. The cause is in the chain:
the evaluator places each piece at the end of the one before it, reads an
intrinsic piece by quadrature and a parametric piece at the parameter of an
arc-length inverse whose `1e-12` tolerance is an estimate. On a chain of a
straight 20 m and `k = 0.01 sin 50 s` over 100 m, `chain_chord_bound` over
`[116.4375, 116.5625]` was 1.95e-5 while the exact curve lies 2.835e-3 from
the chord of the evaluated ends (the ends are 1.72e-3 and 2.84e-3 off), and
the bound of the point span at 116.5 was 0 for a point 1.4e-2 off.

- **`ChainPositionTable` reads the exact chain.** The exact chain is the
  representation's, nothing re-fitted: each piece placed at the exact end of
  the one before. An intrinsic piece's points come from an
  `IntrinsicPositionTable`, built in the piece's identity frame at its share
  of the budget. A line, circle, ellipse or B-spline piece is read at the
  parameter whose certified arc length (the station certificate's chord
  panels, `CurveLengths`) encloses the distance; the evaluator's inverse is
  only a guess, kept when its enclosure holds the distance and bisected
  otherwise, and a piece longer than its curve is refused. A read whose
  enclosure is wider than the piece's share is refused, so `max_error`,
  `max_tangent_error` and `max_curvature_error` are suprema by construction.
  The placements carry a position error `E` and a turn error `Phi` from
  piece to piece (`|R(X) p - R(X*) p*| <= |p - p*| + |X - X*| |p*|`, the
  angles add), mapped by the start frame's stretch.
- **A parametric piece may start off the end of the one before it**, by up
  to `PIECE_FRAME_TOLERANCE`; the exact chain then jumps at the join. A read
  carries a `SeamSide`: `Outgoing` is the piece that starts there (what the
  evaluator reads), `Incoming` the piece that ends there.
- **The chord bound charges the evaluated ends (#414's rule).**
  `chain_chord_bound` is the exact-chord bound plus the larger of the lower
  end's gap (read `Outgoing`) and the upper end's (read `Incoming`), each the
  evaluator's point against the table's point plus that point's error. The
  table is built at an eighth of the exact bound, or a few ulps of the
  placement if more. The harmonic chain now reports a bound within
  [2.8348e-3, 2.9e-3] over the span and [1.39e-2, 1.5e-2] at the point; a
  span that ends where the next piece starts 5e-10 off is charged 5e-10; the
  line/cubic/arc layout's bounds are unchanged to 1%. `chord_bound3` over a
  chain plan takes the same ends from `ElevatedPositionTable`, which now
  accepts a chain plan (`point3_on`, `jet3_on`, `jet_errors`), and
  `certifies_elevated` claims a chain plan again exactly where
  `chain_certifies` does.
- **`chain_certifies` narrows.** A chain with a sinusoid or a polyline
  parametric piece is no longer claimed for flattening: no certified arc
  length exists to charge its ends against, so it stays unbounded by name.
- **The tube reads the table's jets.** Over a chain plan the wall takes its
  centre, tangent and second derivative from `jet3_on` (the plan's from the
  chain table, the profile's and the pivot's closed form), and the bound
  charges `e = e0 + max(W, C)`. With `e1`, `e2` the table's tangent and
  curvature errors, `sigma` the speed floor less `e1`, `k1 = d2 / sigma`: the
  unit tangent is off by `2 e1 / sigma`, its derivative by `(2 e2 + 2 d2 et
  + k1 e1) / sigma`, the reference's projection stays `cos tau` long, and the
  wall's point and `theta`-partial are off by `r (en + eb)`, its `s`-partial
  by `e1 + r (en1 + eb1)` (`wall_jet_errors`); `W` is twice the sum over a
  cell (`h_s <= half`, `h_theta <= pi`), covering the products of two errors
  below `1e-3`, and `C = 2 r e1 / |c'|` the end disks' tilt. Over an
  intrinsic plan `e1 = e2 = 0` and `e = e0`: the amendment on intrinsic
  plans is bitwise unchanged. A tangent known to worse than a thousandth, a
  speed floor that vanishes or a piece that turns the tangent by more than
  1.5 rad leaves the tube unbounded. A corner at a join is read as the
  piece that ends (or starts) there. Intrinsic profiles stay unbounded by
  name (`TODO(#481)`).
- **Not given here.** A frame law for chain plans (stations, placements and
  folds along a chain, banked curves over one) is #394's; the station
  paths still report them unbounded by name.

| Option | Why not |
| --- | --- |
| Certify the evaluator's eight-point Gauss read and its Newton inverse | The #361 amendment rejected certifying the quadrature itself; a uniform supremum cannot be certified for the evaluator's parameter either, and the tube needs one. |
| Fold all-intrinsic chains into one piecewise `Intrinsic2` | It does not cover a cubic or an ellipse piece; kept as a test oracle (the table against `IntrinsicPositionTable` of the equivalent piecewise curve). |
| Take the evaluator's parameter unchecked | Its tolerance is an estimate; the table checks that the enclosure holds the distance and bisects otherwise. |
| A per-cell error method on the `Patch` trait | A uniform charge per wall suffices, as in the #393 amendment. |
| A rigid-motion argument for the jets | It fails on a parametric piece under a profile: the plan is read at a shifted distance while the profile is read at `s`. |

Evidence: `crates/algorithms/parametric/evaluate/src/bound/chain_positions.rs`
(the derivations in its `//!` docs, the solver's unit tests),
`crates/algorithms/parametric/evaluate/tests/chain_position_table.rs`
(points, tangents and curvature against compensated Simpson, the cubic's
binomial series and the piecewise `Intrinsic2`; the offset join; refusals),
`tests/arc_length_chain.rs` (the harmonic chain's bound against the exact
arc), `tests/elevated_positions.rs` (the same over elevated and banked
curves, the jet against the closed forms),
`crates/execution/compile/src/deviation/paths.rs` (`wall_jet_errors` against
adversarially perturbed jets), `crates/execution/compile/tests/elevated_directrix.rs`
(`a_disk_along_a_cubic_chain_plan_is_certified`) and
`scripts/probe_chain_positions_mutants.py`.
