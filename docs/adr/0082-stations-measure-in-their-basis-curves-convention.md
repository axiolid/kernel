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
- ADR 0081 (banked section frame), #239 (`arc_parameter`).

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
