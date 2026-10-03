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
