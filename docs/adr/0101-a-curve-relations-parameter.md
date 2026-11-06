# 0101 — A curve relation's parameter

- **Status:** Accepted
- **Date:** 2026-11-06
- **Deciders:** Friedrich Schrödter
- **Supersedes:** —
- **Amends:** [0082](./0082-stations-measure-in-their-basis-curves-convention.md)

## Context

Issue #424 (consumer: openbimrs/ifc#693). A trim reads its selectors in the
parameter of the curve it trims, and an IFC sweep range on an
`IfcTrimmedCurve` circle directrix is angles of that circle. The reference
polyline compiler (#311, #317) trimmed an atomic curve only: a trim of a
trim, of a composite, of an offset or of a placed curve was refused, so a
directrix sampled as a polyline could not be compared with the swept one.

Three readers of the same relation each had a rule of their own:

- the sweep directrix (`crate::directrix`) reads a trim of a trim in its
  basis's parameter, a sweep range on a trimmed circle being angles of the
  circle, a closed conic wrapped across its seam;
- a station basis (ADR 0082, #285) read a trim of a relation as arc length
  along the relation, which differs on a trim of a trim;
- the polyline compiler had none.

## Decision

**D1. One parameter per relation.** Every reader takes the same one:

| Relation | Its parameter |
| --- | --- |
| a trim | its basis's parameter (ISO 10303-42 inheritance) |
| a composite | its 3D arc length from its start (ADR 0082) |
| an offset | its basis's parameter |
| an instance, a curve placed at a station | its source's parameter |
| an atomic curve | its own |

A trim of a relation keeps a span of that parameter. Where the carrier is a
closed conic or a trim of one the selectors are read modulo the period into
the inner trim's arc; elsewhere they must lie within the inner span, to the
linear tolerance (a station basis, which has none at hand, uses `1e-9`).
Only parameter selectors are read on a relation. The outer trim's sense
reverses the polyline the span flattens to; the span itself is the inner
trim's.

**D2. A trim of a composite is cut through the station flattening.** That is
the one arc-length measure the kernel has, so a polyline cut lands on the
point a station reads at the same distance (within the arc-length
tolerance of the inverse). The pieces between the two distances are read
from `CompositeBasis::pieces_between`; each is flattened between the
native parameters the kernel reads at its start and end measures
(`StationPiece::native_span`) with the existing certified chord bounds. A
cut inside an ellipse or a B-spline piece lies off the distance asked by the
quadrature of its length; the bound adds that certified error
(`CompositeBasis::distance_interval_at`, width at most an eighth of the
chord budget): the part of the trimmed curve a polyline does not reach,
`max(0, upper(start) - start, end - lower(end))`.

**D3. The station basis reads a trim of a trim the same way.** This amends
ADR 0082 (its amendment of 2026-11-06): a trim of a trim, through curves
placed at stations, is read in the inner trim's basis's parameter instead of
as arc length along the inner trim. A trim of a trim of a composite reads
the composite's arc length.

**D4. Refused by name.** Each is an `UnsupportedInput` for the curve
evaluation, naming the case:

- a point or arc-length selector, or a Cartesian preference, on a trim of a
  relation (the directrix, the polyline and the stations alike): #432;
- a trim of a closed composite whose selectors run against its sense: it may
  cross the closing joint like a conic's wrap, which is not decided: #433;
- a trim of a composite measured in plan distance (elevated or banked
  pieces): its parameter would be a plan distance;
- a trim of a composite that ends inside an offset segment: an offset
  segment is measured by its own length, which is not read back to its
  basis's parameter;
- an offset of a trimmed composite cut inside a piece whose length is a
  quadrature: the certified error is along the basis, not along the offset;
- a trim of a composite holding an instanced curve, a parameter curve or a
  surface curve its p-curve governs: no arc length is measured along it.

A trimmed composite needs its joints within `JOINT_TOLERANCE` (`1e-9`
relative), as a station basis does, while an untrimmed one takes the
linear tolerance: a composite whose joints are a hair apart compiles
untrimmed and is refused, as a gap, when trimmed.

## Alternatives considered

| Option | Why not |
| --- | --- |
| Measure arc length a second time inside `polyline.rs` | Duplicates the certified lengths and would disagree with the stations on a quadrature piece. |
| Read a trim of a trim as arc length along the inner trim, as stations did | A 0..90 sweep range on an `IfcTrimmedCurve` circle directrix means 0..90 degrees of the circle, as the directrix reads it. |
| Keep the station reading and refuse the nested trim in a trimmed composite | Two readings of one relation: silent wrong geometry in whichever consumer reads it the other way. The maintainer decided one convention (F1). |

## Consequences

**Positive**

- A sampled directrix, a polyline and a station on one trim of a trim agree.
- A trim of a composite compiles with a certified bound, the cut included.

**Negative / costs**

- A station on a trim of a trim moves where the outer selectors were
  distances along the inner trim (a behaviour change in `axiolid-mesh-compile`,
  recorded in its changelog).

**Follow-ups / risks to watch**

- #432 (point and arc-length selectors on a relation) and #433 (a trim across
  a closed composite's closing joint).

## Relation to existing code

- `crates/execution/compile/src/station/trim.rs`: the carrier and the span
  of a trim, shared by the polyline compiler and the station basis.
- `crates/execution/compile/src/polyline.rs`: `node`, `composite_window`.
- `crates/execution/compile/src/station/basis.rs`: `flatten_in`.
- `crates/algorithms/parametric/evaluate/src/station/composite.rs`:
  `StationPiece::native_span`, `CompositeBasis::distance_interval_at` and
  `end_points`.
- Tests `polyline_relation_trims.rs`, `station_trim_of_trim.rs`,
  `station_cut.rs`; probe `scripts/probe_polyline_relation_trims_mutants.py`.
