# 0092 — Mesh sections report the solid's cut as closed regions

- **Status:** Accepted
- **Date:** 2026-10-29
- **Deciders:** Friedrich Schrödter
- **Supersedes:** —
- **Amends:** [0091](./0091-mesh-sections-report-in-plane-regions-and-open-chains.md)

## Context

Issue #372 (consumer: rheonda/cad). A wall whose sill lies exactly in the
cutting plane has piers that the plane cuts across. ADR 0091 reports those as
open `SectionChain`s that end on the sill's `InPlaneRegion`, and the sill
itself as a region with a side. A consumer that wants to hatch "the wall cut
just below the plane" has to close each chain along the right region edges by
itself: it must know which region edges the cut lies across, and which
vertices are shared. That is topology the kernel already holds, and the
consumer cannot rebuild it from rounded coordinates.

The kernel also must not choose "just below" or "just above" for the
consumer, which ADR 0091 decided. The question here is what the kernel can
report without choosing.

## Decision

An opt-in request returns the **cut** of a solid as closed regions.
`SectionRequest::with_cut_regions(CutRegionLimits)` makes
`PlaneSection::cut_regions` a `Some(Vec<CutRegion>)`; without it the field is
`None`. All types are additive (`#[non_exhaustive]`), so no breaking change
and no `semver-exceptions` entry.

- **Semantics.** The cut is the closure of the plane points that have the
  solid's material on **both** sides: where the solid's interior meets the
  plane. In-plane faces are not part of it. They stay `InPlaneRegion`s with
  their side (ADR 0091). A cut just below the plane is therefore the cut
  regions together with the `Below` regions, and a cut just above is the cut
  regions together with the `Above` regions. The rings run through
  bit-identical vertices, so a consumer joins them without topology of its own
  and still makes the choice.
- **Which region edges bound the cut.** Take an edge on the boundary of an
  in-plane region and the off-plane triangle next to it. If that triangle
  leaves the plane on the side **opposite** the region's material (a `Below`
  face with the neighbour rising, an `Above` face with the neighbour falling),
  the cut lies across the edge. If it leaves on the material side, the edge
  borders the outside, as at the top of a box or a pedestal. The test is an
  exact sign comparison on the source vertices; the edge is directed with the
  cut on its left.
- **Rings.** The directed cut segments and the bounding region edges form a
  graph in which each node has as many outgoing as incoming edges. Rings are
  traced on source topology, so which edge links to which never depends on a
  coordinate. A node is dropped when both of its edges are cut segments and
  they are exactly collinear (the rule of ADR 0091), or when both are region
  edges and the region ring itself dropped that vertex. A junction of a cut
  segment and a region edge stays, which keeps chain ends bit-identical to
  region vertices. Outer rings run counter-clockwise and holes clockwise, so
  the material lies on the left, and every closed chain of a solid is a cut
  ring bit for bit. An open chain with both ends on a region is a contiguous
  run of exactly one ring, forward or reversed.
- **Exactness, amended.** ADR 0091 decides every step on source data. Two
  decisions here cannot: whether a ring is an outer ring or a hole, and how
  rings nest. An edge-crossing node has no exact binary64 position, so these
  are decided exactly **on the returned coordinates**: the sign of the
  polygon area in expansion arithmetic, and `orient2d` for crossing and
  containment tests. This is the precedent of the closed chain of a sheet
  (ADR 0091). The coordinates are the ones the consumer receives, so the
  decision is auditable on the output.
- **Nesting.** A ring is inside another when every vertex of it that is not on
  the other is strictly inside. The innermost container decides: an outer ring
  has no container or lies in a hole; a hole lies in an outer ring, which owns
  it. Candidate edge pairs of a bounding-box sweep and ring pairs of the
  containment test share one budget (`CutRegionLimits::max_pairs`). Where
  two rings touch, the directions in which one leaves the contact point are
  classified, with exact `orient2d`, against the material side of the other;
  directions on both sides mean the rings cross there and are refused, while
  neighbours that share an edge with opposite directions pass. Two
  non-adjacent edges of one ring that touch (a repeated vertex) and a spike
  are refused.
- **Refusals by name.** `REFUSED_PINCHED_CUT_REGION` where the cut boundary
  touches itself at a vertex, so the pairing at that vertex is not unique.
  `REFUSED_CROSSING_CUT_RINGS` where rings cross, overlap along an edge with
  the material on the same side, or nest inconsistently; this is how
  overlapping shells that the mesh audit cannot see are refused (a
  self-intersecting mesh is still not certified, ADR 0033). A provider that
  reports a `Both` region refuses a cut-region request with
  `REFUSED_IN_PLANE_FACES_ON_BOTH_SIDES` until the meaning of `Both` is
  decided. A sheet has no inside: a cut request with `SectionRole::Sheet` is
  `InvalidInput`. A ring that collapses under projection or has zero area is
  `Degenerate`. Budgets get their own names: cut regions, cut region
  vertices, cut pairs.
- **Registry.** A provider must answer cut regions exactly when asked: an
  answer without them, or with unrequested ones, is a contract violation. The
  rings get the checks of in-plane regions (cardinality, finiteness, adjacent
  duplicates, canonical start and order), the budgets and the evidence counts
  (`cut_regions`, `cut_region_vertices`, `cut_pairs`). A provider from before
  this feature ignores the request and answers with `cut_regions == None`; the
  registry treats that as `Unsupported` and tries the next provider, so the
  request still fails closed when none traces the cut. A sheet request is
  refused before dispatch. Orientation is not checked there (#274).
- **Conformance.** The region suite runs every fixture again with cut regions
  requested. A provider whose first such request returns `Unsupported` does
  not claim them (`cut_regions_claimed == false`) and passes. New fixtures with
  hand-derived cuts: an L-step, a wall at a sill, a wall at a head, a stepped
  slab, a pocket floor and nested tubes. The suite checks that the answer
  equals the plain one, the cut regions and their evidence, and the three
  relations between chains, regions and rings above.

## Alternatives considered

| Option | Why not |
| --- | --- |
| Close the open chains in place | Changes the chain output of ADR 0091 and mixes non-crossings into chains. ADR 0091 already rejected a closed loop around an in-plane face. |
| Return the regularised plane ∩ solid, in-plane faces included | No real cut equals it: it merges `Below` and `Above` and hides the side, which is the consumer's decision. |
| Merge each side's regions with the cut into one region per side | A planar union along shared arcs needs a boolean on the output. It can be added later as an option if wanted. |
| Always return cut regions | The new refusals (pinch, crossing rings) would reach callers that never asked. Opt-in keeps them away. |
| Decide outer or hole from segment directions | Directions say which side the material is on, not whether a ring encloses its neighbours. Nesting needs the geometry anyway. |

## Consequences

**Positive**

- A floor plan at a sill or head level gets the piers as closed regions, and
  the consumer composes "just below" or "just above" without any topology.
- Every output ring is traceable to chains and region rings by identity.

**Negative / costs**

- Outer or hole and nesting use the returned coordinates, not source data;
  overlapping shells that the audit cannot see are refused when detected and
  otherwise stay uncertified.
- The crossing, contact and containment tests are argued complete for polygon
  overlap, not proven. Contacts are decided locally and exactly: rings that
  cross through a shared vertex or a vertex on an edge are refused, as are a
  repeated vertex and a spike in one ring. What remains uncertified is an
  overlap of shells that no ring edge reveals (ADR 0033).
- The nesting is quadratic in the worst case; the reference is an oracle, not
  a performance tier.

**Follow-ups / risks to watch**

- Pair the boundary edges at a pinch vertex by angle instead of refusing.
- `Both` regions (ADR 0091, landing 3) must say what the cut means there.

## Relation to existing code

- `crates/contracts/operations/mesh-section/src/cut.rs`: the types and
  refusal names; `request.rs`: the request, result and evidence fields;
  `fixtures.rs` and `region_conformance.rs`: the suite.
- `crates/execution/dispatch/src/section_request.rs`: the answer checks.
- `crates/algorithms/reference/src/section/`: `request.rs` (which region
  edges bound the cut), `cut_regions.rs` (ring tracing), `ring_nesting.rs`
  (nesting).
- `scripts/probe_section_cut_regions_mutants.py`: the mutation probe.
