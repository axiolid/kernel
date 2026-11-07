# 0093 — Mesh sections merge in-plane faces on an exact arrangement

- **Status:** Accepted
- **Date:** 2026-10-29
- **Deciders:** Friedrich Schrödter
- **Supersedes:** —
- **Amended by:** [0094](./0094-visible-face-regions-share-the-hidden-line-view-and-an-exact-arrangement.md): the arrangement moves to `axiolid-exact-arrangement`. [0100](./0100-exact-arrangements-are-snap-rounded-to-the-binary64-grid.md): the sub-ulp collapse that this ADR refuses is snap-rounded for visible faces, and for sections once #435 lands.
- **Amends:** [0091](./0091-mesh-sections-report-in-plane-regions-and-open-chains.md), [0092](./0092-mesh-sections-report-the-solids-cut-as-closed-regions.md)

## Context

Issue #366 (consumer: rheonda/cad). ADR 0091 refused two input families by
name: in-plane faces of one side that overlap, and in-plane faces of opposite
sides that overlap. Both are ordinary in building elements. A column on a slab
puts the column's bottom (facing down) on the slab's top (facing up); stacked
wall parts have flush faces; two bodies of one element overlap in a floor
plan. ADR 0091 left `Both` unused and "the merge of overlapping same-side
faces" to a later landing (through `axiolid-overlay` or a focused provider),
and ADR 0092 left open what the cut means where `Both` appears.

The maintainer decided (comment on #366): an exact 2D arrangement of the
in-plane triangles inside the reference section module, on the exact plane
points of #377; a `Both` region is part of the cut; faces of one side that
only share a boundary segment merge; the two overlap refusal constants stay;
sub-ulp slivers are reported exactly or refused as `Degenerate`.

## Decision

**Semantics.** Let `Up` be the union of the in-plane faces whose material lies
below the plane (`InPlaneSide::Below`) and `Down` the union of the others.
Every point of the plane gets one class: `Below` = `Up` minus `Down`, `Above`
= `Down` minus `Up`, `Both` = `Up` and `Down`, otherwise empty. A region is a
maximal edge-connected part of one class. Regions of all sides have pairwise
disjoint interiors. Faces of one side that overlap, or only share a boundary
segment (including faces of different shells), form one region, so the
regions depend on the area covered and not on the triangulation. Pinched
shapes keep ADR 0091's rules: parts that meet only at a vertex are separate
regions, and a hole touching the outer boundary at a vertex gives one weakly
simple outer ring.

**`Both` is part of the cut** (amends 0092). Its points have the solid's
material on both sides, which is 0092's own definition of the cut. With cut
regions requested, every `Both` region equals one cut region, outer ring and
holes, bit for bit. This keeps 0092's identity: just below the plane is the
cut together with the `Below` regions, just above is the cut together with the
`Above` regions. A box on a slab: `Below` is the slab top minus the footprint,
`Both` is the footprint and is the cut. Two boxes stacked face to face: the
footprint is `Both` and the cut, so the wall does not vanish from a drawing.
Where a `Both` region meets transversal cut area they are two cut regions that
share an edge in opposite directions, as touching shells already are (0092,
#377).

**Algorithm** (`axiolid-reference`, `section/arrangement.rs` and
`regions.rs`), on `PlanePoints` of ADR 0092's #377 amendment and without a
tolerance:

- The unique in-plane source edges are cut at every proper crossing, at every
  end of another edge lying strictly inside, and so at T-junctions and
  collinear overlaps. A crossing is a homogeneous point from the cross
  product of the two lines (degree 3, 3, 2 in the coordinates, so `orient`
  among them has degree at most 8 and the cut rings' degree 11 is not
  exceeded); a binary64 position takes the fast path. Equal points are
  joined (three edges through one point, a crossing at a source vertex, source
  vertices at one position), so each point is one vertex.
- A piece is the stretch between two consecutive vertices of an edge. A
  triangle that carries it covers the side its third vertex is on; any other
  triangle covers both sides if the piece lies in its interior, which is
  decided by exact signs because a triangle is convex and no edge crosses a
  piece. The class of a side follows.
- A piece whose sides differ is a boundary half-edge of each non-empty class.
  The pieces at a vertex are ordered by direction with exact signs, and the
  half-edge after an incoming one is the first boundary of its class met
  turning clockwise. This resolves pinches fan by fan. A cycle is an outer
  ring or a hole by the exact turn at its extreme vertex; a hole belongs to the
  innermost outer ring of its class that holds it. A collinear vertex is
  dropped only if exactly two boundary pieces meet there and it is not a cut
  segment node, so adjacent regions keep their shared vertices.
- The coordinates are rounded once per vertex, so the regions of all classes
  and the cut share bit-identical points.
- The rings of `Both` regions join the cut rings in the nesting of ADR 0092. A
  cut request also refuses a `Both` ring that visits a vertex twice as
  `REFUSED_PINCHED_CUT_REGION`, and a bounding region edge with in-plane
  faces on the cut's side of it as `REFUSED_CROSSING_CUT_RINGS` (the shells
  overlap in volume).

**Refusals that remain.** Saddle vertices, non-manifold input, budgets and the
legacy `section` path. `Degenerate` where an exact arrangement vertex rounds
to the same output point as its neighbour on a ring ("in-plane region
collapsed under projection"): a sliver thinner than an ulp that a consumer's
binary64 world transform makes between faces that were exactly collinear is
reported as its own region when its points stay distinct, and refused when
they collapse. Any tolerance that would read such slivers as contact is a
separate decision (#377). Self-intersection that no in-plane piece reveals
stays uncertified (ADR 0033).

**Contract** (`axiolid-mesh-section-contract`), additive: `PlaneSectionEvidence`
gains `in_plane_crossings` (distinct points other than source vertices where
in-plane edges cross) and `with_in_plane_crossings`. `in_plane_pairs` and
`RegionLimits::max_in_plane_pairs` count the candidate pairs examined for the
in-plane faces: segment pairs, piece-triangle pairs and hole containment pairs.
`InPlaneSide::Both` documents the cut rule. `REFUSED_OVERLAPPING_IN_PLANE_FACES`
and `REFUSED_IN_PLANE_FACES_ON_BOTH_SIDES` stay, not deprecated, as the names
a provider that does not merge may use; the reference stops returning them.
The two overlap fixtures stay `OutputOrRefused`, and the new stacked fixtures
are too. The registry drops its refusal of cut regions beside a `Both` region
and checks instead that each `Both` region equals one cut region, found by
binary search on the sorted outer rings; it makes no geometric overlap check
on the returned rings (they are rounded: #274, #377).

**Budgets.** One budget, `max_in_plane_pairs`, covers all pair work. The
crossings are bounded by the segment pairs and charged against the memory
budget as they are created (ADR 0088); the rest of the working set is
declared per element, 2560 bytes (was 1536), tallied beside
`request_scratch_requirement`.

## Alternatives considered

| Option | Why not |
| --- | --- |
| (a) `axiolid-overlay` | Its tolerance snapping (#222) moves vertices, which contradicts "no tolerance takes part in a decision" and breaks the bit-identity of chain ends with region vertices. It works on binary64 points in the frame, which is the #377 failure. It returns regions without the provenance the cut needs, and it adds `axiolid-overlay` and `i_overlay` to every consumer of `axiolid-reference`. |
| (b) A focused section provider package | Packaging, not an algorithm: it still needs (a) or this inside. Moving `ScalarSection` out behind a re-export is not breaking and can be its own issue. |
| (c) Keep the refusals | Every element made of bodies with flush faces stays refused. |
| Deprecate the two refusal constants | It would force every provider to merge. |
| Leave `Both` out of the cut | It breaks 0092's identity: stacked bodies would vanish from the hatching. |

## Consequences

**Positive**

- A column on a slab, stacked wall parts and overlapping bodies of one element
  section into regions without a refusal; the cut is correct where bodies are
  flush.
- Every decision is an exact sign; the outcome does not depend on the frame,
  and a mesh rotated in binary64 gives the verdict of its rounded input.
- No new dependency and no closure change: `axiolid-reference` already
  depends on `axiolid-exact` since #377.

**Negative / costs**

- Every section with in-plane faces runs the arrangement, with no fast path.
  The reference stays an oracle, not a performance tier. Measured in wall
  clock (`cargo bench -p axiolid-reference --bench section_in_plane`, best of
  five or three, the same bench on the commit before the merge, 8e7f939; its
  `section_in_plane.rs` does not compile there because it reads
  `evidence.in_plane_crossings`, so add the `[[bench]]` entry to its manifest and drop that one field from the `format!` in
  `time` to reproduce the baseline): a conforming grid of
  8,192 in-plane triangles takes 31 to 35 ms, against 15 to 21 ms before (1.5
  to 2.4 times). A first version took 61 ms; pruning the segment and piece
  pairs (a side test that stops most pairs after two signs, skipping segments
  of one source vertex, holding a piece to the triangles whose box can contain
  it) halved it. The cases the merge newly answers were refused before and
  have no baseline: the grid with a square footprint of the other facing, 8,194
  triangles, 32 ms; two skew grids of 2,304 triangles with 3,242 crossings, 29
  ms; a box on a slab, 24 triangles, 15 microseconds.
- The arrangement is exact and its crossing points are big integers, so faces
  that cross cost more than conforming ones: the skew grids take about as long
  as a grid of four times the triangles.
- Elements whose bodies overlap in volume are still not certified (ADR 0033);
  the cut refuses the overlaps that an in-plane piece reveals.

**Follow-ups / risks to watch**

- rheonda/cad adopts `InPlaneSide::Both`: hatch `cut_regions` (they now
  include `Both`), draw `Below` as the visible outline and `Above` as hidden,
  and treat `Both` as cut area when cut regions are not requested.
- A request-level `max_in_plane_crossings` limit, if the run-time memory charge
  is not wanted.
- Moving `ScalarSection` into a focused package with a re-export (ADR 0036).

## Relation to existing code

- `crates/algorithms/reference/src/section/`: `arrangement.rs` (the pieces and
  their sides), `regions.rs` (rings from the arrangement), `plane_points.rs`
  (`segment_crossing`, `segment_parameter`, point location), `request.rs`
  (the cut wiring), `ring_nesting.rs` (nesting of the `Both` rings).
- `crates/contracts/operations/mesh-section/`: `request.rs` (classes,
  evidence), `cut.rs` (the invariant), `fixtures.rs` and
  `region_conformance.rs` (the suite).
- `crates/execution/dispatch/src/section_request.rs`: the answer checks.
- `crates/algorithms/reference/benches/section_in_plane.rs`: the bench.
- `scripts/probe_section_merge_mutants.py`, with
  `probe_section_regions_mutants.py` and `probe_section_cut_regions_mutants.py`.
