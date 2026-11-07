# 0100 — Exact arrangements are snap-rounded to the binary64 grid

- **Status:** Accepted
- **Date:** 2026-11-06
- **Deciders:** Friedrich Schrödter
- **Supersedes:** —
- **Amends:** [0093](./0093-mesh-sections-merge-in-plane-faces-on-an-exact-arrangement.md), [0094](./0094-visible-face-regions-share-the-hidden-line-view-and-an-exact-arrangement.md)

## Context

Issue #410 (found while landing #383; #398 was closed as its duplicate). The
visible face regions of ADR 0094 are decided on an exact arrangement and
their ring vertices are rounded once to binary64. In a generic rotated scene
two distinct exact vertices can round to the same binary64 point, and the
region then collapses: the request was refused as `Degenerate`, typed and
naming the mesh. That is fail-closed, but it refuses valid scenes. The
independent ray-casting oracle of #383 hit it in 2 of 168 seeded scenes
(seed 1 oblique and seed 9 isometric, each with a near plane).

ADR 0093 refused a sub-ulp collapse of a section region in the same way and
called any tolerance for slivers "a separate decision (#377)", and ADR 0094
pointed its known limit at "the snapping decision of #377". No such decision
was made: #377 and the amendment of ADR 0092 decided cut rings on exact
points and left the returned coordinates unchanged. This ADR makes the
decision, for the faces now and for the sections when they follow (#435).

A first attempt dropped each rounded vertex that equals its predecessor and
validated the rounded rings (at least three vertices, simple, no crossing
rings, the exact orientation). It did not fix either seeded scene: in both
the whole ring of a region lies inside one rounding cell, so nothing remains
to validate, and the only choices are to refuse or to drop. Dropping without
a record is a silent degradation, so the decision below carries the record.

## Decision

**Single-pass snap rounding on the binary64 tensor grid, decided exactly,
followed by an exact re-arrangement and re-trace of the snapped boundary,
then certification.**

- **Hot cells.** A hot cell is the set of points that
  `ExactPoints::round` rounds to the rounded position of an end of an input
  edge. The rounding is to nearest with ties to even, so a cell is a box
  `[(prev + c) / 2, (c + next) / 2]` per axis whose bounds belong to the
  cell exactly when the centre `c` has an even significand: the cells are
  half-open and partition the plane, whatever the binade, the subnormals or
  zero.
- **Edges and routing.** An *edge* is one exact ring edge: the straight run
  between two consecutive kept vertices, once for all sources. Each edge is
  routed through the centres of every hot cell that some point of the closed
  exact segment rounds to, in order along the segment. Where nothing but an
  edge's own end cells is met, the fragment is the edge with its ends
  rounded once: the output is today's.
- **Re-arrangement.** The routed fragments are arranged exactly with
  `arrange`. A rounded vertex that lies on another fragment splits it and
  becomes a shared vertex.
- **Per-source re-trace.** Each source's chain of exact ring edges is mapped
  onto the snapped pieces and netted: opposite traversals of a piece cancel,
  so a part that became zero-width vanishes and a zero-width gap merges its
  two sides. The net half-edges are traced with the existing `trace`, whose
  rule for kept vertices now decides on the snapped boundary.
- **Certification.** Any crossing left among the fragments, a piece with net
  count two or more, or a point of winding other than 0 or 1 for a class
  (`check_winding`) is refused as `Degenerate`. None has occurred on the
  seeded scenes or in the property test.
- **What a region is.** The regions are the exact regions snap-rounded. The
  rings are weakly simple: they may touch at a vertex, never cross, never run
  along each other, and no ring vertex lies inside the edge of any ring. A
  region may be split at a pinch. `VisibleRegion::triangles` lists the
  triangles frontmost somewhere in an exact region this region was snapped
  from, so a pinch-split region lists the whole set in each part.
- **Nothing collapses silently.** A region or hole whose area vanishes is
  not returned, and the evidence counts it: `collapsed_regions` and
  `collapsed_holes`. `rerouted_edges` counts the edges routed through a cell
  other than those of their ends, and `snap_displacement`, a new
  `DisplacementBound` (finite, non-negative, `Eq`), bounds how far any
  boundary moved in the maximum norm. All are zero when every ring vertex is
  a binary64 point and no edge was rerouted. A region is never refused
  because rounding collapses it.

**Displacement and winding.** The bound is the largest half-extent of a cell
that a route visits, per axis the larger gap to a neighbouring binary64
number halved, and zero for an end that is already a binary64 point. Claim:
an image point farther than the bound, in the maximum norm, from the exact
boundary of a source lies in that source's regions exactly when it lies in
its exact regions. Proof: take consecutive centres `c_i, c_{i+1}` of a
route. Each lies within the bound of a point `p_i`, `p_{i+1}` of the exact
edge `e` in its own cell, and the points follow each other along `e`. The
point `(1-s) c_i + s c_{i+1}` is within the bound of `(1-s) p_i + s p_{i+1}`,
a point of `e`, because the maximum norm is convex. So the snapped chain lies
in the tube of the bound around the exact one, and so does every point of
the straight-line homotopy between them. The point is not on the homotopy,
so the winding number of the class about it does not change; netting only
removes opposite traversals, which add nothing to a winding. The proof does
not depend on the grid being uniform.

## Alternatives considered

| Option | Why not |
| --- | --- |
| Drop rounded vertices that equal their predecessor, then validate the rounded rings (the first attempt) | Does not fix the seeded scenes: the whole ring of the collapsed region is in one cell. It can only refuse or return a non-region. |
| Iterated snap rounding (Halperin and Packer) | The drift is unbounded (O(n) cells). The case it targets, a vertex rounding onto an edge, is already exact here because the fragments are re-arranged. |
| Stable snap rounding (Hershberger) | Bounded, but considerably more complex. A candidate for the general API of #135. |
| A uniform coarse grid (the ulp of the largest coordinate) | Moves every vertex, including binary64 inputs, and breaks bit identity with the hidden-line endpoints. |
| Drop collapsed regions without a count | A silent degradation; the counts make the dropped area auditable. |
| Exact rational output | An API break; the consumers need binary64. |

## Consequences

**Positive**

- Seeds 1 and 9 are answered: each loses one region of no width
  (`collapsed_regions` 1; seed 9 also reroutes two edges), with a
  displacement of 2^-52 and 2^-53. The oracle's 168 scenes all pass the
  ray-casting comparison and an independent topology checker (exact
  `Dyadic` signs: no crossing, no vertex inside an edge, exact orientation).
- Where nothing is rerouted the output is today's: vertices are the exact
  points rounded once and equal the hidden-line endpoints of the same points.
  A sub-ulp bend that is straight in binary64 is now dropped from a ring.
- The new items are additive: `DisplacementBound`, four fields on the
  `#[non_exhaustive]` `VisibleFaceEvidence`, `snap_round` and its types, and
  `check_winding`. No `semver-exceptions` entry. The contract wording
  tightens an obligation (a region is not refused for collapsing); the
  versions are the maintainers' call.

**Negative / costs**

- Every request re-arranges and re-traces its boundary. Wall clock of
  `scenario/visible-faces`, the base `675f09e2` against this change, best
  of three alternating runs, in milliseconds (mean per iteration):

  | case | base ms | change ms | ratio |
  | cut-plan/20x20/bvh | 102 | 86 | 0.85 |
  | cut-plan/20x20/scan | 150 | 111 | 0.74 |
  | cut-plan/20x20/unpruned | 112 | 83 | 0.74 |
  | cut-plan/50x50/bvh | 993 | 1141 | 1.15 |
  | cut-plan/50x50/scan | 4172 | 4420 | 1.06 |
  | cut-plan/50x50/unpruned | 1060 | 872 | 0.82 |
  | isometric/10x10/bvh | 191 | 201 | 1.05 |
  | isometric/10x10/scan | 190 | 204 | 1.07 |
  | isometric/10x10/unpruned | 701 | 670 | 0.96 |
  | isometric/20x20/bvh | 1363 | 864 | 0.63 |
  | isometric/20x20/scan | 1334 | 885 | 0.66 |
  | isometric/20x20/unpruned | 3553 | 3190 | 0.90 |
  | overlap-isometric/6x6/bvh | 330 | 222 | 0.67 |
  | overlap-isometric/6x6/scan | 325 | 224 | 0.69 |
  | overlap-isometric/6x6/unpruned | 560 | 480 | 0.86 |
  | overlap-plan/6x6/bvh | 33 | 27 | 0.80 |
  | overlap-plan/6x6/scan | 32 | 26 | 0.80 |
  | overlap-plan/6x6/unpruned | 39 | 32 | 0.81 |

  The machine was shared with other jobs and the run-to-run spread is large
  (up to 50%), so read the ratios as noise-limited: the snap costs nothing
  measurable on most cases (the ratios below 1 are noise, not a speed-up),
  and the largest rise is 15% on the 50x50 cut plan with the BVH. No gain
  and no loss beyond that is claimed. The bench's assertion that the three
  variants return identical regions held on the change.

- The no-crossing theorem of snap rounding is proved for uniform grids. It
  carries over to any box where both axes stay within one binade (an axis
  scaling of a square grid). It is unproven for fragments that cross a power
  of two or zero. Every such failure is refused, never returned; the
  property test of `snap_round` (2,000 seeded cases of non-binary64 triangles
  near zero, 0.5, 1 and 2 at sizes of 2^-40 to 2^-70) aims at those cases and
  found none. A counterexample would be a refusal, not a wrong answer: the
  displacement and winding guarantee does not depend on the grid.
- Fills and hidden lines still meet exactly except along rerouted edges and
  collapsed parts, which the evidence reports. Sharing the hot points of the
  two requests would need a combined request, a later decision.
- A no-event fast path (skip the re-arrangement when nothing is rerouted) is
  sound only under the uniform-grid argument, so it would be an issue of its
  own.

**Follow-ups / risks to watch**

- **Interim gap.** Sections still refuse "in-plane region / cut region /
  section chain collapsed under projection" where visible faces now snap, so
  the two operations disagree until #435. The sections' arrangement lives in
  projected coordinates while their output is in plane-frame coordinates, so
  their binary64 grid is not axis-aligned in the arrangement, and their
  chains must share the snapped points with the regions. They need exact
  plane-frame points and one `snap_round` over region and chain edges.
- Patch versus minor for `axiolid-contracts`, the contract crate,
  `axiolid-exact-arrangement` and `axiolid-hidden-lines` is the maintainers'
  call.

## Relation to existing code

- `crates/contracts/common/base`: `DisplacementBound`.
- `crates/algorithms/planar/exact-arrangement`: `src/snap.rs` (`snap_round`,
  `Snapped`, `ClassBoundary`, the private cell and meeting decisions),
  `check_winding` in `src/trace.rs`, `ExactPoints::frame_sign`; tested by
  `tests/snap.rs` and probed by `scripts/probe_snap_rounding_mutants.py`.
- `crates/contracts/operations/mesh-hidden-lines`: the four evidence fields,
  the contract wording and the frozen snap cases of the conformance suite.
- `crates/algorithms/discrete/hidden-lines`: `src/faces/regions.rs` builds
  the exact regions, then the snapped ones; `tests/faces_snap.rs`,
  `tests/common/topology.rs` and `scripts/probe_visible_faces_mutants.py`.
