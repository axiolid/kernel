# 0108 — Minkowski sums are read from the term edges that can bound them

- **Status:** Accepted
- **Date:** 2026-11-07
- **Deciders:** Friedrich Schrödter
- **Supersedes:** —
- **Amends:** [0087](./0087-disc-morphology-has-an-exact-and-a-fast-mode.md)

## Context

Issue #344, with #292 as its first consumer: a downstream rule engine
searches furnished rooms for a place a rectangle fits, by summing furniture
outlines with the rectangle and eroding the free floor by it. A sum or
erosion by a convex shape is the union of many convex terms: the shape's hull
at both ends of each boundary edge, a copy of the set moved by the shape's
anchor, and so on. `sum_parts`, `erode_parts` and the fast disc morphology
put every term ring into one `ArcArrangement` and read each face with
`Rule::holds`. With a shape larger than the region's features every hull
overlaps nearly every other, so the arrangement has `O(n^2)` crossings and
every piece lies in `O(n)` hulls. The committed benchmark
`benches/minkowski_furniture.rs` measures it.

The results are pinned bit for bit (digests in `tests/minkowski.rs` and the
benches, ADR 0087's "`Exact` bit for bit"), so the fix must not change which
terms are summed, nor how a vertex is rounded.

## Decision

Keep every term and the rule that reads the union; change only how the union
is evaluated (`src/term_union.rs`).

1. **Certify and drop interior term edges.** An edge of a term ring is left
   out when some convex *piece* ring (a hull, chord hull, trimmed hull or
   parallelogram), other than the edge's own, holds its open segment strictly
   inside: both ends in the closed piece, and no side line of the piece
   holding both ends. The test is the exact orientation of `f64` points.
   Pieces *absorb*: the rule's answer is the same for every list of rings
   holding one, since a sum holds such a point and an erosion leaves it out.
   So a dropped edge never bounds the result.
2. **Cut only the kept edges, at every term vertex.** The exact core gains a
   partial-ring mode of its subdivision (`build_partial`), sharing its body
   with the ring build: each ring contributes a chosen subset of its edges,
   every given point inside a kept edge becomes a stop, and no winding
   numbers are counted.
3. **Read each piece with the same rule.** The rings holding a side of a
   piece are the carriers whose inside it is that side; a piece ring other
   than a carrier holding the piece's sample strictly is added to both sides
   (absorbing); otherwise the set rings with a nonzero winding number at the
   sample are. `Rule::holds` then decides each side, unchanged.
4. **Link, round, present and settle as `region_of` does** (`region_from`).

The result equals what the arrangement route gives for the same terms, bit
for bit: vertex values, the start of every ring, the order of rings and
holes, for `Exact` and `Fast` alike. The argument is in the module docs of
`term_union.rs`. The arrangement route is kept only as a test oracle
(`by_arrangement`); tests compare `polygons()` with `assert_eq!`. Refusals
are unchanged and no error variant is added.

## Alternatives considered

| Option | Why not |
| --- | --- |
| Textbook reduced convolution (Behar-Lien; CGAL `Minkowski_sum_by_reduced_convolution_2`) on vertex sums, settled by winding | It defines a different set at the ulp level: the terms are hulls of vertex sums rounded once, and #323's trimmed hulls, so segments between rounded sums would not reproduce today's boundary, and the digests and ADR 0087's "`Exact` bit for bit" would break. Its winding of the full convolution equals the Euler characteristic of `R ∩ (x − K)`, which is 0 when `x − K` covers a hole of `R` while meeting `R`; the erosion's complement (a frame minus a small room) has exactly that case, so it still needs a per-face collision test. Exact vertex sums need non-`f64` ring vertices in the exact core (ADR 0070). |
| Divide-and-conquer union of the hulls | Intermediate unions have crossing-point vertices. Rounding them changes the result, and keeping them exact needs non-`f64` input rings. |
| Convex decomposition of `R`, one slope-merged convex sum per piece | Different terms, so a different rounding. Concave tessellated outlines cut into many pieces that all overlap under a large shape. |
| `axiolid-exact-arrangement` (ADR 0094) for the segments | A second exact core inside overlay against ADR 0070, and a new internal dependency edge (a closure change). Its rounding path would also have to be proved identical. |

## Consequences

**Positive**

- The measured costs are in `docs/research/minkowski-sums-from-bounding-term-edges.md`.
- Results, refusals and both morphology modes are unchanged bit for bit, so
  no digest moves.
- The exact core's partial mode and `Windings` are reusable by any caller
  that knows which rings can hold a piece.

**Negative / costs**

- Two subdivision modes in the core (one shared body), and per-piece
  queries (a box-tree search per kept edge and per piece).
- Teeth and needle walls under a disc keep their quadratic `Exact` cost,
  because their arcs and slivers at every tip are not inside one convex term.
  `Fast` remains the answer for those walls (ADR 0087).
- Where a side normal of the shape falls in a reflex turn, that side of the
  shape at that corner stays, since no single term holds it.

**Follow-ups / risks to watch**

- Certification by the union of two terms, if a measured case needs it.

## Relation to existing code

- `crates/algorithms/planar/overlay/src/term_union.rs` (`region`, the proof
  and the cost argument).
- `crates/algorithms/planar/overlay/src/exact_arc/arrangement.rs`
  (`subdivide`, `build_partial`, `Windings`, `Raw::sample`),
  `exact_arc/boxes.rs` (`BoxTree::first`), `exact_arc.rs` (`inside_convex`).
- `crates/algorithms/planar/overlay/src/minkowski.rs` (`sum_parts`,
  `erode_parts`, the fast disc branch, `region_from`, `Rule::first_piece`,
  the test oracle `by_arrangement`) and `src/arrangement.rs`
  (`oriented_rings`).
- `scripts/probe_term_union_mutants.py`;
  `tools/benchmark/benches/minkowski_furniture.rs`.
