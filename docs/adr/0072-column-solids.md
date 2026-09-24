# 0072 — Column solids over one planar arrangement

- **Status:** Accepted
- **Date:** 2026-09-24
- **Deciders:** Friedrich Schrödter
- **Supersedes:** —

Part of axiolid/kernel#120.

## Context

Two coaxial booleans and one plane cut were still refused by name after
ADR 0071:

- a union of prisms with different spans, and a difference whose tool stops
  inside the subject (a counterbore, a blind pocket, a slot through the
  middle heights) — the result is *stepped*;
- a sloped plane that crosses a prism's top or bottom cap inside the
  section — part of the old cap survives next to the cut.

None of these is one prism, but all of them are the same shape: the plan is
cut into cells, and above each cell the solid occupies a stack of height
intervals, each bounded by a plane (flat, or sloped as in ADR 0071). Which
intervals a cell carries depends only on which input rings contain it.

`boolean_stepped::union_prisms_stepped` already returned such a result as
bands, but left the ledge faces between bands unbuilt. Building each band
with its own overlay and gluing them fails for a structural reason: two
overlays round crossing points independently, so the ledge and the band
walls would disagree on vertices by an ulp and the shell would not close.

## Decision

Build every such result from **one exact arrangement of all input rings**
and **one column builder** over it.

- `axiolid-overlay` gains `ArcArrangement`: every ring's edges split at
  every crossing, each piece labelled with the rings containing its left
  and its right side, shared boundary pieces merged and naming every ring
  that carries them (`EdgeSource`: ring, edge, direction). `regions(pred)`
  links the pieces bounding any membership predicate into outer and hole
  rings. It reuses the exact predicates and ring assembly of `arc_overlay`
  (ADR 0070), which now shares that assembly.
- `axiolid-construct::column` (crate-private) takes the arrangement, a list
  of planes and a function from membership mask to height blocks, and
  emits: a cap per (plane, facing) over the cells where a block ends there;
  a wall per arrangement piece over each height run where exactly one side
  is solid, facing the empty side; vertical edges split at every height any
  face meets at a vertex, so walls meeting there share edges. Faces are
  grouped into shells through shared edges; each connected shell with
  positive enclosed volume is one `ExactBRep`.
- Curved walls stay `Cylinder` faces; sloped rims keep the ellipse edge and
  `Sinusoid2` pcurve of ADR 0071, built by the same helpers.

The public entry points change behaviour, not signature:

- `boolean_prisms_exact`, `boolean_arc_prisms_exact` and their `_solids`
  variants build stepped results instead of refusing them. The equal-span
  cases keep their previous path.
- `clip_arc_prism_exact` builds a plane crossing a cap.
- Names: a wall is named after the operand edge it lies on
  (`Side(n)` counted across that operand's rings, as its own extrusion
  counts them), a cap after the operand whose own cap lies on that plane
  (subject first), and the floor or ceiling a difference leaves after the
  tool cap that made it. Only the tool can leave one: the subject's own
  opposite cap never bounds a coaxial result.

**Still refused, by name:**

- A result enclosing a cavity (a tool buried inside the subject). The
  builder can build the void shell, but `tessellate` reads only
  `solids()[0].outer` and treats voids as boolean intent, so a void would
  be lost silently downstream. The refusal names the cavity.
- Pieces of solid that touch only along an edge: not a manifold solid.
- Non-coaxial curved booleans and cones, spheres, tori and NURBS. Those
  need general surface/surface B-rep booleans (OCCT's `BOPAlgo` scale); the
  column reduction does not apply and nothing here pretends it does.

## Alternatives considered

| Option | Why not |
| --- | --- |
| One overlay per band, glue the bands | Crossing points round independently per overlay; shared vertices disagree by an ulp and the shell cannot close without a tolerance-based sewing step. |
| Snap per-band vertices within tolerance | Exactly the tolerance guessing ADR 0070 removed from the planar path. |
| Keep returning bands and let callers assemble | Every caller would rebuild the ledge logic; the band list is also not a solid, so it cannot enter the exact B-rep pipeline. |
| General B-rep boolean now | Correct end state for non-coaxial cases, but orders of magnitude more work; the vertical-column family covers the building-model cases #120 names (openings, counterbores, stepped footings, sloped roofs over columns). |

## Consequences

**Positive**

- Every vertical-column result shares one vertex table, so faces agree on
  every vertex by index and the audit checks them against each other.
- One builder replaces three special cases (single prism, sloped cap,
  bands) for everything it covers; the single-prism paths remain for their
  existing inputs.
- Found and fixed a measurement bug on the way: `exact_properties` ignored
  face orientation, which only showed on solids off `z = 0`.

**Negative / costs**

- The builder is crate-private and specialised to vertical columns; a
  general boolean will need its own topology construction.
- Cavities stay refused until tessellation and measurement handle void
  shells.

**Follow-ups / risks to watch**

- Non-coaxial curved booleans and cones/spheres/tori need their own issue.
- `boolean_stepped::union_prisms_stepped` is now the lighter alternative,
  kept for callers who want bands; its volumes are checked against the
  solid's.

## Relation to existing code

- `crates/algorithms/planar/overlay/src/arrangement.rs`,
  `src/exact_arc/arrangement.rs`: `ArcArrangement`.
- `crates/algorithms/construction/construct/src/column.rs`: the builder.
- `crates/algorithms/construction/construct/src/boolean_column.rs`: the
  coaxial and clip adapters.
- `crates/algorithms/construction/construct/src/boolean_exact.rs`: entry
  points, now routing stepped and crossing cases to the adapters.
- `crates/algorithms/query/measure/src/exact.rs`: the orientation fix.

## Verification

- 11 column tests (`column_booleans.rs`), every volume derived by hand
  from the inputs: a slab-and-tower union (48), a corner notch (114), a
  blind pocket (104), a round union with the lens subtracted, a counterbore
  (11 pi), a strip that splits a disc only below z = 1 (one solid, not
  two), the single- and multi-solid paths agreeing, a hole's wall numbering
  after the outer ring, and both non-manifold refusals. Every solid passes
  the geometric audit.
- The clip tests that pinned the crossing-cap refusal now check the built
  solid against a closed form for the kept volume (`disc_excess`, the
  integral of `max(0, x - a)` over the unit disc).
- `boolean_stepped` bands and the stepped solid agree on volume.
- 7 arrangement tests: labels against `arc_overlay`, inclusion-exclusion
  over three rings, shared pieces carrying both sources, source indices
  surviving ring reversal, crossings shared exactly.
- A 16-fault mutation probe, `scripts/probe_column_mutants.py`, catches all
  16. Its first run caught 14 of 17: nothing split a wall where it changes
  sides, and nothing pinned hole-wall numbering, so tests for both were
  added. The seventeenth fault targeted a naming branch no boolean can
  reach; the branch was removed. Flipping cap facing by hand fails the
  column tests, as a check on the probe itself.
- The measurement fix has its own test (a raised unit cube measured 7/3),
  which fails without the fix.
