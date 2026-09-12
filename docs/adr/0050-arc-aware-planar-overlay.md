# Arc-aware planar overlay

Status: proposed
Date: 2026-09-12

## Context

The exact solid boolean (`boolean_prisms_exact`) reduces a 3D boolean to a 2D
cross-section boolean crossed with a height-interval boolean. Its exactness
comes entirely from the planar stage. That stage is polygon-only:

```rust
pub struct Ring { pub points: Vec<Point2> }
```

So a cylindrical operand cannot enter the exact path at all, and
`Prism.rings: Vec<Vec<Point2>>` cannot represent a disc. This is the root
cause of the reviewer's finding that "no curved-surface boolean exists".

## Measured starting point

The instruction that opened this work described the task as a rewrite of a
"971-line planar boolean engine". Measurement contradicts that premise:

- `axiolid-overlay` is **971 lines of wrapper**, not an engine. `overlay()`
  validates inputs, maps our enums onto the backend's, and calls out.
- The actual boolean algorithm is the third-party crate **`i_overlay`
  (4.5.2, MIT OR Apache-2.0), 24,098 lines**.
- `i_overlay` has no arc support. Grep hits for "arc" in its sources are the
  substring inside "se**arc**h".

Writing an arc-aware boolean from scratch therefore means replacing ~24k lines
of mature, integer-predicate-backed sweep-line code, not ~1k.

## Alternatives measured, not assumed

`cavalier_contours` 0.9.0 (MIT OR Apache-2.0, `#![forbid(unsafe_code)]`,
16,017 lines, MSRV 1.88 — matching our pinned toolchain) implements planar
booleans over polylines with **native circular-arc segments** (bulge
encoding), including arc/arc and arc/segment intersection and curved
containment classification.

Probe results, each checked against an independent closed-form area:

| case | measured | expected | error |
|------|----------|----------|-------|
| arc/arc lens (intersection) | 1.22836969860876 | 1.22836969860876 | 1.1e-15 |
| arc/arc union | 5.05481560857083 | 5.05481560857083 | 1.8e-15 |
| disjoint union (2 results) | 6.28318530717959 | 6.28318530717959 | 0.0 |
| wall minus round opening | pos 30 − neg pi | 30 − pi | exact |

Arcs are preserved as arcs through the boolean (`arcs=4`, `lines=0` on the
disc cases) — the output is not tessellated. Disjoint operands correctly
produce two result loops rather than a silently merged one. A fully contained
hole is returned in `neg_plines`, not `pos_plines`; an integration that reads
only the positive loops would silently drop holes.

## Decision

Adopt `cavalier_contours` as the arc-aware planar backend behind our existing
`axiolid-overlay` contract, rather than writing a new engine.

Rationale:

- The house rule is to check proven libraries, standards, and prior art
  before building custom. An arc-aware boolean with machine-precision
  agreement against closed-form areas already exists, is permissively
  licensed, is unsafe-free, and matches our MSRV.
- Our current design already treats the planar boolean as a swappable backend
  behind a validated contract. This is a second backend, not a new coupling.
- A from-scratch arc sweep-line is a multi-week numerically delicate project
  whose failure modes (arc/arc tangency, near-coincident arcs, containment of
  an arc inside an arc) are exactly the ones that take longest to get right.

Deliberately NOT decided here: replacing `i_overlay`. The polygon path is
mature and gated by existing tests. The arc backend is additive.

## Consequences

- `Ring` gains an arc-capable sibling; the polygon-only type stays, so every
  existing caller keeps compiling and the current tests keep their meaning.
- `Prism` gains an arc-capable cross-section, which is what finally lets
  `boolean_prisms_exact` accept a cylinder.
- Exactness claims must be re-stated honestly for the arc path: agreement
  with closed-form area was measured at 1e-15, which is machine-precision
  agreement, NOT the bit-exact integer predicate story `i_overlay` provides
  for polygons. Any doc that claims "exact" for the arc path must say which
  of the two it means.
- A second backend is a second set of numerical conventions. The contract
  layer must not leak `cavalier_contours` types.

## Plan

1. Contract types: arc-capable ring/segment in `axiolid-overlay`, polygon
   types untouched. **Done.** `ArcRing`/`ArcVertex` carry a bulge per
   departing edge (`tan(theta/4)`, DXF convention), with `arc_ring_area`,
   `arc_edge_radius`, `reverse_arc_ring` and `validate_arc_ring`. The
   existing `Ring`/`Polygon` types are unchanged, so nothing on the
   polygon path shifted.
2. Backend adapter behind the existing validation, with the neutral contract
   preserved (no library types in our public API).
3. Verification harness: every case checked against an independent
   closed-form area, plus the hole case that lives in `neg_plines`.
4. Mutation testing on the adapter: a dropped `neg_plines` loop, a dropped
   disjoint result, and a tessellated-instead-of-arc output must all be
   caught by the tests before the path is trusted.
5. Only then: widen `Prism` and wire `boolean_prisms_exact`.

## Open questions

- Does the arc backend's `neg_plines`/`pos_plines` split map cleanly onto our
  outer/holes `Polygon` shape in every operator, or only in the cases probed?
  Must be measured per operator before step 5.

## Degeneracy probe results

The degenerate cases the plan flagged as unverified were probed against
closed-form areas. All behave correctly, with no refusal needed:

| case | result | reading |
|------|--------|---------|
| externally tangent, AND | 0 loops | point contact is not an area |
| externally tangent, OR | 1 loop, 2pi, 4 arcs | both discs kept, joined at the touch |
| identical discs, AND | 1 loop, pi, 2 arcs | coincident boundaries handled |
| identical discs, NOT | 0 loops | empty, not a sliver |
| internally tangent, AND | 1 loop, pi/4 | equals the inner disc exactly |
| internally tangent, NOT | 1 loop, 3pi/4, 4 arcs | crescent, arcs preserved |
| contained disc, NOT | pos 1 loop pi, neg 1 loop pi/16 | hole in `neg_plines` |

Every area matches its closed form to the printed 14 decimal places. Notably
no case produced a degenerate sliver loop or a tessellated boundary, and the
zero-area cases return zero loops rather than an empty-but-present loop.

Remaining unprobed: zero-radius arcs and a full circle expressed as a single
segment. Those are input-validation concerns for the contract layer rather
than backend behaviour, and step 1 should refuse them explicitly.


## Step 1 findings

Two things the plan did not anticipate.

**The segment-area term must keep a signed angle.** A first implementation
using `|theta|` gave the right area for a counter-clockwise disc and `+pi r^2`
for a clockwise one, so a reversed ring failed to negate. Caught by
comparing against a closed form rather than against another call of the
same function. Mutation testing reproduces it: restoring `theta.abs()`
is killed.

**`ZeroRadiusArc` has a provably narrow trigger window.** Since
`R = chord (1 + b^2) / 4|b|` and `(1 + b^2) / 4|b| >= 1/2` with equality at
`b = 1`, every arc satisfies `R >= chord / 2`. A chord above tolerance `t`
therefore forces `R > t/2`, so a sub-tolerance radius is reachable only for
a chord in `(t, 2t]`. The first version of the test used a huge bulge on a
short chord, assuming that shrinks the radius; it does the opposite -- a
large bulge approaches a full circle and grows `R`. The test now targets the
proven window explicitly.

**Self-intersection is deliberately not validated for arc rings.** Arc/arc
and arc/segment crossing tests are part of the overlay algorithm itself
(step 2), and an approximate pre-check would either reject valid input or
admit invalid input. The straight-edge contract keeps its exact check; the
arc contract names the gap instead of implying a guarantee it cannot make.

Mutation evidence for step 1: 7/7 killed -- unsigned segment term, dropped
segment term (the tessellated reading), reversal without bulge negation,
reversal without the index shift, removed zero-radius guard, two-vertex
polygonal ring admitted, and unchecked bulge finiteness.
