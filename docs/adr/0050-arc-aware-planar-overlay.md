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
   preserved (no library types in our public API). **Done.** `arc_overlay`
   in `arc_overlay.rs`, returning `ArcPolygon` with holes. No
   `cavalier_contours` type appears in the public API.
3. Verification harness: every case checked against an independent
   closed-form area, plus the hole case that lives in `neg_plines`.
   **Done.** 7 tests in `tests/arc_overlay.rs`.
4. Mutation testing on the adapter: a dropped `neg_plines` loop, a dropped
   disjoint result, and a tessellated-instead-of-arc output must all be
   caught by the tests before the path is trusted. **Done.** All three
   named mutants are killed, plus a mis-wound hole.
5. Only then: widen `Prism` and wire `boolean_prisms_exact`. **Done.** Added
   `ArcPrism` + `boolean_arc_prisms_exact` alongside the polygon path (not
   replacing it), and `extrude_arc.rs` which sweeps an arc edge into a
   `Cylinder` face and a straight edge into a `Plane` face. Height logic is
   now shared via `resolve_span` so the two paths cannot disagree about
   which spans are representable.

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

## Step 2-4 findings

The `neg_plines` trap is real and the wall case proves it. A 10x3 wall
minus a unit-radius opening returns the wall in `pos_plines` and the
opening in `neg_plines`. An adapter reading only positives returns a
SOLID wall with a plausible area and no error -- the opening silently
disappears. That mutant is killed by asserting `holes == 1`, not by
any area check.

Edge counts, not areas, are what catch tessellation. A backend that
approximated arcs by segments would still report an area within any
reasonable tolerance. `ArcOverlayEvidence::arc_edges` and `line_edges`
make that observable: an annulus must have zero straight edges and a
handful of arcs, and a lens must have no straight edge at all.

Mutation results on the adapter:

| mutation | result |
|----------|--------|
| drop `neg_plines` (holes vanish) | KILLED |
| keep only the first positive region | KILLED |
| bulge dropped in conversion (tessellating) | KILLED |
| hole not reoriented (winding lie) | KILLED |
| input normalisation removed alone | survives (redundant) |
| output reorientation removed alone | survives (redundant) |
| ALL orientation handling removed | KILLED |

The two lone survivors are genuine redundancy rather than a test gap:
normalising the output repairs a mis-wound input, so either mechanism
alone still produces the right answer. Both are kept and the source
says why. This was verified by removing them together, which fails.

Areas agree with closed forms to 1e-12 in every case: wall minus
opening `30 - pi`, lens `2 acos(1/2) - sqrt(3)/2`, annulus `3 pi`,
disjoint union `2 pi`.

## Step 5 findings

A curved-surface exact boolean now exists. A unit cylinder intersected
with a box returns an `ExactBRep` carrying `Cylinder` wall faces and
`Plane` caps -- not a fan of planar strips.

Three things worth recording.

**A surface-kind count is not enough.** A wall can be a `Cylinder` of the
correct radius and still sit in the wrong place: the arc centre comes from
the bulge by way of a sagitta offset, and dropping that term leaves the
radius intact while moving the axis to the chord midpoint. That mutant
survived a radius check and a face-kind check; only asserting the wall axis
passes through the disc centre killed it.

**The cap pcurve has to be the arc, not its chord.** A cap loop built from
chords would disagree with the wall loop along the same edge while still
closing, still validating, and still reporting a plausible area. The cap
pcurve is a `Circle2` for arc edges for that reason.

**A hole is refused, not filled.** A plate minus a centred disc has an
interior opening, and the arc extruder does not build a cap face with two
bounds. Returning a solid plate would be a silent geometry change, so the
case refuses with `interior hole` in the reason.

Mutation testing, 4/4 killed: every wall built planar, wrong
bulge-to-angle constant, centre ignoring the sagitta, hole refusal removed.

Remaining gap: results with interior holes, and sections whose arcs survive
into more than one disconnected region. Both refuse by name.

## Follow-up: stepped union and general-polygon fillet

Two gaps named in the thesis feedback are now closed.

### Stepped union

`boolean_prisms_exact` refuses a union of differing spans because one
prism cannot hold a stepped solid. `union_prisms_stepped` returns the
band decomposition instead: the result is a list of constant-section
prisms, which is what the shape actually is.

The cut heights are the operand bounds, deduplicated within tolerance so
near-equal heights give one cut rather than a sliver band no solid can
carry. Verified exhaustively over 34 span configurations: bands tile the
full height with no gap and no overlap. Volume matches inclusion-
exclusion.

Mutation, 4/4 killed: a both-active band copying only the subject, cuts
not merged within tolerance, membership leaking at band boundaries, and
disjoint operands no longer refused.

The both-active mutant initially SURVIVED. The tests used a tower
standing on a slab, where the tower's section is contained in the slab's,
so 'union of both' and 'subject only' coincide. A crossing-bars case,
where neither section contains the other, kills it.

### Fillet beyond rectangles

`blend_corner` summed the two tangent offsets to find the arc centre.
That is correct only at a right angle, which is why the feature was
rectangle-only. The general corner uses the interior angle `theta`:

    setback along each edge = r / tan(theta/2)
    centre along bisector   = r / sin(theta/2)

Both reduce to the old expressions at `theta = pi/2`, so rectangles are
bit-for-bit unchanged. Verified before implementing: tangency and radius
exact to 1.1e-16 across right, obtuse, and oblique corners.

`fillet_polygon_corner` takes a ring directly, so an L-shape, a hexagon,
or a boolean result can be filleted. Refusals: a reflex corner (the arc
would bulge into the material, a different surface), a radius whose
setback overruns an adjacent edge (would swallow a neighbouring corner),
and a degenerate corner with no bisector.

Mutation, 4/4 killed. The right-angle-shortcut mutant initially
SURVIVED: the tangency test checked only the arc centre, and a wrong
setback still leaves the centre on the bisector at the correct distance
while the arc meets the walls in the wrong place. Asserting the tangent
point positions -- against the closed-form setback -- kills it.

Still open: variable-radius fillets, filleting several corners at once,
and fillets on arc-bounded profiles.
