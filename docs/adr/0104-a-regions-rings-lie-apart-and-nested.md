# 0104 — A region's rings lie apart and nested, decided on the exact arrangement

- **Status:** Accepted
- **Date:** 2026-11-07
- **Deciders:** Friedrich Schrödter
- **Supersedes:** —

## Context

Issue #448, found while reviewing #345. `Region::new` validated each ring on
its own (`validate_ring`), so it accepted a hole crossing or leaving its outer
ring, holes overlapping and polygons overlapping. The crate then read such a
set three ways: the boolean by winding numbers (outer `+1`, hole `-1`, so the
part of a hole outside its outer ring counted as filled under `NonZero`), the
morphology by ring membership (`minkowski`'s `Layout`: a polygon holds a point
in its outer ring and in none of its holes), and the fast morphology's notch
fill (`fill`, ADR 0087) assumed rings that are apart and nested. Fuzz seed
21351: a toothed hole whose needles reach y = 2.2955, past its outer ring's
top edge at y = 2.1202; `dilate_outer_with(0.9046, fast(1e-2))` left out
(2.0172, 3.0256), 0.8993 from the region, with a stated deviation of 4.7e-3.
`area` summed overlapping polygons twice and `component_count` counted them
apart. A consumer (rheonda/cad `check_polygon`) already documented
`Region::new` as checking "the holes inside the outer ring and apart".

Outputs must stay inputs (#188, #191): `tests/roundtrip.rs` feeds every
operation's result back into `Region::new`. A union of subdivided cells there
returns two polygons that share an edge, run in opposite directions.

## Decision

**`Region::new` refuses a polygon set unless every hole lies inside its own
outer ring and in no other of its holes, and no point lies in two polygons,
decided exactly on one arrangement of all its rings (`src/nesting.rs`).**

- **The rule is about the point set.** Each ring is read as the region it
  encloses, whatever its winding; a polygon is its outer ring's inside less
  its holes'. Refused, by name: `HoleOutsideOuter` (part of a hole's inside
  outside its outer ring, a crossing included), `HolesOverlap` (crossing,
  nested or repeated holes of one polygon), `PolygonsOverlap` (a point in two
  polygons). Rings may touch at points, and along a stretch two rings run
  with the polygons' insides on opposite sides; both occur in settled outputs and
  every operation reads them as the point set they bound. A ring whose
  vertices all lie on one line encloses nothing and is left out (#219).
  Contacts that cut a polygon's inside apart are not refused.
- **Exact, tolerance-free.** As in `check_arc_polygon` (ADR 0095), every
  decision is an exact sign; the tolerance only validates rings one by one.
- **How.** A box-tree pass over the edges of different rings refuses a proper
  crossing whose defect the rings' roles alone name (a hole with its own
  outer ring or another hole of its polygon; the outer rings of two polygons
  without holes), so a soup given as a region is refused before its
  crossings are cut. Then all rings go into one exact arrangement
  (`exact_arc::arrangement`), and each side of each piece is read from the
  rings holding it, by the morphology's membership rule. This finds what a
  per-ring probe cannot: a crossing at a vertex, a hole leaving through two
  vertices on its outer ring's edge, a stretch two rings run the same way.
- **The constructor and `translate` check.** Booleans and the morphology return
  the faces of an exact arrangement, settled, which are apart and nested by
  construction; they are not checked again. Evidence: with the check run on
  every internal region of the crate's test suite, none failed, nor did
  3,480 outputs of the exact and fast morphology on 478 valid regions of
  #345's fuzz family. `Region::translate` is the exception: it moves every
  vertex by its own rounding, which can turn a contact into a crossing, so
  it checks its result and refuses by name. A hole sharing a vertex with its
  outer ring, outer `[(-0.9, 3.8000000000000003), (-0.5, 3.5), (0.1, 4.3),
  (-0.30000000000000004, 4.6000000000000005)]`, translated by
  `(0.0223, 0.1729)`, ends with the hole crossing the outer ring; a hole
  and an outer ring that both lie on one line (left out, so accepted)
  translate to rings that no longer do. Both are tests in
  `tests/region_nesting.rs`. `Region::sweep` translates, so it can refuse
  the same way.
- **Operands of `overlay` and `union_soup` are unchanged**: their rings may
  overlap, read by the fill rule.
- **Errors.** `OverlayError` becomes `#[non_exhaustive]` and gains
  `HolesOverlap` and `PolygonsOverlap`; axiolid-overlay 0.4.0, the facade
  0.5.0. Input accepted before is now refused, so downstream crates opt in
  by updating.

## Alternatives considered

| Option | Why not |
| --- | --- |
| Fast morphology routes such regions to the exact mode | Leaves the three readings disagreeing (the boolean's winding count, `area`, `component_count`) and accepts input the crate cannot read one way; #448 asks to fail closed. |
| OGC validity: no shared stretch between any rings | Refuses the crate's own outputs: a union of subdivided cells returns polygons sharing an edge (`tests/roundtrip.rs`). |
| One probe point per ring, located through a box tree | About 4 to 7 times cheaper, but a hole can leave its outer ring through two vertices on its edge with no proper crossing; catching that needs a side analysis at every contact, which the arrangement already does exactly. |
| Contacts decided within the tolerance | The fast fill's proofs are exact, a sub-tolerance crossing still breaks them, and ADR 0095 keeps contact decisions tolerance-free. |
| Check every output too | Costs every operation for a property their construction gives; the suite and the fuzz family found no output that fails. |

## Consequences

**Positive**

- Seed 21351 is refused as `HoleOutsideOuter`; every region an operation
  sees has one reading. `area` and `component_count` are true for inputs.
- The check names its defect; consumers validating polygons through
  `Region::new` get what they assumed.

**Negative / costs**

- `Region::new` of a region of more than one ring builds an arrangement:
  a floor of 100,000 corners with 200 holes
  takes 0.34 s instead of 0.083 s, one of 10,000 corners 0.032 s instead of
  0.007 s (`benches/ring_validation.rs`, case `holes`, median of 3).
  Rooms standing on one corridor edge cost 0.05 s for 1,000 and 0.22 s for
  4,000; a rectangle with holes standing on its bottom edge 0.05 s and
  0.23 s (cases `corridor` and `bays`). The first build of this check
  searched every edge sharing a stretch with a piece's edge, the square of
  the number standing on one edge (4,000 rooms 35 s, 4,000 holes 87 s); the
  arrangement now sweeps each edge's pieces against the edges covering
  them. A region of one ring costs what it did.
- Inputs accepted before are refused: overlapping polygons given as one
  region (axioval's crowded-sliver test fixtures) must be united first, and
  so must adjacent polygons computed in floating point (rounded room
  footprints), which are now refused as `PolygonsOverlap`.
- An input of polygons with holes that cross each other many times is cut
  in full before it is refused.

**Follow-ups / risks to watch**

- Contacts that cut a polygon's inside apart (`DisconnectedInterior` in ADR
  0095) are accepted (follow-up #467), so `component_count` can count one polygon whose
  inside is in pieces.
- Outputs are not snap-rounded (ADR 0100 covers another crate): a rounding
  that turns a contact into a crossing of an ulp would make an output fail
  this check when fed back. None was found; a failure is a refusal, not a
  wrong answer.

## Relation to existing code

- `crates/algorithms/planar/overlay/src/nesting.rs`: `check`.
- `crates/algorithms/planar/overlay/src/region.rs`: `Region::new`.
- `crates/algorithms/planar/overlay/tests/region_nesting.rs`,
  `tests/roundtrip.rs`; `scripts/probe_region_nesting_mutants.py`.
- `crates/algorithms/planar/overlay/benches/ring_validation.rs` (`holes`).
