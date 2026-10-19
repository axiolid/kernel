# axiolid-overlay

Validated, deterministic planar booleans (intersection, union, difference, xor) and offsets over regions with holes, plus the planar operations built on them: arc-aware booleans and arrangements, polyline strokes, Minkowski morphology bounds, minimum enclosing circles and rectangles, visibility, and exact intersection of many segments by a sweep. Inputs are validated and refused with a typed error rather than repaired. It answers a query and keeps no structure; editable subdivisions with persistent identity live in `axiolid-arrangement`.

```bash
cargo add axiolid-overlay
```

- API documentation: [docs.rs/axiolid-overlay](https://docs.rs/axiolid-overlay)
- Reference page: [axiolid.github.io/kernel](https://axiolid.github.io/kernel/reference/crates/axiolid-overlay)
- Source and issues: [github.com/axiolid/kernel](https://github.com/axiolid/kernel)

## Design notes

- Straight-edged booleans (`Region`, `overlay`, `union_soup`) and arc-aware ones (`arc_overlay`,
  `ArcArrangement`) share one exact core in `src/exact_arc.rs` (ADR 0070, #173): every
  topological decision is an exact sign, and output is rounded once, so an input vertex comes
  back bit-identical. `i_overlay` remains only for offsets. The core's maintenance rules and
  verification commands are in that module's docs.
- Operand validation and output settling share one exact notion of when two boundary edges
  meet (`src/contact.rs`, #272): an end on the other edge within the tolerance, or a proper
  crossing, so every output is a valid operand wherever it lies.
  `scripts/probe_edge_contact_mutants.py` lists the faults its tests must catch.
- A degenerate input ring refuses the operation by default. `overlay_with` and
  `union_soup_with` with `OverlayOptions::dropping_degenerate_rings` leave out a ring that
  is too short, has a repeated vertex, or has zero area under the exact `ZeroArea` rule,
  and report each one with its place, reason and area (`src/degenerate.rs`, #256). A
  self-intersecting ring is never dropped. `scripts/probe_degenerate_drop_mutants.py` lists
  the faults its tests must catch.
- `segment_intersections` (`src/segment_sweep.rs`, #146) is a Bentley-Ottmann sweep over
  many segments with its own exact rational points; its reporting rule and degeneracy
  handling are in that module's docs, and `scripts/probe_segment_sweep_mutants.py` lists
  the faults its tests must catch.
- Minkowski sums and erosions take one hull per boundary edge. Exactly collinear corners
  are merged first, and a chord proved exactly to lie in the summed set trims the hulls it
  spans to what its own hull does not hold, so a wall of many short edges costs about its
  teeth, not their square, and the union stays the same set (`src/cover.rs`, #323).
  A chord whose hull or any trimmed hull would be a ring of zero area is not used, so
  trimming never adds a ring the arrangement refuses (#331). `scripts/probe_minkowski_cover_mutants.py`
  and `scripts/probe_minkowski_degenerate_mutants.py` list the faults its tests must catch.
- The disc morphology has two modes (`MorphologyMode`, ADR 0087, #330). `Exact`, the
  default and the behaviour of `erode_inner` and its siblings, is bit for bit the route
  above; a CAD kernel takes it. Its cost still grows with the square of the corners where
  teeth stand closer than their tips' slivers are wide. `Fast`, through
  `erode_inner_with` and its siblings, first fills the summed set's notches with chords
  proved to lie outside it, each used only when the excess it adds is within the caller's
  budget, and keeps the inner or outer side: the excess goes into `MorphologyBound`'s
  deviation and is reported by `Region::excess`. A wall of teeth then costs about its
  corners; a rule engine takes it (`src/fill.rs`).
  `scripts/probe_fast_morphology_mutants.py` lists the faults its tests must catch.
