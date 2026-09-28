# axiolid-overlay

Validated, deterministic planar booleans (intersection, union, difference, xor) and offsets over regions with holes, plus the planar operations built on them: arc-aware booleans and arrangements, polyline strokes, Minkowski morphology bounds, minimum enclosing circles and rectangles, and visibility. Inputs are validated and refused with a typed error rather than repaired. It answers a query and keeps no structure; editable subdivisions with persistent identity live in `axiolid-arrangement`.

```bash
cargo add axiolid-overlay
```

- API documentation: [docs.rs/axiolid-overlay](https://docs.rs/axiolid-overlay)
- Source and issues: [github.com/axiolid/kernel](https://github.com/axiolid/kernel)

## Design notes

- Straight-edged booleans (`Region`, `overlay`, `union_soup`) and arc-aware ones (`arc_overlay`,
  `ArcArrangement`) share one exact core in `src/exact_arc.rs` (ADR 0070, #173): every
  topological decision is an exact sign, and output is rounded once, so an input vertex comes
  back bit-identical. `i_overlay` remains only for offsets. The core's maintenance rules and
  verification commands are in that module's docs.
