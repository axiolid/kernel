# axiolid-overlay

Planar booleans and offsets over validated regions.

## Two backends

- Straight edges: `Region` / `overlay` on the `i_overlay` integer backend.
- Arcs: `arc_overlay` on the in-tree exact core `src/exact_arc.rs` and `src/exact_arc/`
  (ADR 0070). `point.rs` holds exact points `(a + b*sqrt(d)) / w` and their
  filtered predicates; `edge.rs` holds segments and bulge arcs (circle in
  exact dyadic conic form, half-angle parametrization for rational samples);
  `mod.rs` splits, classifies, keeps, links and nests.

## Rules

- No decision in `exact_arc` may read a tolerance. The tolerance only
  validates operands and cleans up output rounding (`presented` in
  `arc_overlay.rs`).
- Every predicate goes through `point::sign`, which tries cached boxes,
  then the `axiolid-exact` interval tier, then exact arithmetic. `Sign` is
  `non_exhaustive`: match `Positive`/`Negative` and treat the rest as zero.
- Output vertices are rounded once. Never feed rounded output back into a
  decision.

## Verification

- `tests/arc_exact_oracle.rs`: area identities and point membership against
  tessellated operands on random grid-snapped and decimal scenes.
- `python3 scripts/probe_arc_overlay_mutants.py`: every listed fault must
  fail the suite.
- `cargo bench -p axiolid-overlay --bench arc_overlay`: per-call cost;
  `SCALE=1` runs the edge-count scaling scenes instead. Quote it when
  claiming a speed change.
- Bounding boxes (`edge.rs::Bounds`) only skip work and must contain the
  whole edge. Changing how they are built needs the mutation probe: the
  oracle tests include major arcs and many-edge rings for this.
