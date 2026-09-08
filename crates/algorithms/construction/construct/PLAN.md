# axiolid-construct — known limits and open work

Findings that outlived the session that produced them. Each entry records
what was measured, what was tried, and what a fix would actually involve,
so a later pass does not re-derive it or repeat a reverted approach.

## Repeated grid-aligned subtraction refuses (open)

`boolean_polyhedra_exact` refuses partway through a long chain of
axis-aligned differences. The depth-2 Menger sponge void set (147
consecutive subtractions from a unit cube) refuses at subtraction 82.
Reproduce with the benchmark harness at `axiolid/benchmarks`:
`AXIOLID_MENGER_DEPTH=3 cargo run --release -- 1`.

### What was measured

The refusal is `every probe direction met a vertex or edge exactly`. The
cause is NOT an unlucky ray: 4 directions and 12 directions both fail at
exactly step 82. The blocking face was identified as a collapsed quad --

```
(1.0, 1.0, 0.3333333333333333 )   (0.6666666666666666, 1.0, 0.333...3)
(1.0, 1.0, 0.33333333333333326)   (0.6666666666666666, 1.0, 0.333...3)
```

-- spanning a third of the model, vertices paired, the pairs ONE ULP
apart. It encloses no area, so its normal is meaningless and every ray
meets it edge-on. No probe direction can classify a face with no plane.

Such rings arise because `plane_crossing` builds intersection coordinates
in f64 (ADR 0045). A vertex that should be shared between operands lands
a few ULPs apart after several operations, and splitting through it emits
a face that corresponds to nothing in the modelled solid.

### What was tried and reverted

Commit `40b5069` dropped split fragments enclosing no area at f64
precision, and was reverted. It DID clear the blockage -- all 147
subtractions completed -- but the answers were wrong: depth-2 volume was
short by 4.57e-4 against a cell volume of 1.37e-3, roughly a third of a
cell. Deleting the collapsed faces opens holes in the shell, and an open
shell integrates to a wrong volume.

That is worse than the refusal it replaced. A refusal is actionable; a
plausible wrong volume is silent. **Do not reintroduce a drop-based fix.**

### Why exact construction is not the answer either

ADR 0045 declines exact constructions, and the benchmark data supports
that for this case. Against CGAL's exact-construction kernel on the
thin-overlap sweep, Axiolid tracks within ~1.5x down to 1e-12, and at
1e-15 both are catastrophically wrong because the ambiguity is in the
input representation rather than the arithmetic.

### The actual shape of a fix

The collapsed rings are a SYMPTOM. Two candidate directions, neither
attempted:

1. Make coincident split points exactly coincident, so the drift never
   arises. This is snapping, which ADR 0045 rejects by name -- it would
   need a superseding ADR.
2. Refuse on a collapsed ring instead of deleting it, keeping the
   diagnosis without the wrong answer. Strictly better than the current
   state, and does not touch ADR 0045.

`boolean_polyhedra_exact` has no production consumers today, so neither
is urgent. The blast radius of `plane_crossing` is one call site.

## Thin overlap below 1e-12 produces an open shell (open)

Found by the contact matrix, not by the sponge. Two unit boxes
overlapping by `eps` along +X:

```
eps 1e-3 : 14 faces, 28 tris, 28 usable, 0 degenerate, 0 boundary edges
eps 1e-6 : same
eps 1e-9 : same
eps 1e-12: 14 faces, 28 tris, 20 usable, 8 degenerate, 8 boundary edges
eps 1e-15: same
```

From 1e-12 the union and intersection come back with a hole, so they
cannot be measured. This is NOT a precision floor: a 1e-12 slab on unit
boxes is ~4503 ULPs wide, comfortably resolvable in f64. It is a defect
in how a very thin overlap is split.

`tests/contact_matrix.rs` asserts this limit explicitly rather than
loosening its expectations. When it is fixed, delete the
`Outcome::Unmeasurable(_) if eps <= 1e-12` branches and the sweep will
hold the fix in place.
