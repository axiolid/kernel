# axiolid-field implementation plan

Status: coverage, representation, morphology, clearance, and opt-in traversal
are implemented and mutation-verified. This is planning context, not standing
agent instruction.

## Standing invariants

- Explicit `Frame3` / `FieldBounds` / cell size / `Tolerance` /
  `FieldResourceBudget` configuration, validated together.
- Two-channel `LayeredCell`: ordered `SurfaceHit`s and disjoint occupancy.
- Deterministic scalar CPU triangle coverage emitting surface hits only.
- Separate closed-shell occupancy derivation with typed failure evidence.
- Planar masks, metric dilation/erosion, connected components.
- Clearance queries along the local layering axis.
- Geometry-only traversal behind the `navigation` feature.
- 10/10 mutation probes killed by `scripts/probe_field_gate.py`.

## Deliberately not done

**No GPU provider.** The constraint is explicit: a GPU batch provider requires
representative benchmarks demonstrating a meaningful advantage over the CPU
provider on an agreed workload. No such workload or threshold has been agreed,
so adding GPU code now would be unjustified complexity. The CPU sampler is
structured for it — per-cell independent work over a flat triangle slice — but
the seam stays unimplemented until evidence exists.

**No navigation promotion to a shared contract.** The rule is at least two
consumers needing the same neutral contract. One exists today. Until a second
appears, `navigation` stays an opt-in feature of this crate rather than a
kernel-level trait.

## Shape of the work

Only when driven by a real consumer:

- BVH-accelerated candidate triangle lookup for large scenes (measure first).
- Parallel per-row sampling via an explicit context-local pool.
- Sparse cell storage if profiling shows empty-cell waste dominates.

## Exit evidence

Targeted tests, feature-isolated compile, mutation-verified gates, and
benchmarks before any performance claim.
