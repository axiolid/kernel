# Selecting a package

The facade is convenient; leaf packages are the enforceable boundaries. Pick the smallest set that covers the use case. The [crate reference](/reference/) lists every package by layer.

- Core scalar, vector, transform and tolerance values: [`axiolid-core`](/reference/crates/axiolid-core).
- Mesh or sampled-field values without algorithms: [`axiolid-mesh`](/reference/crates/axiolid-mesh) or [`axiolid-field`](/reference/crates/axiolid-field).
- Sampling, morphology and navigation over fields: [`axiolid-field-ops`](/reference/crates/axiolid-field-ops).
- Certified exact-arithmetic predicates: [`axiolid-predicates`](/reference/crates/axiolid-predicates), the focused substrate.
- Broad reference oracles: [`axiolid-reference`](/reference/crates/axiolid-reference). It is a convenience umbrella; a narrow package does not depend on it.
- Linear values without the curve aggregate: [`axiolid-linear`](/reference/crates/axiolid-linear).
- NURBS analysis and exact shape-preserving transformations: [`axiolid-nurbs`](/reference/crates/axiolid-nurbs).
- Neutral authored graph storage: [`axiolid-model`](/reference/crates/axiolid-model).
- Exact analytic B-rep results: [`axiolid-brep`](/reference/crates/axiolid-brep). It does not tessellate.
- Mesh Boolean or plane-section portability: depend on the operation contract, then choose a provider or dispatch policy explicitly.
- Graph-to-mesh execution: [`axiolid-mesh-compile-contract`](/reference/crates/axiolid-mesh-compile-contract) for the seam, [`axiolid-mesh-compile`](/reference/crates/axiolid-mesh-compile) for the reference implementation.

Representation-only facade use:

```bash
cargo add axiolid --no-default-features --features model
```

This must not resolve compiler, field-operation, mesh-Boolean-provider, source-format or GPU dependencies. The [closure profiles](/architecture/closure-profiles) are the gate-checked examples of narrow dependency sets.
