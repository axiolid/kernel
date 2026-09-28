# axiolid-curve

Exact, format-neutral curve values: lines and polylines (from
`axiolid-linear`), conics, rational and polynomial B-splines,
natural-equation (intrinsic) curves, elevated alignment curves, and the
curves where surfaces meet. Knots, multiplicities, weights and domains are
kept as authored. It declares the `CurveEvaluator` seam but evaluates
nothing; `axiolid-evaluate` does that. Composite, trimmed, offset and
surface-bound curves are relations in `axiolid-model`, which keeps curves
and surfaces free of a dependency cycle.

```bash
cargo add axiolid-curve
```

- API documentation: [docs.rs/axiolid-curve](https://docs.rs/axiolid-curve)
- Source and issues: [github.com/axiolid/kernel](https://github.com/axiolid/kernel)
