# axiolid-curve

Exact, format-neutral curve values: lines and polylines (from
`axiolid-linear`), conics, rational and polynomial B-splines,
natural-equation (intrinsic) curves, arc-length chains of pieces placed
end to end, elevated alignment curves, banked
centrelines carrying a cant law under a named convention, and the
curves where surfaces meet. Knots, multiplicities, weights and domains are
kept as authored. It declares the `CurveEvaluator` seam but evaluates
nothing; `axiolid-evaluate` does that. Composite, trimmed, offset and
surface-bound curves are relations in `axiolid-model`, which keeps curves
and surfaces free of a dependency cycle; `CurvePath` carries the pieces of
atomic curves such a relation runs along, reversed or placed, or offset
beside one (`PathOffset`), as a neutral value an evaluator measures.

The optional `serde` feature (off by default) derives `Serialize` and
`Deserialize` for the curve values the geometry graph's wire format carries
(ADR 0085). The format, its version and its refusal rules belong to
`axiolid-model`; this crate only supplies the derives.

```bash
cargo add axiolid-curve
```

- API documentation: [docs.rs/axiolid-curve](https://docs.rs/axiolid-curve)
- Reference page: [axiolid.github.io/kernel](https://axiolid.github.io/kernel/reference/crates/axiolid-curve)
- Source and issues: [github.com/axiolid/kernel](https://github.com/axiolid/kernel)
