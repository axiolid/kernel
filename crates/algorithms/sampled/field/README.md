# axiolid-field-ops

Deterministic algorithms over `axiolid-field` layered fields: scalar CPU triangle coverage, planar masks with metric dilation, erosion and connected components, clearance along the layering axis, and, behind the opt-in `navigation` feature, geometry-only route finding under an explicit agent envelope. It reports coverage, spans, components and route existence, never an application verdict such as accessibility or compliance.

```bash
cargo add axiolid-field-ops
```

- API documentation: [docs.rs/axiolid-field-ops](https://docs.rs/axiolid-field-ops)
- Reference page: [axiolid.github.io/kernel](https://axiolid.github.io/kernel/reference/crates/axiolid-field-ops)
- Source and issues: [github.com/axiolid/kernel](https://github.com/axiolid/kernel)

## Design notes

There is no GPU coverage provider. The CPU sampler is per-cell independent work over a flat
triangle slice, so a batch provider can be added once a benchmark on an agreed workload shows it
pays; without that evidence it would be complexity with no measured benefit.
