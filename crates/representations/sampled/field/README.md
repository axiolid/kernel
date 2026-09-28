# axiolid-field

Frame-neutral, deterministic layered spatial-field values: row-major cells
in an explicit frame, with surface crossings and positive-length occupancy
spans kept in separate channels so a zero-thickness facet never reads as
filled space. It owns the values, their validation and caller-supplied
configuration and budgets. Sampling, morphology, clearance and navigation
live in `axiolid-field-ops`.

```bash
cargo add axiolid-field
```

- API documentation: [docs.rs/axiolid-field](https://docs.rs/axiolid-field)
- Source and issues: [github.com/axiolid/kernel](https://github.com/axiolid/kernel)
