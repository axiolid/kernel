# axiolid-fixtures

A shared corpus of adversarial and degenerate mesh fixtures (sliver and
zero-area triangles, open shells, extreme scale disparity, coplanar contact
between operands, and the like), each with a `Provenance` saying where the
case comes from, its licence, and what an implementation must do with it.
Fixtures are built in code rather than stored as files, so the exact bit
patterns that make them degenerate cannot be rounded away by an exporter.
Differential tests, for example in `axiolid-mesh-boolean-boolmesh`, iterate
`corpus()`. Every fixture is original work under the repository licence.

```bash
cargo add axiolid-fixtures
```

- API documentation: [docs.rs/axiolid-fixtures](https://docs.rs/axiolid-fixtures)
- Reference page: [axiolid.github.io/kernel](https://axiolid.github.io/kernel/reference/crates/axiolid-fixtures)
- Source and issues: [github.com/axiolid/kernel](https://github.com/axiolid/kernel)
