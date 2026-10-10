# axiolid-model

The format-neutral geometry graph that source adapters lower into and
kernels consume: an immutable, append-only DAG of typed nodes preserving
exact curves, surfaces, profiles, primitives, topology, curve and surface
relations, CSG instructions and instancing, alongside source meshes. Handles
are branded per graph and references must point to earlier nodes of the
right family, so cycles, dangling references and cross-graph handles cannot
be built. It evaluates, tessellates and compiles nothing, and keeps source
identifiers outside the graph.

The optional `serde` feature (off by default) adds the graph's versioned
wire format (ADR 0085): `GeometryGraph::to_json`, `from_json`, `to_cbor`
and `from_cbor`, and the `wire` module. Every payload names its format and
`MAJOR.MINOR` version; additions are minor versions and anything else is
major. A reader refuses by name a version it does not read, a kind,
variant or field it does not know, a non-finite number and a graph that
fails the builder's validation, and never returns part of a payload.

```bash
cargo add axiolid-model
```

- API documentation: [docs.rs/axiolid-model](https://docs.rs/axiolid-model)
- Reference page: [axiolid.github.io/kernel](https://axiolid.github.io/kernel/reference/crates/axiolid-model)
- Source and issues: [github.com/axiolid/kernel](https://github.com/axiolid/kernel)
