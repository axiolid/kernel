# axiolid-model

The format-neutral geometry graph that source adapters lower into and
kernels consume: an immutable, append-only DAG of typed nodes preserving
exact curves, surfaces, profiles, primitives, topology, curve and surface
relations, CSG instructions and instancing, alongside source meshes. Handles
are branded per graph and references must point to earlier nodes of the
right family, so cycles, dangling references and cross-graph handles cannot
be built. It evaluates, tessellates and compiles nothing, and keeps source
identifiers outside the graph.

```bash
cargo add axiolid-model
```

- API documentation: [docs.rs/axiolid-model](https://docs.rs/axiolid-model)
- Reference page: [axiolid.github.io/kernel](https://axiolid.github.io/kernel/reference/crates/axiolid-model)
- Source and issues: [github.com/axiolid/kernel](https://github.com/axiolid/kernel)
