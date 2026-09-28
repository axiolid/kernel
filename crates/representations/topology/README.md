# axiolid-topology

Exact B-rep topology with typed handles: vertices, edges, edge uses, loops,
faces, shells and solids, each with its own handle type and explicit
orientation, plus a structural audit. Geometry is linked through a
caller-chosen handle type (`BRep<G>`), so the graph does not depend on any
curve or surface model and serves exact kernels, mesh converters and import
adapters alike. `axiolid-brep` binds it to Axiolid's own curves and
surfaces.

```bash
cargo add axiolid-topology
```

- API documentation: [docs.rs/axiolid-topology](https://docs.rs/axiolid-topology)
- Reference page: [axiolid.github.io/kernel](https://axiolid.github.io/kernel/reference/crates/axiolid-topology)
- Source and issues: [github.com/axiolid/kernel](https://github.com/axiolid/kernel)
