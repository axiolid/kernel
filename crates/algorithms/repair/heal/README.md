# axiolid-heal

Explicit diagnosis and opt-in repair of triangle meshes. `diagnose` reports non-manifold edges, inconsistent winding, boundary edges, duplicate vertices, degenerate triangles and exact self-intersections without touching the mesh. Repairs (weld vertices, drop degenerate elements, unify orientation, orient outward) run only when a caller names them in a `RepairPlan`, and the `RepairReport` records what was applied, what was skipped and what happened to each attribute channel. There is no repair-everything mode, and no other operation heals implicitly.

```bash
cargo add axiolid-heal
```

- API documentation: [docs.rs/axiolid-heal](https://docs.rs/axiolid-heal)
- Reference page: [axiolid.github.io/kernel](https://axiolid.github.io/kernel/reference/crates/axiolid-heal)
- Source and issues: [github.com/axiolid/kernel](https://github.com/axiolid/kernel)
