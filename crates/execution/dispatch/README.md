# axiolid-dispatch

Runtime registries for operation providers: `MeshBooleanRegistry`,
`MeshPlaneSectionRegistry` and `PointcloudReconstructionRegistry`, each
behind its own feature. A registry owns provider ordering, device matching,
fallback to the next provider, and memory-budget admission. It defines no
request or result types: those live in the operation-contract crates, and
providers implement them without depending on this crate. The `parallel`
feature scopes each dispatched call to a caller-owned CPU pool.

```bash
cargo add axiolid-dispatch --features mesh-boolean
```

- API documentation: [docs.rs/axiolid-dispatch](https://docs.rs/axiolid-dispatch)
- Reference page: [axiolid.github.io/kernel](https://axiolid.github.io/kernel/reference/crates/axiolid-dispatch)
- Source and issues: [github.com/axiolid/kernel](https://github.com/axiolid/kernel)
