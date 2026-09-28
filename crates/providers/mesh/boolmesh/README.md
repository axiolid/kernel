# axiolid-mesh-boolean-boolmesh

A `MeshBoolean` provider for closed, outward-oriented triangle meshes: union,
intersection, difference, and symmetric difference composed from those. The
algorithm is absorbed from the `boolmesh` crate into a private module (ADR
0014, ADR 0047); this crate adds the conversion, an orientation gate on every
input, result checks, and batch overrides (`subtract_many` fuses disjoint
cutters, `union_many` reduces as a balanced tree). It also offers
`subtract_boxes_analytic`, an opt-in closed-form path for axis-aligned box
cutters in an axis-aligned box. Register it with `axiolid-dispatch` or call it
directly.

```bash
cargo add axiolid-mesh-boolean-boolmesh
```

- API documentation: [docs.rs/axiolid-mesh-boolean-boolmesh](https://docs.rs/axiolid-mesh-boolean-boolmesh)
- Reference page: [axiolid.github.io/kernel](https://axiolid.github.io/kernel/reference/crates/axiolid-mesh-boolean-boolmesh)
- Source and issues: [github.com/axiolid/kernel](https://github.com/axiolid/kernel)

## Design notes

- Results carry attribute channels but no normals; derive normals from the
  topology you want.
- The general path is `Determinism::Topological`. Use the analytic box path
  when you need byte-identical output across processes.
- The `parallel` feature threads inside one solve and is off by default;
  `parallel-batch` runs independent `union_many` pairs concurrently with
  output identical to the sequential path.
