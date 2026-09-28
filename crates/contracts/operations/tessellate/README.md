# axiolid-tessellation-contract

The contract for turning exact geometry from an `axiolid-model` graph into
triangles under an explicit tolerance: `TessellationOptions` (which has no
default chord error), the `Tessellator` trait, and a `TessellatedMesh` that
carries the tolerance it was built to. Adjacent faces must share one
discretisation of each topological edge, because tessellating faces
independently is not watertight. It is a contract only; providers implement
it.

```bash
cargo add axiolid-tessellation-contract
```

- API documentation: [docs.rs/axiolid-tessellation-contract](https://docs.rs/axiolid-tessellation-contract)
- Source and issues: [github.com/axiolid/kernel](https://github.com/axiolid/kernel)
