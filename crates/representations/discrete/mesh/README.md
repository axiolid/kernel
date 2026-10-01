# axiolid-mesh

Triangle and polygon mesh values: `TriMesh` as the compact exchange type,
`PolygonMesh` keeping n-gons and holes until explicit triangulation, u32
indices, named per-vertex and per-corner attribute channels, derived edge
adjacency, connected components, and a deterministic structural audit that
reports defects instead of rejecting dirty input. `HalfedgeMesh` is the
editable surface mesh for algorithms that rewrite connectivity in place:
O(1) adjacency, vertex circulators, boundary loops, and edge flip, split and
collapse (with the link condition) and face split that keep it a valid
oriented 2-manifold; non-manifold input is refused by name. `MeshView` and
`TriangleMeshView` let a foreign mesh be read without copying it. Mesh
operations such as booleans, sections and repair live in other crates.

```bash
cargo add axiolid-mesh
```

- API documentation: [docs.rs/axiolid-mesh](https://docs.rs/axiolid-mesh)
- Reference page: [axiolid.github.io/kernel](https://axiolid.github.io/kernel/reference/crates/axiolid-mesh)
- Source and issues: [github.com/axiolid/kernel](https://github.com/axiolid/kernel)

## Design notes

Rendering appearance (materials, shaders, colours as styling) is not part of a
mesh value. Data that must follow the geometry, such as a material id per
vertex, travels as an attribute channel.
