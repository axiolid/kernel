# axiolid-project

Projection of triangle meshes onto a plane and intersection with prisms: the bridge between the kernel's 3D meshes and its planar booleans. `project_mesh` folds a mesh onto a plane, unions the result and keeps holes, and reports how many edge-on triangles it dropped instead of hiding them. It computes geometry, not a footprint: choosing the mesh, the reference plane and what to include is left to the consumer (ADR 0066).

```bash
cargo add axiolid-project
```

- API documentation: [docs.rs/axiolid-project](https://docs.rs/axiolid-project)
- Reference page: [axiolid.github.io/kernel](https://axiolid.github.io/kernel/reference/crates/axiolid-project)
- Source and issues: [github.com/axiolid/kernel](https://github.com/axiolid/kernel)
