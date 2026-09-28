# axiolid-ray-mesh

Narrow-phase ray/triangle-mesh intersection: the nearest hit with its parameter, barycentric coordinates and a certified front/back/coplanar side. It composes with a broad phase such as `axiolid-spatial` by taking candidate triangle indices, but does not depend on one. Degenerate triangles are refused rather than silently missed. It owns the intersection only, not what a ray means to the caller.

```bash
cargo add axiolid-ray-mesh
```

- API documentation: [docs.rs/axiolid-ray-mesh](https://docs.rs/axiolid-ray-mesh)
- Reference page: [axiolid.github.io/kernel](https://axiolid.github.io/kernel/reference/crates/axiolid-ray-mesh)
- Source and issues: [github.com/axiolid/kernel](https://github.com/axiolid/kernel)
