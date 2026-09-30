# axiolid-spatial

Deterministic, callback-based spatial acceleration: a median-split BVH over bounded objects and a uniform grid for point KNN and radius search, both behind the `SpatialIndex` query contract. They return candidates only, never exact intersections. It also provides barycentric, mean-value, Wachspress and discrete harmonic coordinates for interpolating values given at triangle, tetrahedron, polygon and closed-triangle-mesh corners.

```bash
cargo add axiolid-spatial
```

- API documentation: [docs.rs/axiolid-spatial](https://docs.rs/axiolid-spatial)
- Reference page: [axiolid.github.io/kernel](https://axiolid.github.io/kernel/reference/crates/axiolid-spatial)
- Source and issues: [github.com/axiolid/kernel](https://github.com/axiolid/kernel)
