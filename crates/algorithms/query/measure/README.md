# axiolid-measure

Metric properties of geometry: surface area, signed volume, centroids and second moments of triangle meshes, closest points and distances between segments, triangles and meshes, certified Hausdorff distance between meshes, Frechet distance between polylines, and winding numbers. Undefined quantities are refused: an open or non-manifold mesh gets an error, not a plausible volume. The optional `exact` feature adds mass properties and certified boundary distance for exact B-reps without imposing them on mesh-only consumers.

```bash
cargo add axiolid-measure
```

- API documentation: [docs.rs/axiolid-measure](https://docs.rs/axiolid-measure)
- Reference page: [axiolid.github.io/kernel](https://axiolid.github.io/kernel/reference/crates/axiolid-measure)
- Source and issues: [github.com/axiolid/kernel](https://github.com/axiolid/kernel)
