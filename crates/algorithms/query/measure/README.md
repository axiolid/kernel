# axiolid-measure

Metric properties of geometry: surface area, signed volume, centroids and second moments of triangle meshes, closest points and distances between segments, triangles and meshes, certified Hausdorff distance between meshes, Frechet distance between polylines, winding numbers, and the area and extent of the faces selected by their normal (a floor, a wall side). Undefined quantities are refused: an open or non-manifold mesh gets an error, not a plausible volume. `cycle_volume_properties` also accepts meshes closed up to shells that touch along an edge, which are not two-manifold, and still refuses open surfaces. The optional `exact` feature adds mass properties, certified boundary distance and certified boundary Hausdorff distance for exact B-reps, single solids or bodies of several, without imposing them on mesh-only consumers.

```bash
cargo add axiolid-measure
```

- API documentation: [docs.rs/axiolid-measure](https://docs.rs/axiolid-measure)
- Reference page: [axiolid.github.io/kernel](https://axiolid.github.io/kernel/reference/crates/axiolid-measure)
- Source and issues: [github.com/axiolid/kernel](https://github.com/axiolid/kernel)
