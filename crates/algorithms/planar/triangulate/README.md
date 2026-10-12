# axiolid-triangulate

Constrained Delaunay triangulation with bounded quality refinement. The output covers the exact convex hull of the input, exactly coincident points are merged into one vertex and the merge is reported, every constraint edge survives as a union of output edges, the result is Delaunay away from the constraints (decided by the certified `incircle` predicate), and optional Ruppert refinement drives interior angles toward a caller-chosen minimum. Refinement carries an explicit Steiner budget and reports when it was capped, so an unmet angle bound is never returned silently.

```bash
cargo add axiolid-triangulate
```

- API documentation: [docs.rs/axiolid-triangulate](https://docs.rs/axiolid-triangulate)
- Reference page: [axiolid.github.io/kernel](https://axiolid.github.io/kernel/reference/crates/axiolid-triangulate)
- Source and issues: [github.com/axiolid/kernel](https://github.com/axiolid/kernel)
