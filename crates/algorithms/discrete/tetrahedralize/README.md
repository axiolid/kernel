# axiolid-tetrahedralize

Exact 3D Delaunay tetrahedralization. `Delaunay3` inserts points one at a
time (Bowyer-Watson, located by a visibility walk) and keeps the convex hull
with an infinite vertex. Every orientation and in-sphere decision is an
exact sign from `axiolid-predicates`, and ties -- cospherical points,
points in the plane of a hull face -- are broken by a symbolic perturbation
that depends only on the coordinates, so a point set has one result
whatever the insertion order. Duplicates are merged; inputs spanning fewer
than three dimensions are refused. Constraints (input segments and faces)
are not recovered.

```bash
cargo add axiolid-tetrahedralize
```

- API documentation: [docs.rs/axiolid-tetrahedralize](https://docs.rs/axiolid-tetrahedralize)
- Reference page: [axiolid.github.io/kernel](https://axiolid.github.io/kernel/reference/crates/axiolid-tetrahedralize)
- Source and issues: [github.com/axiolid/kernel](https://github.com/axiolid/kernel)
