# axiolid-nurbs

Format-neutral algorithms over polynomial and rational B-spline curves and
surfaces: differential geometry, exact shape-preserving transforms (knot
insertion, reversal, splitting, Bezier decomposition, degree elevation),
tolerance-bounded knot removal and degree reduction, interpolation,
least-squares curve and surface approximation, lofting, certified
projection, inversion and intersection queries, exact analytic curve and
surface intersection, and verified periodic seams.
Lossy operations measure their deviation and refuse above the caller's
tolerance. It owns no importer, tessellator or file-format vocabulary, and
it uses `axiolid-evaluate` for evaluation rather than reimplementing it.

```bash
cargo add axiolid-nurbs
```

- API documentation: [docs.rs/axiolid-nurbs](https://docs.rs/axiolid-nurbs)
- Reference page: [axiolid.github.io/kernel](https://axiolid.github.io/kernel/reference/crates/axiolid-nurbs)
- Source and issues: [github.com/axiolid/kernel](https://github.com/axiolid/kernel)
