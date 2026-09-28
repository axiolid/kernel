# axiolid-minkowski

Minkowski sum and difference of closed planar-faced (triangle-mesh) solids.
The sum of two convex solids is computed exactly as the hull of pairwise
vertex sums; non-convex operands are decomposed into convex parts and the
pairwise sums unioned through a caller-supplied mesh Boolean provider,
under a budget. The difference is computed as an erosion, not as a hull of
pairwise differences, and refuses a non-convex subject rather than return a
result that is too large. Curved operands are refused.

```bash
cargo add axiolid-minkowski
```

- API documentation: [docs.rs/axiolid-minkowski](https://docs.rs/axiolid-minkowski)
- Source and issues: [github.com/axiolid/kernel](https://github.com/axiolid/kernel)
