# axiolid-refine

Mesh refinement and Laplacian smoothing with bounded, reported deviation.
Refinement splits triangles; when the source surface of a tessellated
B-rep is supplied, each new vertex is placed on that surface instead of at
the edge midpoint, so refinement converges on the real geometry rather than
subdividing the facets. Smoothing keeps boundary vertices bit-identical by
default. It does not reduce triangle counts (see `axiolid-decimate`) and
does not implement limit-surface subdivision schemes.

```bash
cargo add axiolid-refine
```

- API documentation: [docs.rs/axiolid-refine](https://docs.rs/axiolid-refine)
- Reference page: [axiolid.github.io/kernel](https://axiolid.github.io/kernel/reference/crates/axiolid-refine)
- Source and issues: [github.com/axiolid/kernel](https://github.com/axiolid/kernel)
