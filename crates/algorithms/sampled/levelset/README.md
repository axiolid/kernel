# axiolid-levelset

Level-set extraction: a closed two-manifold triangle mesh from a scalar field sampled on a grid. Cells are split into Kuhn tetrahedra rather than marched as cubes, so shared faces always split the same way and the result is watertight by construction; exact grid tangency is resolved by simulation of simplicity. The mesh interpolates the field linearly along cell edges, so it is an approximation whose error shrinks with the grid, not a certified surface.

```bash
cargo add axiolid-levelset
```

- API documentation: [docs.rs/axiolid-levelset](https://docs.rs/axiolid-levelset)
- Reference page: [axiolid.github.io/kernel](https://axiolid.github.io/kernel/reference/crates/axiolid-levelset)
- Source and issues: [github.com/axiolid/kernel](https://github.com/axiolid/kernel)
