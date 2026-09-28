# axiolid-linear

Format-neutral linear values: lines, rays, segments and polylines in 2D and
3D. It holds data only (no evaluation, tolerance policy or algorithms) and
depends only on `axiolid-core`, so an application that needs lines alone
does not compile curves, surfaces, meshes or topology. `axiolid-curve`
re-exports these types unchanged; use this crate when lines are all you
need.

```bash
cargo add axiolid-linear
```

- API documentation: [docs.rs/axiolid-linear](https://docs.rs/axiolid-linear)
- Reference page: [axiolid.github.io/kernel](https://axiolid.github.io/kernel/reference/crates/axiolid-linear)
- Source and issues: [github.com/axiolid/kernel](https://github.com/axiolid/kernel)
