# axiolid-pointcloud

Point-sampled geometry: an ordered set of 3D points with optional normal,
colour and intensity channels, each validated to have exactly one entry per
point. A pointcloud is a sample of a surface, with no topology or
adjacency. It is not a file format (LAS, E57 and similar are parsed outside
the kernel) and not an algorithm: queries live in `axiolid-spatial`
and reconstruction behind `axiolid-pointcloud-reconstruction-contract`. See
ADR 0044.

```bash
cargo add axiolid-pointcloud
```

- API documentation: [docs.rs/axiolid-pointcloud](https://docs.rs/axiolid-pointcloud)
- Reference page: [axiolid.github.io/kernel](https://axiolid.github.io/kernel/reference/crates/axiolid-pointcloud)
- Source and issues: [github.com/axiolid/kernel](https://github.com/axiolid/kernel)
