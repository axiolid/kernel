# axiolid-decimate

Edge-collapse decimation of triangle meshes with a bounded, reported
deviation. The caller asks for a triangle budget or a maximum deviation;
either way the result never moves a vertex further than the caller's
bound, and `DecimateReport` states the collapses performed, the refusals
by cause and the largest distance any vertex actually moved. Collapses that
would invert a triangle or create a non-manifold edge are refused. Output
is deterministic. It does not remesh isotropically, detect sharp features
or use quadric error metrics: the cost is edge length and the new vertex
is the edge midpoint. For adding triangles instead, see `axiolid-refine`.

```bash
cargo add axiolid-decimate
```

- API documentation: [docs.rs/axiolid-decimate](https://docs.rs/axiolid-decimate)
- Source and issues: [github.com/axiolid/kernel](https://github.com/axiolid/kernel)
