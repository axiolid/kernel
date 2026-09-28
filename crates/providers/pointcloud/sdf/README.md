# axiolid-pointcloud-reconstruction-sdf

The reference `PointcloudReconstruction` provider. `SdfReconstruction` builds
a signed-distance field from the samples (nearest neighbours through
`axiolid-spatial`) and extracts its zero level set with `axiolid-levelset`.
With normals the surface passes through the samples; without them it wraps
around them, and the evidence says which. It is not a hole filler: where the
capture has no data the surface is extrapolated, and those triangles are
counted. It has no external dependency, so a better reconstruction can
replace it behind the same contract.

```bash
cargo add axiolid-pointcloud-reconstruction-sdf
```

- API documentation: [docs.rs/axiolid-pointcloud-reconstruction-sdf](https://docs.rs/axiolid-pointcloud-reconstruction-sdf)
- Reference page: [axiolid.github.io/kernel](https://axiolid.github.io/kernel/reference/crates/axiolid-pointcloud-reconstruction-sdf)
- Source and issues: [github.com/axiolid/kernel](https://github.com/axiolid/kernel)
