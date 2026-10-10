# axiolid-profile

Exact 2D profiles for sweeps and sectioned solids: parameterised
rectangles, circles, ellipses and structural sections, closed contours of
bounded exact curve segments with holes, centre-line profiles, and
transformed or composite profiles. It stores profile intent only. Boolean
cleanup, offsetting and triangulation are algorithms in higher tiers, so a
consumer can read profiles without them.

The optional `serde` feature (off by default) derives `Serialize` and
`Deserialize` for the profile values the geometry graph's wire format carries
(ADR 0085). The format, its version and its refusal rules belong to
`axiolid-model`; this crate only supplies the derives.

```bash
cargo add axiolid-profile
```

- API documentation: [docs.rs/axiolid-profile](https://docs.rs/axiolid-profile)
- Reference page: [axiolid.github.io/kernel](https://axiolid.github.io/kernel/reference/crates/axiolid-profile)
- Source and issues: [github.com/axiolid/kernel](https://github.com/axiolid/kernel)
