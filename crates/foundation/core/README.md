# axiolid-core

The dependency root of Axiolid: points, vectors, frames, transforms,
intervals, bounding boxes, simple 2D and 3D primitives, and the explicit
`Tolerance` policy every tolerance-sensitive operation takes. It holds data
only. There are no algorithms, no serialization format, no source-format
identifiers, and no hardware backends here, and it depends on no other
Axiolid package. Coordinates are `f64` in whatever length unit the caller's
model uses.

The optional `serde` feature (off by default) derives `Serialize` and
`Deserialize` for the values the geometry graph's wire format carries
(ADR 0085). The format, its version and its refusal rules belong to
`axiolid-model`; this crate only supplies the derives.

```bash
cargo add axiolid-core
```

- API documentation: [docs.rs/axiolid-core](https://docs.rs/axiolid-core)
- Reference page: [axiolid.github.io/kernel](https://axiolid.github.io/kernel/reference/crates/axiolid-core)
- Source and issues: [github.com/axiolid/kernel](https://github.com/axiolid/kernel)
