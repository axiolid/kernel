# axiolid-core

The dependency root of Axiolid: points, vectors, frames, transforms,
intervals, bounding boxes, simple 2D and 3D primitives, and the explicit
`Tolerance` policy every tolerance-sensitive operation takes. It holds data
only. There are no algorithms, no serialization, no source-format
identifiers, and no hardware backends here, and it depends on no other
Axiolid package. Coordinates are `f64` in whatever length unit the caller's
model uses.

```bash
cargo add axiolid-core
```

- API documentation: [docs.rs/axiolid-core](https://docs.rs/axiolid-core)
- Reference page: [axiolid.github.io/kernel](https://axiolid.github.io/kernel/reference/crates/axiolid-core)
- Source and issues: [github.com/axiolid/kernel](https://github.com/axiolid/kernel)
