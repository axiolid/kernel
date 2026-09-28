# axiolid-reference

The portable scalar reference implementation that every optimized backend
is differentially tested against (ADR 0012): a solid boolean oracle, an
exact-sign mesh plane-section oracle, triangle/triangle and
segment/triangle relations, clash detection, convex hulls, polygon
triangulation and tessellation. It favours readability over speed: no
intrinsics, threading, feature gates or `unsafe`. It is also a convenience
umbrella that re-exports `axiolid-predicates` and `axiolid-evaluate`
unchanged (ADR 0036); a consumer that needs only certified signs or curve
evaluation should depend on those directly.

```bash
cargo add axiolid-reference
```

- API documentation: [docs.rs/axiolid-reference](https://docs.rs/axiolid-reference)
- Reference page: [axiolid.github.io/kernel](https://axiolid.github.io/kernel/reference/crates/axiolid-reference)
- Source and issues: [github.com/axiolid/kernel](https://github.com/axiolid/kernel)
