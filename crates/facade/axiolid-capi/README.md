# axiolid-capi

A versioned C ABI over the `axiolid::application` facade, for C and C++
applications. Every symbol carries the `axiolid_v0_4_` prefix; results are
reached through scalar handles owned by a context, data is copied into
caller-sized buffers so no Rust allocation crosses the boundary, and no
function unwinds. Exact and triangle-mesh results are distinguishable, and
an unsupported exact operation is refused rather than tessellated. The C
header is generated from the Rust surface into `include/axiolid.h`. This is
the only Axiolid crate that contains `unsafe` code. See ADR 0040.

```bash
cargo add axiolid-capi
```

- API documentation: [docs.rs/axiolid-capi](https://docs.rs/axiolid-capi)
- Reference page: [axiolid.github.io/kernel](https://axiolid.github.io/kernel/reference/crates/axiolid-capi)
- Source and issues: [github.com/axiolid/kernel](https://github.com/axiolid/kernel)
