# axiolid-construct

Solid generation from exact inputs: extrusion, revolution, sweeps, lofts,
centre-line profiles, half-space clipping proxies, offsets, fillets and
chamfers on supported families, and focused exact booleans (planar
polyhedra, coaxial column solids). Every `Profile` variant extrudes to an
exact B-rep, and full-turn revolution covers any profile that lowers to a
contour; sweeps and lofts produce meshes by default. Geometry the kernel cannot
represent exactly is refused, never tessellated in its place. The crate
takes geometry and returns geometry: it owns no operation graph, cache,
execution context or provider dispatch (ADR 0023); `axiolid-mesh-compile`
does those and calls in here.

```bash
cargo add axiolid-construct
```

- API documentation: [docs.rs/axiolid-construct](https://docs.rs/axiolid-construct)
- Reference page: [axiolid.github.io/kernel](https://axiolid.github.io/kernel/reference/crates/axiolid-construct)
- Source and issues: [github.com/axiolid/kernel](https://github.com/axiolid/kernel)

## Design notes

Tests that check a generated mesh is accepted by a mesh Boolean provider
live in `crates/execution/compile/tests/`, not here. A dev-dependency on a
provider or execution crate would pull this algorithms crate's tests above
its tier; the architecture check only enforces the tier edge for normal
dependencies, so the allowlist in `Cargo.toml` is what keeps it out.
