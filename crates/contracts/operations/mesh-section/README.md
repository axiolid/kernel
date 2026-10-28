# axiolid-mesh-section-contract

The portable contract for cutting a triangle mesh with a plane: limits,
section contours, evidence, and conformance suites. It computes nothing
itself; providers implement `MeshPlaneSection`, and `axiolid-dispatch`
selects one. `section` returns closed contours and refuses a face lying in the
plane (ADR 0033). `section_request` is the opt-in alternative: a face in the
plane comes back as a region with a side, cut curves are transversal chains
that may end on a region, and touching edges are counted (ADR 0091).

```bash
cargo add axiolid-mesh-section-contract
```

- API documentation: [docs.rs/axiolid-mesh-section-contract](https://docs.rs/axiolid-mesh-section-contract)
- Reference page: [axiolid.github.io/kernel](https://axiolid.github.io/kernel/reference/crates/axiolid-mesh-section-contract)
- Source and issues: [github.com/axiolid/kernel](https://github.com/axiolid/kernel)
