# axiolid-mesh-section-contract

The portable contract for cutting a triangle mesh with a plane: limits,
section contours, evidence, and a conformance suite. It computes nothing
itself; providers implement `MeshPlaneSection`, and `axiolid-dispatch`
selects one. See ADR 0033.

```bash
cargo add axiolid-mesh-section-contract
```

- API documentation: [docs.rs/axiolid-mesh-section-contract](https://docs.rs/axiolid-mesh-section-contract)
- Source and issues: [github.com/axiolid/kernel](https://github.com/axiolid/kernel)
