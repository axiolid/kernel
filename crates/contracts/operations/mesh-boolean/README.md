# axiolid-mesh-boolean-contract

The portable mesh-boolean contract: the `MeshBoolean` provider trait, the
evidence a provider must report with its result, and a conformance suite.
Operand admissibility comes from `axiolid-mesh-contracts`. It performs no
booleans itself; providers such as `axiolid-mesh-boolean-boolmesh`
implement it, and `axiolid-dispatch` chooses between them.

```bash
cargo add axiolid-mesh-boolean-contract
```

- API documentation: [docs.rs/axiolid-mesh-boolean-contract](https://docs.rs/axiolid-mesh-boolean-contract)
- Source and issues: [github.com/axiolid/kernel](https://github.com/axiolid/kernel)
