# axiolid-mesh-compile-contract

The contract for compiling an `axiolid-model` geometry graph into a
triangle mesh. The outcome says whether the mesh is closed, and the
contract never claims to preserve an exact B-rep. It is a contract only:
`axiolid-mesh-compile` implements it, and `axiolid-exact-compile-contract`
is the exact counterpart.

```bash
cargo add axiolid-mesh-compile-contract
```

- API documentation: [docs.rs/axiolid-mesh-compile-contract](https://docs.rs/axiolid-mesh-compile-contract)
- Reference page: [axiolid.github.io/kernel](https://axiolid.github.io/kernel/reference/crates/axiolid-mesh-compile-contract)
- Source and issues: [github.com/axiolid/kernel](https://github.com/axiolid/kernel)
