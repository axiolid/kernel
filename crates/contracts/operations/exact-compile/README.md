# axiolid-exact-compile-contract

The contract for compiling an `axiolid-model` geometry graph into an exact
B-rep. An implementation either preserves analytic supports and trims or
refuses; there is deliberately no variant that returns a mesh, so a caller
that asked for exactness is never handed an approximation. It is a contract
only; `axiolid-mesh-compile` provides an implementation.

```bash
cargo add axiolid-exact-compile-contract
```

- API documentation: [docs.rs/axiolid-exact-compile-contract](https://docs.rs/axiolid-exact-compile-contract)
- Source and issues: [github.com/axiolid/kernel](https://github.com/axiolid/kernel)
