# axiolid-contracts

Common, provider-neutral contracts shared by every operation: backend
identity and descriptors, cancellation, diagnostics and `GeomError`, output
bounds and execution options, operation plans, and the integration profiles
a downstream application checks against. It defines no operation schema of
its own (those are the sibling `axiolid-*-contract` packages) and does no
provider selection or fallback, which belong to `axiolid-dispatch`.

```bash
cargo add axiolid-contracts
```

- API documentation: [docs.rs/axiolid-contracts](https://docs.rs/axiolid-contracts)
- Reference page: [axiolid.github.io/kernel](https://axiolid.github.io/kernel/reference/crates/axiolid-contracts)
- Source and issues: [github.com/axiolid/kernel](https://github.com/axiolid/kernel)
