# axiolid

The feature-gated entry point to Axiolid's format-neutral geometry stack.
It re-exports the representation, algorithm, contract and provider packages
behind Cargo features and adds a small application layer over them; it owns
no geometry semantics of its own. `default = []`, so a consumer names the
capabilities it needs (for example `mesh`, `brep`, `nurbs`,
`tessellation`) or a bundle (`standard`, `discrete`, `parametric`,
`advanced`, `full`) and compiles nothing else. Exact geometry stays exact:
nothing here converts to a mesh unless the caller asked for a mesh and
supplied a tolerance. Every package behind a feature can also be used
directly.

```bash
cargo add axiolid
```

- API documentation: [docs.rs/axiolid](https://docs.rs/axiolid)
- Reference page: [axiolid.github.io/kernel](https://axiolid.github.io/kernel/reference/crates/axiolid)
- Source and issues: [github.com/axiolid/kernel](https://github.com/axiolid/kernel)
