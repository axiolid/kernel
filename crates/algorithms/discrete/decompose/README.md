# axiolid-decompose

Convex decomposition of a closed triangle-mesh solid. `Strategy::Exact`
splits at reflex features until every part is convex and the union
reproduces the input; `Strategy::Approximate` stops once each part is
within a stated concavity bound, giving far fewer parts. The returned
`Decomposition` always says which it is, and the approximate path reports
the concavity it actually reached. It works on meshes only, not on exact
B-reps.

```bash
cargo add axiolid-decompose
```

- API documentation: [docs.rs/axiolid-decompose](https://docs.rs/axiolid-decompose)
- Source and issues: [github.com/axiolid/kernel](https://github.com/axiolid/kernel)
