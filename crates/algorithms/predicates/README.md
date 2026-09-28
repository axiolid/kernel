# axiolid-predicates

Certified geometric predicates: `orient2d`, `orient3d`, `incircle` and
`insphere`, built on error-free transformations and expansion arithmetic,
with static filters for callers that can bound their coordinates. Every
public predicate is a filtered cascade that escalates to exact arithmetic
instead of comparing against an epsilon, and returns a `Certified` sign.
The crate is deliberately narrow: no curve, surface, mesh, B-rep, provider
or big-integer dependency. `axiolid-reference` re-exports it unchanged
(ADR 0036); for signs of constructed values, see `axiolid-exact`.

```bash
cargo add axiolid-predicates
```

- API documentation: [docs.rs/axiolid-predicates](https://docs.rs/axiolid-predicates)
- Source and issues: [github.com/axiolid/kernel](https://github.com/axiolid/kernel)
