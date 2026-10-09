# axiolid-curve-evaluate-contract

The curve-evaluation capability as a contract: point, tangent and oriented
frame at a place on a `Curve3`, plus a conformance suite every provider must
pass. A caller says whether its number is a distance or a native parameter
(`CurveMeasure`), and a provider says which distance it measures for each
curve (`DistanceConvention`) or that it cannot. Frames are reference-up, not
Frenet: `x` tangent, `y` up, `z` right. The `*_on` queries read a seam side
(`SeamSide`) by the station seam rule of ADR 0082; a provider that does not
implement them refuses `Incoming` by name. This crate evaluates nothing itself; `axiolid-evaluate` provides the
scalar implementation. See ADR 0063.

```bash
cargo add axiolid-curve-evaluate-contract
```

- API documentation: [docs.rs/axiolid-curve-evaluate-contract](https://docs.rs/axiolid-curve-evaluate-contract)
- Reference page: [axiolid.github.io/kernel](https://axiolid.github.io/kernel/reference/crates/axiolid-curve-evaluate-contract)
- Source and issues: [github.com/axiolid/kernel](https://github.com/axiolid/kernel)
