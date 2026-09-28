# axiolid-pointcloud-reconstruction-contract

The portable contract for reconstructing a surface from an
`axiolid-pointcloud`: request, result, evidence, typed refusal, and a
conformance suite. A reconstruction is an estimate, so the contract makes a
provider report interpolated surface and resolved sample spacing, and
refuse rather than return an empty or fabricated mesh. Providers such as
`axiolid-pointcloud-reconstruction-sdf` implement it. See ADR 0044.

```bash
cargo add axiolid-pointcloud-reconstruction-contract
```

- API documentation: [docs.rs/axiolid-pointcloud-reconstruction-contract](https://docs.rs/axiolid-pointcloud-reconstruction-contract)
- Reference page: [axiolid.github.io/kernel](https://axiolid.github.io/kernel/reference/crates/axiolid-pointcloud-reconstruction-contract)
- Source and issues: [github.com/axiolid/kernel](https://github.com/axiolid/kernel)
