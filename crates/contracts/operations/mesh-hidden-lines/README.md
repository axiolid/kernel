# axiolid-mesh-hidden-lines-contract

The portable contract for drawing the edges of triangle meshes as seen in a
view, split into visible and hidden pieces: request, view, edge
selection, limits, segments, evidence, a refusal type that names the mesh at
fault, and a conformance suite. It computes nothing itself; providers
implement `MeshHiddenLines`, and `axiolid-dispatch` selects one. See ADR 0090.

It also holds the sibling `MeshVisibleFaces` contract: the regions of the
same view where each face group, or the near-plane cap of a cut solid, is
frontmost. See ADR 0094.

A view is orthographic by default. A perspective view (`ViewProjection`, ADR 0107)
is part of the contract, and a provider says whether it draws one with
`supports_projection`; the registries never give it a view it does not
support, and `PerspectiveConformanceSuite` checks the claim: a hidden-line
provider that claims it must draw the perspective fixtures, and a visible-face
claim is not yet verifiable.

```bash
cargo add axiolid-mesh-hidden-lines-contract
```

- API documentation: [docs.rs/axiolid-mesh-hidden-lines-contract](https://docs.rs/axiolid-mesh-hidden-lines-contract)
- Reference page: [axiolid.github.io/kernel](https://axiolid.github.io/kernel/reference/crates/axiolid-mesh-hidden-lines-contract)
- Source and issues: [github.com/axiolid/kernel](https://github.com/axiolid/kernel)
