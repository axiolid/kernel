# axiolid-evaluate

The scalar evaluation oracle for parametric geometry (ADR 0012, ADR 0036):
native-domain evaluation of analytic and B-spline curves and surfaces,
derivatives, jets, adaptive flattening, elementary surface inversion,
arc-length evaluation of intrinsic (natural-equation) and elevated curves
and of arc-length chains, arc length and its inverse on any curve,
the section frames of banked curves, stations (a section frame at a
distance along a curve),
and `ReferenceCurveEvaluator`, the reference implementation of the
curve-evaluation contract (ADR 0063). It has no mesh, spatial, measure or
provider dependency, so a parametric consumer gets evaluation without the
`axiolid-reference` umbrella. It favours obvious correctness over speed: no
intrinsics, threading or feature gates.

```bash
cargo add axiolid-evaluate
```

- API documentation: [docs.rs/axiolid-evaluate](https://docs.rs/axiolid-evaluate)
- Reference page: [axiolid.github.io/kernel](https://axiolid.github.io/kernel/reference/crates/axiolid-evaluate)
- Source and issues: [github.com/axiolid/kernel](https://github.com/axiolid/kernel)
