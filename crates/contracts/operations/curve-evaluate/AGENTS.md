# axiolid-curve-evaluate-contract

Names curve evaluation as a capability so a consumer can request it
without depending on an engine. Shaped like the mesh contracts.

- `contract.rs` the `CurveEvaluator: Backend` trait: `point_at`,
  `tangent_at`, `frame_at`, `distance_convention`.
- `convention.rs` `DistanceConvention`: which distance a provider
  measures, or `Unsupported`.
- `conformance.rs` the suite every provider must pass.

## Pitfalls

- Distance is NOT the native curve parameter. Only `Line`, `Circle` and
  `Intrinsic` recover it in closed form; the rest report `Unsupported`
  rather than returning a parameter dressed as a distance.
- `frame_at` is reference-up, not Frenet. The Frenet normal flips at a
  vertical inflection and is undefined on a straight, which would
  silently invert a placement. See ADR 0063.
- `Elevated` reports `PlanDistance`, which is shorter than 3D arc length
  by the grade factor. Do not convert silently.
