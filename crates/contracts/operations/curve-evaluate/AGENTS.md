# axiolid-curve-evaluate-contract

Names curve evaluation as a capability so a consumer can request it
without depending on an engine. Shaped like the mesh contracts.

- `contract.rs` the `CurveEvaluator: Backend` trait: `point_at`,
  `tangent_at`, `frame_at`, `distance_convention`.
- `convention.rs` `DistanceConvention`: which distance a provider
  measures, or `Unsupported`.
- `measure.rs` `CurveMeasure`: whether the caller's number is a length or
  a native parameter. Maps 1:1 to `IfcCurveMeasureSelect`.
- `conformance.rs` the suite every provider must pass.

## Pitfalls

- Distance is NOT the native curve parameter, and the two axes are
  separate: `DistanceConvention` says what a distance measures,
  `CurveMeasure` says whether the value IS a distance. Only `Line`,
  `Circle` and `Intrinsic` recover distance in closed form; the rest report
  `Unsupported` for distance but still answer `CurveMeasure::Parameter`.
- `frame_at` is reference-up, not Frenet. The Frenet normal flips at a
  vertical inflection and is undefined on a straight, which would
  silently invert a placement. See ADR 0063.
- `Elevated` reports `PlanDistance`, which is shorter than 3D arc length
  by the grade factor. Do not convert silently.
