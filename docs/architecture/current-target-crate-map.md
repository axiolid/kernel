# Axiolid restructure: implemented crate map

Implemented from base `372b16f64fa9be962df48751a654f8cda3f3b4a0` under ADR 0035.

## Physical ownership

| Package | Path | Role |
| --- | --- | --- |
| `axiolid-core` | `crates/foundation/core` | dependency root |
| `axiolid-curve`, `axiolid-surface`, `axiolid-primitive` | `crates/representations/analytic/*` | analytic values |
| `axiolid-profile` | `crates/representations/region/profile` | bounded region values |
| `axiolid-topology`, `axiolid-brep` | `crates/representations/{topology,brep}` | topology/exact B-rep values |
| `axiolid-mesh` | `crates/representations/discrete/mesh` | discrete mesh values |
| `axiolid-field` | `crates/representations/sampled/field` | sampled-field values/configuration/evidence |
| `axiolid-model` | `crates/representations/modeling/graph` | authored immutable graph |
| `axiolid-guarantees` | `crates/contracts/guarantees` | certification/escalation/precision vocabulary |
| `axiolid-contracts` | `crates/contracts/common/base` | common backend/execution/diagnostic contracts |
| `axiolid-mesh-contracts` | `crates/contracts/common/mesh` | shared mesh admissibility |
| operation contract packages | `crates/contracts/operations/*` | tessellation, mesh Boolean, mesh section, graph-to-mesh compile schemas |
| focused algorithm packages | `crates/algorithms/*` | reference, NURBS, construction, query, planar, sampled, repair |
| `axiolid-mesh-boolean-boolmesh` | `crates/providers/mesh/boolmesh` | concrete optional provider |
| `axiolid-dispatch` | `crates/execution/dispatch` | registration/fallback/device/budget policy |
| `axiolid-mesh-compile` | `crates/execution/compile` | reference graph-to-mesh execution |
| CPU/GPU packages | `crates/execution/{cpu,gpu}` | execution contexts/adapters |
| `axiolid` | `crates/facade/axiolid` | additive public feature facade |

The generated [crate map](./crate-map.md) and [dependency graph](./dependency-graph.md) are authoritative and freshness-checked by `cargo xtask architecture check`.

## Public migration table

| Before | After | Reason |
| --- | --- | --- |
| `axiolid-scalar` / `axiolid_scalar` | `axiolid-reference` / `axiolid_reference` | reference algorithm role, not scalar storage |
| `axiolid-generate` / `axiolid_generate` | `axiolid-construct` / `axiolid_construct` | construction semantics instead of vague verb |
| `axiolid-boolmesh` | `axiolid-mesh-boolean-boolmesh` | operation and concrete provider are explicit |
| `axiolid-kernel` aggregate | guarantees/common/mesh/operation contract packages plus `axiolid-dispatch` | portable contracts no longer own runtime policy |
| `axiolid-tessellate` | `axiolid-tessellation-contract` | package already defines a portable seam, not an implementation |
| `axiolid-compile` | `axiolid-mesh-compile` | output result domain is explicit |
| `GeometryCompiler::compile` | `MeshCompiler::compile_mesh` | exact B-rep and discrete mesh results cannot be silently conflated |
| `ScalarCompiler` | `ReferenceMeshCompiler` | reference implementation and output domain are explicit |
| combined `axiolid-field` | `axiolid-field` values + `axiolid-field-ops` algorithms | representation-only consumers remain light |

Facade feature migration:

- `kernel` is replaced by `contracts`.
- `field` is value-only; add `field-ops` or `field-navigation` for algorithms.
- `pointcloud` is value-only: `axiolid` + `axiolid-core` + `axiolid-pointcloud`
  and nothing else. Add `pointcloud-queries` for KNN/radius search,
  `pointcloud-reconstruction` for the contract,
  `dispatch-pointcloud-reconstruction` to select a provider at runtime, or
  `pointcloud-provider` for the bundled reference implementation.
- `mesh-boolean`, `mesh-section`, and `graph-compile` expose portable contracts.
- provider selection is opt-in through `dispatch-mesh-boolean` or `dispatch-mesh-section`.

## Resolved conflicts

- Exact B-rep ownership from ADRs 0020/0024 is preserved. Mesh compilation was renamed rather than generalized falsely.
- The existing default facade remains `mesh + cpu`; default-feature changes are outside this restructure.
- Source schemas, wire formats, and vendor interpretation remain outside Axiolid geometry packages.
- The neutral packages remain evidence targets for external `openbim.geometry` capability claims; no Pkl runtime or schema types are carried into Rust contracts.
- Dev-only upward edges are permitted only when explicitly allowlisted for integration/conformance tests; production/build edges still obey the role DAG.

## Downstream

`openbimrs/ifc`'s `ifc-geometry` crate must pin the landed Axiolid commit, retain source-format lowering outside Axiolid, and use the format-neutral authored `OpenProfile` graph declaration. A facade-only model consumer must continue to work with:

```toml
axiolid = { default-features = false, features = ["model"] }
```

without compiler, provider, field-operation, source-format, or GPU dependencies.

## Measurement

Fresh-target `cargo check -p axiolid --no-default-features` measurements on the same host and Rust 1.88.0:

| Feature | Before packages | After packages | Before elapsed | After elapsed |
|---|---:|---:|---:|---:|
| `model` | 17 | 17 | 1.915 s | 1.915 s |
| `field` | 5 | 5 | 1.781 s | 1.815 s |
| `pointcloud` | 3 | 3 | 3.54 s | 3.54 s |
| `pointcloud-reconstruction` | 10 | 10 | 4.04 s | 4.04 s |

The package split changes ownership and compilation units, not external normal-dependency count. This single cold run is noise-equivalent and is **not** evidence of a performance improvement.
