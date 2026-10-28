# `openbim.geometry` capability boundary

Status: authoritative Axiolid-side boundary for the application- and
vendor-neutral [`openbim.geometry` Pkl
package](https://github.com/openbimrs/pkl/tree/openbim.geometry%400.1.0/packages/openbim.geometry).
The released package URI is
`package://openbimrs.github.io/pkl/openbim.geometry@0.1.0`.

Axiolid implements geometry capabilities. The Pkl package supplies a neutral
catalog and schema for declaring which capabilities an implementation claims.
It does not define an Axiolid runtime adapter, execution plan, provider order,
or product policy.

> [!IMPORTANT]
> The `openbim.geometry` 0.1.0 Axiolid manifest pins evidence to Axiolid revision
> [`ed07679ee698d960380ec4913f62ffd8c7a653a5`](https://github.com/axiolid/kernel/tree/ed07679ee698d960380ec4913f62ffd8c7a653a5). It describes claims for that
> immutable revision; it is not evidence for newer Axiolid source. A successor
> manifest must pin a newly reviewed Axiolid revision before claiming current
> conformance.

## Ownership

| Owner | Owns | Must not own |
| --- | --- | --- |
| `openbim.geometry` | catalog capability IDs; manifest, claim, requirement, and conformance-report schemas | Axiolid provider selection, execution requests, or product policy |
| Axiolid | Rust representation and operation contracts; implementations; tests; bounded evidence and typed refusals | Pkl runtime/schema types or downstream requirement policy |
| Consumer or checker | requirement profiles and matching a manifest to those requirements | rewriting implementation claims or interpreting absence as success |

The package remains outside the Rust dependency graph. Axiolid source may be
evidence for a Pkl claim, but geometry packages do not import Pkl modules or
carry Pkl values through public APIs.

## Two distinct capability identities

`openbim.geometry` uses a closed `CapabilityId` catalog such as
`openbim.geometry:tessellation.surface`. Axiolid operation contracts use their own
versioned IDs such as `org.axiolid.geometry.tessellate.v1`.

These IDs are intentionally **not interchangeable**:

- the Pkl ID names a neutral capability that implementations can claim;
- the Rust ID versions one concrete Axiolid operation contract;
- one Pkl claim can require evidence from multiple Rust contracts or
  representation packages;
- the external Axiolid manifest owns the mapping and its evidence locators.

String equality, shared prefixes, or a provider name never establish that
mapping.

## Axiolid operation evidence anchors

| Rust contract ID | Rust contract | Requirements | Guarantee / evidence boundary |
| --- | --- | --- | --- |
| `org.axiolid.geometry.tessellate.v1` | `axiolid_tessellation_contract::Tessellator` | typed graph/root and validated `TessellationOptions` | triangle mesh under explicit request parameters; no exact B-rep claim |
| `org.axiolid.geometry.mesh-boolean.v1` | `axiolid_mesh_boolean_contract::MeshBoolean` | oriented closed mesh operands, operator, and execution options | `BooleanOutcome` and `BooleanEvidence`; unsupported cases remain explicit |
| `org.axiolid.geometry.mesh-section.v1` | `axiolid_mesh_section_contract::MeshPlaneSection` | oriented closed mesh, plane frame, and typed limits/options | closed plane-local contours and `SectionEvidence`; unsupported topology fails closed |
| `org.axiolid.geometry.mesh-hidden-lines.v1` | `axiolid_mesh_hidden_lines_contract::MeshHiddenLines` | meshes with roles, an orthographic view, edge selection, and typed limits/options | visible and hidden edge pieces and `HiddenLineEvidence`; the exact scalar provider refuses non-manifold input and undecidable input |
| `org.axiolid.geometry.graph-to-mesh.v1` | `axiolid_mesh_compile_contract::MeshCompiler` | geometry graph/root and execution options | explicitly mesh-valued output; never an implicit exact B-rep result |

This table identifies possible evidence anchors; it does not add claims to the
Pkl manifest. The released manifest is the source of truth for its own declared
catalog capabilities and scope.

## Claim-to-evidence mapping

| `openbim.geometry` field | Axiolid evidence rule |
| --- | --- |
| `implementation.id` and `implementation.version` | identify the released implementation; claim evidence must still pin an immutable source revision |
| capability map key | use the package's closed catalog ID; map it explicitly to relevant Rust contracts and representations |
| `CapabilityClaim.version` | version the claim itself; do not reuse a crate or wire-format version implicitly |
| `level` and `exactness` | state only what the cited implementation and tests demonstrate within the declared scope |
| `scope` | list dimensions, input/output representations, and supported cases narrowly enough to verify |
| `limitations` | preserve unsupported, approximate, bounded, and partial cases explicitly |
| `evidence` | point to stable source, test, documentation, or benchmark locations at the claimed revision |

Missing geometry, unsupported representations, unavailable evidence, and
provider absence are non-success outcomes. An external conformance evaluator
may report `satisfied`, `unsatisfied`, or `indeterminate`; it must not turn
missing evidence into a successful claim.

## Provider obligations

A provider-backed claim is valid only when all of the following are true:

1. the claimed neutral capability and scope map to implemented Rust contracts;
2. requirements are admitted before execution;
3. unsupported inputs or missing evidence produce an explicit typed refusal or
   indeterminate result;
4. evidence identifies immutable source and the relevant suite/vector version;
5. device, memory, scheduling, retry, fallback, and provider ordering remain
   execution metadata—not capability semantics.

A provider name, registry entry, package manifest, or benchmark alone is never
capability evidence.

## Update protocol

When Axiolid capability behavior changes:

1. land and independently review the Rust contract, implementation, tests, and
   documentation;
2. identify the exact Axiolid revision and evidence locations;
3. update `openbim.geometry/conformance/axiolid.pkl` in the Pkl repository;
4. run the Pkl package's tests and release the successor package immutably;
5. update consumers only after that released package is verifiable.

This keeps Axiolid format-neutral while giving `openbim.geometry` precise,
inspectable implementation evidence.
