# 0044 — Pointcloud representation and reconstruction capability

- **Status:** Accepted
- **Date:** 2026-09-03 (accepted 2026-09-05)
- **Deciders:** GeneralPawz, Hermes
- **Supersedes:** —

## Context

Axiolid currently expresses two independent discrete/sampled families: connected
triangle meshes (`axiolid-mesh`, `representations/discrete/mesh`) and gridded
layered fields (`axiolid-field`, `representations/sampled/field`). It has no way
to represent point-sampled geometry at all — laser scans, photogrammetry
captures, and other unconnected 3D point sets with optional per-point
attributes (normal, color, intensity).

A pointcloud is neither a mesh (no topology, no adjacency) nor a sampled field
(no grid, no layering axis). It is its own discrete representation family:
an unordered or ordered set of positions plus optional per-point channels. It
needs its own value type, its own spatial-query support (KNN/radius rather
than mesh BVH traversal), and — for consumers who need it — a reconstruction
path that produces a mesh or exact result from the point set.

ADR 0035 already fixes the ownership DAG (`foundation ← representations ←
contracts ← algorithms/providers ← execution ← facade`) and the rule that a
contract is only added with a real implementation or independent consumer.
This ADR extends that DAG with one new representation family and one new
operation family, rather than introducing a parallel structure.

Source-format parsing (LAS, LAZ, E57, PCD, COPC) is explicitly out of scope:
ADR 0035's format-neutrality rule that keeps STEP/IFC out of `crates/` applies
identically here. Axiolid owns the in-memory point-set value; ingestion is a
sibling concern.

## Decision

Add a pointcloud representation and an optional reconstruction capability,
following the existing role DAG exactly:

- **Representation** (`representations/discrete/pointcloud/`, package
  `axiolid-pointcloud`): a validated point-set value with optional per-point
  attribute channels (normal, color, intensity). No topology, no adjacency, no
  algorithms. Depends only on `axiolid-core`. Sibling to `discrete/mesh`, not a
  child of it — a pointcloud is not a degenerate mesh.
- **Query** (extend `algorithms/query/spatial/`, package `axiolid-spatial`):
  add KNN and radius queries over `axiolid-pointcloud`, callback-based and
  deterministic-order like the crate's existing BVH/mesh queries. Broad-phase
  candidates remain explicitly distinct from exact adjacency, matching the
  crate's current invariant.
- **Contract** (new `contracts/operations/pointcloud-reconstruction/`,
  package e.g. `axiolid-pointcloud-reconstruction-contract`): a portable
  request/result/evidence/conformance schema for point-set → mesh
  reconstruction, mirroring `axiolid-mesh-boolean-contract` and
  `axiolid-mesh-section-contract`. Depends on representation values only,
  never on a concrete provider.
- **Provider** (new `providers/pointcloud/<backend>/`): a concrete adopted
  implementation of the reconstruction contract, mirroring
  `axiolid-mesh-boolean-boolmesh` in `providers/mesh/boolmesh`.
- **Execution** (extend `execution/dispatch`, package `axiolid-dispatch`): a
  `dispatch-pointcloud-reconstruction` feature registering the provider with
  the same fallback/failure policy as `dispatch-mesh-boolean` /
  `dispatch-mesh-section`.
- **Facade** (extend `crates/facade/axiolid`): two additive features —
  `pointcloud` (representation-only, zero algorithm/provider dependency, same
  contract `model` and `field` already honor) and `pointcloud-reconstruction`
  (contract + dispatch). No change to default features.

No universal point/mesh/field result enum is introduced. Reconstruction
failure, unsupported density, and provider absence are typed non-success
outcomes, consistent with every other operation contract.

## Alternatives considered

| Option | Why not |
| --- | --- |
| Represent a pointcloud as a mesh with no faces | Overloads `axiolid-mesh`'s topology invariants for a representation that has none; every mesh consumer would need to special-case "no faces means not actually a mesh." |
| Fold pointcloud storage into `axiolid-field` | A field is a dense grid on a layering axis; a pointcloud is a scattered, ungridded set. Forcing one into the other's shape loses the scattered case and adds grid-only invariants no pointcloud consumer needs. |
| Ship reconstruction as a single crate (representation + algorithm + provider fused) | Violates ADR 0035's role separation — collapses trust/provider/licensing boundaries the DAG exists to keep visible and swappable. |
| Vendor a source-format parser (LAS/E57) directly into the pointcloud crate for convenience | Breaks format-neutrality (ADR 0035); would repeat the STEP/IFC mistake the DAG was built to avoid. |

## Consequences

**Positive**

- Fills a real representation gap without inventing a new architectural
  pattern — every new package slots into an existing DAG layer.
- Reconstruction gets the same typed-refusal, evidence-carrying seam as
  mesh-boolean/mesh-section, so provider swaps and external
  `openbim.geometry` capability claims stay stable.
- `axiolid-pointcloud` remains independently useful (e.g. a consumer that only
  wants to store/query scan points, with zero reconstruction/provider
  dependency), matching the existing `field`/`model` leaf-crate pattern.

**Negative / costs**

- Six new or extended packages to land, gate, and document
  (`axiolid-pointcloud`, `axiolid-spatial` extension, the reconstruction
  contract, the provider, the dispatch feature, the facade features).
- `cargo xtask architecture check`'s allowlist and generated docs
  (`crate-map.md`, `dependency-graph.md`) need new entries; skipping this
  blocks every dependent sub-issue at merge time.
- No implementation exists yet; per ADR 0035 the contract must not be merged
  ahead of a real provider or independent consumer — this ADR authorizes the
  work, it does not preempt the ADR-0035 sequencing rule.

**Follow-ups / risks to watch**

- Source-format ingestion (LAS/LAZ/E57/PCD/COPC) is intentionally excluded
  from this ADR's scope and must land as a sibling crate outside `crates/`,
  not as a "just this once" exception inside the pointcloud representation.
- Reconstruction algorithm choice (Poisson vs. ball-pivot vs. alpha-shape) is
  a provider decision, not fixed by this ADR or its contract.

## As built

Landed 2026-09-05 across #60–#65. The shape matches this ADR; two points are
worth recording because they were decided during implementation:

- **The provider is not an adopted dependency.** The ADR anticipated
  vendoring an external reconstruction library. Instead
  `axiolid-pointcloud-reconstruction-sdf` composes two capabilities the
  kernel already owns: `PointIndex` neighbour queries (#61) and
  `axiolid-levelset` extraction (#87). No third-party numerics are adopted,
  so there is no new licensing surface and nothing unauditable. A Poisson or
  ball-pivot provider can still be added later — that is what the contract is
  for — but the contract is verifiable *now*, satisfying ADR 0035's rule that
  a contract lands with a real implementation.
- **Refusal is not a fallback trigger in dispatch.** When a provider says the
  data cannot support the request, that is an answer about the *data*, not a
  provider failure. Falling through would search for a provider willing to
  guess, which is exactly the silent degradation the contract forbids. Only
  `Unsupported`/`Unavailable` fall through.

The ingestion boundary is unchanged and enforced by the architecture gate:
`axiolid-pointcloud` declares `allowed-internal-dependencies = ["axiolid-core"]`
and no crate under `crates/` names a source format.

## Relation to existing code

- `crates/representations/discrete/mesh/` — sibling pattern for the new
  `crates/representations/discrete/pointcloud/`.
- `crates/representations/sampled/field/` — the other existing discrete/sampled
  sibling this decision distinguishes pointcloud from.
- `crates/algorithms/query/spatial/` (`axiolid-spatial`) — extended with
  point-set KNN/radius queries.
- `crates/contracts/operations/{mesh-boolean,mesh-section}/` — structural
  template for the new pointcloud-reconstruction contract.
- `crates/providers/mesh/boolmesh/` — structural template for the new
  pointcloud provider.
- `crates/execution/dispatch/` (`axiolid-dispatch`) — gains the
  `dispatch-pointcloud-reconstruction` feature.
- `crates/facade/axiolid/` — gains `pointcloud` and `pointcloud-reconstruction`
  features.
- Tracking issue: axiolid/kernel#59, decomposed into #60–#65.
