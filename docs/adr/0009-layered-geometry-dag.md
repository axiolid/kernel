# 0009 - Layered geometry DAG and operation providers

- **Status:** Accepted
- **Date:** 2026-08-19
- **Deciders:** Friedrich, Hermes
- **Supersedes:** 0002 and 0004 backend/capability decisions; 0008 in part
  (neutral profile and shape vocabulary)
- **Package ownership amended by:** [0035](./0035-nested-ownership-and-capability-contracts.md)

> The neutral-DAG, runtime-selection, and provider-separation decisions remain
> active. Package names and paths below are historical; ADR 0035 decomposed
> `axiolid-kernel` into portable contracts, guarantees, operation contracts,
> dispatch, and execution-owned packages.

## Context

IFC translation needs exact curves, surfaces, B-rep, CSG, sweeps, mapped
instances, and tessellations without turning the geometry kernel into an IFC
implementation. Small mesh consumers must not pay for NURBS, topology, Rayon,
or GPU dependencies. Optimized implementations must remain replaceable on
x86_64, AArch64, discrete GPUs, and integrated GPUs.

Prior implementations also expose recurring failure modes: monolithic format
and kernel types, serialization/rendering leakage, global epsilon policy,
per-face tessellation cracks, non-deterministic hash iteration, compile-time
native CPU assumptions, and backend feature unification.

## Decision

We will use a downward-only crate graph, an immutable append-only neutral DAG,
small operation traits, and physically separate execution/adaptor crates.

- `axiolid-core` owns scalar, transform, bounds, and tolerance values.
- representation crates own exact values only; `axiolid-model` composes them by
  typed `NodeId` in a backward-reference DAG.
- algorithm crates own narrow extension traits; `axiolid-kernel` owns backend
  identity/policy/errors plus opt-in operation contracts and executable,
  operation-specific registries. Implementing the trait is the capability proof.
- `axiolid-backend-cpu` is an execution context; `axiolid-backend-gpu` contains narrow
  operation adapters. Neither claims an algorithm it cannot execute.
- `axiolid` is an opt-in facade; direct leaf-crate use remains supported.
- format adapters resolve schema semantics and units into neutral values. They
  cannot depend on concrete backends.
- CPU optimization uses runtime feature detection. GPU APIs implement an open,
  operation-specific seam such as `GpuGraphExecutor`; API and precision facts
  remain outside the stable geometry representation.

## Alternatives considered

| Option | Why not |
| --- | --- |
| One geometry crate with backend features | Feature unification hides coupling and makes minimal builds difficult to prove. |
| IFC-local primitive/request vocabulary | Duplicates the neutral model and forces every future format adapter to translate again. |
| Recursive boxed shape tree | Makes cycles, sharing, mapped instances, stable diagnostics, and bounded traversal harder. |
| `target-cpu=native` release builds | Produces host-specific binaries and cannot serve portable x86 or AArch64 users. |
| CUDA, Metal, Vulkan, or WebGPU types in kernel traits | Couples the stable contract to one API and its dependency weight. |
| Preserve a vendor kernel API for migration | Vendor adapters and ports must fit the neutral contract, not distort it. |

## Consequences

**Positive**

- IFC, future STEP/CityGML adapters, and algorithms migrated from prior art
  share one exact representation vocabulary.
- `axiolid --no-default-features`, capability bundles, CPU/GPU adapters, and leaf
  crates provide measured compile/dependency choices.
- Backends can be implemented out of tree under Rust's orphan rules because the
  implementer owns its backend type.
- Graph construction rejects forward references, so cycles cannot reach a
  compiler accidentally.
- Provider identity is observable, while operation support follows directly
  from executable trait registration rather than contradictory booleans.

**Negative / costs**

- More crates and explicit conversions increase initial scaffolding work.
- The GPU crate is an adapter contract, not a bundled GPU algorithm suite.
- Many IFC declarations are represented or assigned but not evaluated yet;
  coverage status must not be presented as implementation status.

**Follow-ups / risks to watch**

- Add exact evaluation and tessellation in small vertical slices, retaining a
  portable scalar oracle and differential tests for every optimized path.
- Benchmark before setting dispatch thresholds or claiming SIMD/GPU wins.
- Keep `crates/PLAN.md` and the executable 163-declaration IFC
  manifest current as capabilities land.
- Reject orphan scaffold files and files that become multi-responsibility.

## Relation to existing code

- `crates/AGENTS.md` and per-crate `AGENTS.md` files
- `crates/PLAN.md`
- `crates/representations/modeling/graph/src/graph.rs`
- `crates/contracts/common/`
- `crates/facade/axiolid-backend-{cpu,gpu}/`
- `packages/ifc/ifc-geometry/references/ifc4-add2-tc1-geometry-declarations.tsv`
- `packages/ifc/ifc-geometry/tests/declaration_manifest.rs`

## Amendment 2026-09-28: `crates/AGENTS.md` and `crates/PLAN.md` retired

[ADR 0078](./0078-context-lives-beside-the-code.md) removed both files. The
layer DAG and exact dependency allowlists are enforced by
`cargo xtask architecture check`. The standing rules are in the root
`AGENTS.md`, and open work is tracked in GitHub issues.
