# 0017 — Solid boolean semantics are defined before an implementation is chosen

- **Status:** Accepted (all sections landed)
- **Date:** 2026-08-26
- **Deciders:** Friedrich, Hermes
- **Supersedes:** — (extends [0003](./0003-pure-rust-mesh-boolean.md) and [0014](./0014-adopt-boolmesh-mesh-boolean.md))
- **Package ownership amended by:** [0035](./0035-nested-ownership-and-capability-contracts.md)

> Boolean semantics and conformance obligations remain active. Point-in-time
> paths below are historical; portable Boolean contracts now live in
> `axiolid-mesh-boolean-contract`, shared admissibility in
> `axiolid-mesh-contracts`, and provider policy in `axiolid-dispatch`.

## Context

CSG is deferred until its contract exists. The risk is not that we lack a
backend — we have one — but that the backend we already adopted has begun
defining what the operation *means*.

Measured against the tree at the time of writing, that leak has started:

| # | Leak | Evidence |
| --- | --- | --- |
| 1 | `BooleanOperator` has exactly three variants — `Union`, `Intersection`, `Difference` — which is exactly `boolmesh::OpType` (`Add`, `Intersect`, `Subtract`), a 1:1 map in `provider.rs:97`. `axiolid-overlay`, designed contract-first, has **four** (adds `Xor`). The 2D and 3D operation sets disagree, and the 3D one has the shape of its backend. | `axiolid-core/src/operation.rs:5`, `axiolid-overlay/src/lib.rs:16` |
| 2 | Preconditions are enforced in the **L3 adapter**, not the L2 contract: `to_manifold` in `axiolid-mesh-boolean-boolmesh/src/convert.rs:45` decides closedness, orientation, and zero-volume. A second provider brings a second interpretation of "valid input". | `convert.rs:45-79` |
| 3 | `boolean()` returns a bare `TriMesh`. Overlay returns `OverlayResult { polygons, evidence }`; field returns `FieldEvidence`. The operation most in need of evidence has none, because the backend returns none. | `axiolid-kernel/src/boolean.rs:28` |
| 4 | `GeomError::Cancelled` is declared and produced **nowhere**; `ExecutionOptions` carries no token or deadline. The cancellation contract is fictional. | `grep Cancelled` → only the `enum` definition |
| 5 | The only provider declares `ScratchRequirement::Unbounded`, so the memory budget is unenforceable for every real call. | `provider.rs:70` |
| 6 | `axiolid-reference` implements no `MeshBoolean`. The single most consequential operation has **no oracle**, in direct violation of ADR 0012's ordering rule. | `grep MeshBoolean crates/algorithms/reference/src/` → no match |
| 7 | All five `axiolid-mesh-boolean-boolmesh` test files bind to the concrete `BoolmeshBoolean`. They test *boolmesh*, not *the contract*. A second provider inherits zero obligations. | `tests/{winding,batch,registry,conservation,fixture_issue_2019}.rs` |

ADR 0014 measured `boolmesh` honestly and adopting it was right. This ADR is
not a reversal. It says the API must be corrected **around** that provider
rather than **by** it, before a second provider — especially a native/C++ one
whose error model, lifetimes, and degeneracy conventions are far more
opinionated — makes the current shape permanent.

## Decision

We will define the six contracts below, land them with executable tests, and
choose no further CSG implementation until they exist.

### 1. Public operation semantics

- The operation set is **`Union`, `Intersection`, `Difference`,
  `SymmetricDifference`** — aligned with `axiolid-overlay`. A provider that
  cannot do symmetric difference composes it or reports `Unsupported`; it does
  not shrink the vocabulary.
- Semantics are **regularized**: the result is the closure of the interior of
  the set-theoretic result. Lower-dimensional residue — dangling faces,
  isolated edges, zero-thickness sheets — is never output. This is stated by
  Axiolid because backends disagree about it.
- `Difference` is **ordered**: subject minus tool, matching overlay.
- An **empty result is a value, not an error** (a tool containing the subject).
- **Coincident-face policy is Axiolid's**, not the backend's. Coplanar
  overlapping faces with agreeing normals belong to the boundary for `Union`
  and are removed for `Difference`; opposing normals cancel. Whatever is
  chosen, it is specified here and reported in evidence — this is the single
  largest source of cross-kernel divergence.
- `subtract_many` must be **observationally equal** to the sequential fold. It
  is a scheduling hint, never a different answer.

### 2. Input / topology requirements

Preconditions move **into L2** and are validated before dispatch, so every
provider receives identical, already-admissible input:

```text
SolidValidation::Structural            index range, finite coords, no degenerate triangles
SolidValidation::Closed                every edge has exactly two incident faces
SolidValidation::Oriented              consistent winding, positive signed volume
SolidValidation::SelfIntersectionFree  pairwise non-degenerate (opt-in: O(n log n)+ cost)
```

The level is caller-chosen because the last one is expensive and not always
needed. Failing a level is an `InvalidInput`/`NotManifold`/`Degenerate` error
naming the level — never a silent repair. A provider **may not** widen or
narrow admissibility; if it needs more, that is a capability gap, not a
precondition.

### 3. Diagnostics

`boolean` returns an outcome, not a mesh:

```text
BooleanOutcome { mesh, evidence }

BooleanEvidence {
    subject_triangles, tool_triangles, output_triangles,
    intersection_curves,        // topology actually computed
    coplanar_faces_resolved,    // where the policy above fired
    vertices_merged_by_tolerance,
    degenerate_configurations,  // tolerated, with what rule
    provider, precision, validation_level,
    result_verified,            // did the contract check run
}
```

This mirrors `OverlayEvidence` and `FieldEvidence`. Consistency across the
three is the point: one mental model for "what did the kernel actually do".

### 4. Resource and cancellation contracts

- **Budget:** `ScratchRequirement` stays, but `Unbounded` is a declared
  deficiency, not a default to live with. A provider seeking a bounded budget
  must publish a real bound. Refusal-before-allocation is already correct.
- **Cancellation becomes real.** `ExecutionOptions` gains a cooperative
  `CancellationToken` (no async runtime). Providers poll at a defined
  granularity, and cancellation is **safe**: `GeomError::Cancelled` or a
  complete result, never a partial mesh. A conformance test proves the provider
  actually polls; today `Cancelled` is unreachable and therefore untrue.
- Long operations report progress through evidence counters, not callbacks.

### 5. Scalar correctness oracle

Per ADR 0012, `axiolid-reference` owns a reference boolean **before** any further
provider is adopted. It is judged on correctness, not speed: exact predicates
via the existing `Certified`/`Sign` ladder, no threading, no intrinsics, never
feature-gated off. Quadratic is acceptable for a reference.

Implementation-independent invariants it anchors:

```text
vol(A \ B) + vol(A ∩ B) == vol(A)
vol(A ∪ B) + vol(A ∩ B) == vol(A) + vol(B)
(A \ B) \ B == A \ B                       idempotent
A ∪ ∅ == A,  A ∩ A == A                    identity
A △ B == (A ∪ B) \ (A ∩ B)                 symmetric difference consistency
```

Plus an **independent** point-membership Monte-Carlo cross-check. ADR 0014 ran
exactly this by hand, once, and it caught what conservation alone cannot — a
wrong-but-self-consistent result. It becomes an executable gate.

### 6. Provider conformance tests

A harness in `axiolid-kernel`, generic over `impl MeshBoolean` and **exported**
so out-of-tree providers can run it — the pattern already proven by
`axiolid-backend-gpu/tests/out_of_tree_executor.rs`. It asserts: all four
operations, empty results, precondition-rejection parity, evidence presence and
plausibility, budget refusal, cancellation honoured, bit-identical repeated
runs, agreement with the scalar oracle, and no inside-out output.

**A provider that has not passed the conformance suite is not registrable.**

## Status of each leak

Sections 1-4 landed on 2026-08-26; the table records what closed and how.

| # | Leak | Status | Landed as |
| --- | --- | --- | --- |
| 1 | 3D operation set mirrors the backend | **Closed** | `BooleanOperator` is `#[non_exhaustive]` with four regularized operands matching `axiolid-overlay`. `SymmetricDifference` is composed via `symmetric_difference_via_composition` where a provider lacks it natively. |
| 2 | Preconditions owned by the L3 adapter | **Closed** | `axiolid-kernel::solid::SolidRequirements` (`Structural` / `Enclosing` / `Oriented`), validated by the registry **before** dispatch, so admissibility cannot vary by provider. |
| 3 | Boolean returns a bare mesh | **Closed** | `BooleanOutcome { mesh, evidence }` with `BooleanEvidence` counters, mirroring `OverlayEvidence` and `FieldEvidence`. |
| 4 | Cancellation is fictional | **Closed** | `CancellationToken` + `ExecutionOptions::with_cancellation`; providers declare `CancellationGranularity` honestly rather than claiming responsiveness they lack. |
| 5 | Provider declares `Unbounded` scratch | **Closed** | Measured with a counting allocator (`axiolid-mesh-boolean-boolmesh/src/bin/scratch_probe/`): linear, under 1 KiB/triangle asymptotically, 2,184 B/triangle worst at small N after a warmup call (#110). On a rayon pool each worker adds 5 to 7 KB however small the input (#226), so a per-triangle bound is not an upper bound there: declared `Affine { base_bytes: 65536, bytes_per_element: 1536, bytes_per_worker: 16384 }` (`parallel`), checked against peaks measured on pools of 1 to 64 workers under a forced worst-case schedule by `tests/scratch_bound.rs`. |
| 6 | No scalar boolean oracle | **Closed** | `axiolid_reference::ScalarBoolean`: exact `orient3d` classification and ray parity, `O(n·m)`, no shared code path with `boolmesh`. Validated against analytic volumes in `axiolid-reference/tests/oracle.rs`. |
| 7 | Provider tests bind the concrete type | **Closed** | `axiolid_contracts::conformance` is generic over `impl MeshBoolean` and exported. `MeshBooleanRegistry::register_conformant` makes passing it a precondition of registration. |

### Consequences of the landed work

- `ExecutionOptions` is **no longer `Copy`**, because it carries a shared
  cancellation handle. An implicitly copied cancellation token is a footgun.
  Accessors now borrow instead of consuming, which is the better shape anyway.
- A latent bug surfaced and was fixed: `boolmesh` reports the intersection of
  disjoint solids as an *error* (`empty pos matrix`). The contract says an
  empty result is a **value**, so the adapter now translates it. Without the
  four-operand contract forcing a disjoint `A ∩ B`, this would have stayed
  hidden until a user hit it.
- The retired guards behaved as designed: three fired the moment their gap
  closed. A fourth, `precondition_ownership_gap_is_still_open`, kept passing
  because it probed for the name `SolidValidation` while the landed type is
  `SolidRequirements` — a false negative, and precisely why guards are
  mutation-probed rather than trusted.

### Closing sections 5 and 6

**The oracle is a separate algorithm, not a second boolean.** It classifies
whole operands by exact containment rather than cutting along intersection
curves. That makes it genuinely independent of `boolmesh` — the point of a
differential reference — at the cost of only answering non-interpenetrating
cases. It returns `Unsupported` for the rest rather than guessing, because an
oracle that approximates certifies wrong answers as correct.

That limit is narrower than it sounds. Disjoint, nested, identical, and
face-contact arrangements already pin identity, annihilation, idempotence,
commutativity, and containment across all four operations. The provider is
checked against the oracle on 16+ operation/arrangement pairs, and separately
required to handle the interpenetrating case the oracle refuses — so the
oracle's gap is exactly where the provider must earn its place.

**Conformance is enforced at registration.** A suite that must be remembered
is a suite that will be forgotten, so `register_conformant` returns the failing
report instead of admitting the provider. `register` remains for tests and
deliberately partial providers.

**Skips are reported, never counted as passes.** A provider cannot reach
"conformant" by declining everything: `ConformanceReport::exercised()` reports
what was actually proven, and the suite's own tests assert a minimum.

### Two bugs found by doing this

1. The oracle initially tested the *open* triangle interior for ray crossings,
   so a hit landing exactly on the diagonal shared by two triangles of a quad
   was missed by both — and interpenetration went undetected. Fixed by testing
   the closed triangle. This is the same shared-edge degeneracy the field
   sampler hit; it is the characteristic failure of this geometry family.
2. The oracle panicked on an empty operand by indexing `positions[0]`. A
   reference implementation that crashes takes down the harness meant to be
   judging correctness, so it now refuses with `InvalidInput`.

Both were found by the conformance suite running against the oracle — the
suite earning its keep before it ever judged a production provider.

### Retiring the deferral guard

`csg_deferral.rs` and `probe_csg_deferral.py` are deleted. Every gap they
tracked is closed, and a guard asserting "this is still missing" would now
assert the opposite of reality. The one durable rule they carried — no native
CSG backend — was always enforced more broadly by
`layering.rs::geometry_crates_do_not_declare_native_cpp_bridges`, which covers
`bindgen`, `cmake`, `cxx`, `manifold3d`, and `opencascade`.

## Alternatives considered

| Option | Why not |
| --- | --- |
| Bind a native/C++ CSG kernel now | The exact failure mode this ADR prevents: its error model, object lifetimes, tolerance conventions, and degeneracy handling would become Axiolid's public semantics by default. ADR 0011 already keeps native backends out of tree. |
| Keep `MeshBoolean` as-is and add providers | Locks in seven measured leaks. Each new provider raises the cost of fixing them. |
| Define the contract but skip the scalar oracle | Then conformance has nothing to compare against and reduces to self-consistency, which a wrong-but-consistent kernel passes. |
| Write the spec as prose only | Prose does not fail a build. Every clause above is testable and must be tested. |
| Treat `boolmesh` output as the reference | Makes one adopted crate the definition of correctness — the leak, formalised. |

## Consequences

**Positive**

- The 2D (`axiolid-overlay`), sampled (`axiolid-field`), and 3D boolean
  contracts share one vocabulary: explicit tolerance, validated input,
  structured evidence, no silent repair.
- A future native or GPU provider is a *conformance candidate*, not an author
  of semantics.
- `GeomError::Cancelled` and the memory budget stop being decorative.

**Negative / costs**

- The scalar reference boolean is real work and is the schedule's long pole.
- `boolean()` changing to return `BooleanOutcome` is a breaking change to an
  L2 trait; `axiolid-mesh-boolean-boolmesh` and the registry move with it.
- Precondition validation in L2 costs a pass the adapter was doing anyway.

**Follow-ups / risks to watch**

- **Open decision for Friedrich:** ADR 0012 says the scalar path lands *before*
  any optimized one. It did not, for booleans. Either (a) write the oracle now
  and keep 0012 intact, or (b) record a scoped exception naming `boolmesh` as
  a pre-oracle adoption. Recommendation: (a) — the oracle need not be fast, and
  without it conformance is unfalsifiable.
- The coincident-face policy must be pinned by fixtures before it is claimed.
- `subtract_many` grouping is currently gated only by volume equality; under
  the new contract it needs evidence equality too.

## Relation to existing code

- `crates/contracts/common/src/boolean.rs` — trait, registry, dispatch; the
  surface this ADR redefines.
- `crates/contracts/common/src/execution.rs` — `ExecutionOptions`,
  `ScratchRequirement`; gains cancellation.
- `crates/contracts/common/src/error.rs` — `GeomError::Cancelled`, currently
  unreachable.
- `crates/contracts/common/src/certainty.rs` — `Certified`/`Sign`/
  `EscalationLadder`, the arithmetic the oracle builds on.
- `crates/providers/mesh/boolmesh/src/convert.rs` — precondition logic to be lifted
  into L2.
- `crates/algorithms/reference/` — owner of the reference boolean; today has none.
- `crates/algorithms/planar/overlay/src/lib.rs` — the 2D contract this aligns with.
