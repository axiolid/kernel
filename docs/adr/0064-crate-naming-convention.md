# Crate names are derived from architecture metadata

Status: accepted

## Context

Two crate pairs in this workspace read as the same pattern and are not:

```
axiolid-mesh-compile-contract   ->  trait MeshCompiler: Backend
axiolid-mesh-compile            ->  an implementation?

axiolid-mesh-boolean-contract   ->  trait MeshBoolean: Backend
axiolid-mesh-boolean-boolmesh   ->  an implementation
```

The obvious reading is that the suffix is optional. It is not. The two
crates have different ROLES:

- `axiolid-mesh-boolean-boolmesh` is `role = "provider.mesh"`. It wraps one
  engine and implements the capability directly.
- `axiolid-mesh-compile` is `role = "execution.orchestration"`. It is
  `impl<B: MeshBoolean> MeshCompiler for ReferenceMeshCompiler<B>` -- it
  TAKES a boolean engine and drives it. It is not a provider at all.

So the pair that looks like `<capability>` + `<capability>-<engine>` is
really `<capability>-contract` paired with an orchestrator whose name
happens to collide with the bare capability.

## Why the bare name must stay reserved

Both capabilities have MORE THAN ONE implementation in-tree:

```
MeshBoolean   -> BoolmeshBoolean (providers/mesh/boolmesh)
              -> ScalarBoolean   (algorithms/reference)
MeshCompiler  -> ReferenceMeshCompiler (execution/compile)
              -> GpuCompiler<E>        (execution/gpu)
```

A crate named `axiolid-mesh-boolean` would claim to be THE mesh boolean
while being one of two. The engine suffix is what keeps that honest, so
dropping it is the wrong direction.

## Decision

A crate name is DERIVED from metadata the manifest already declares,
`role` and `domain`, rather than agreed by convention:

| Role | Name |
| --- | --- |
| `contract.operation` | `axiolid-<domain>-contract` |
| `provider.*` | `axiolid-<domain>-<engine>` |
| everything else | named for what it is; may not end in `-contract`, may not take a name a contract reserves |

where `<domain>` is the declared domain with `.` replaced by `-`. The rule
is checked by `cargo xtask architecture check`, so it fails in the gate
rather than in review.

Keeping the `-contract` suffix rather than moving it to implementations
is deliberate. It is a WEIGHT label: resolving
`axiolid-mesh-compile-contract` alone pulls 20 crates against 51 for the
implementation, and a consumer reading the name in a dependency list can
tell before resolving it that naming the capability is cheap.

## Grandfathered, not silently allowed

`axiolid-mesh-compile` violates the third rule: it takes the bare name
`axiolid-mesh-compile-contract` reserves. The honest name is
`axiolid-graph-compile` -- it compiles a GRAPH, which is what its own
`domain = "graph.compile"` already says.

It is not renamed. It is published on crates.io at 0.1.0 and 0.2.0, and
`openbimrs/ifc` builds against it. A rename is a breaking change for
every downstream consumer in exchange for a name that reads better, and
48 internal references would move in the same commit.

Instead it is one EXPLICIT entry in the checker's exception table, with
its reason and its exit condition recorded there. The exception list is
itself checked: a stale entry -- one naming a crate that no longer
violates anything -- fails the gate, so the exemption cannot outlive the
problem quietly.

## What was fixed, and what was not

Two crates carried a `domain` that disagreed with their own name, which
is metadata drift rather than a naming problem, and both were corrected
in place -- no crate was renamed, no consumer affected:

```
axiolid-pointcloud-reconstruction-contract
    domain: operation.pointcloud-reconstruction -> pointcloud.reconstruction
axiolid-tessellation-contract
    domain: operation.tessellate                -> tessellation
```

`axiolid-tessellation-contract` remains the one contract whose trait is
not `: Backend` (`Tessellator: Debug + Send + Sync`). That is a contract
SHAPE question, not a naming one, and is left for its own decision.

`axiolid-backend-gpu` implements `MeshCompiler` but is named by hardware.
It is `role = "execution.context"`, not a provider, so the rule does not
claim it. Naming execution contexts by hardware is coherent on its own
terms.

## The gate was decoration until the probe caught it

The first version of `scripts/probe_naming_gate.sh` mutated crate NAMES
and reported all five mutations caught. Every one was a false pass: a
renamed crate breaks `cargo metadata` resolution, so the gate exited
non-zero before the naming rule ran at all. The probe was measuring
cargo, not the rule.

The mutations now move `domain` and `role`, which cargo ignores entirely,
and each mutation asserts on the RULE'S OWN message rather than on exit
status. This is why the probe exists: a naming gate that can never fail
is worse than none, because it advertises a guarantee it does not hold.

## Consequences

- A new crate whose name disagrees with its declared role and domain
  fails `cargo xtask architecture check`.
- The rule reads from metadata that already existed; no new manifest
  field was introduced.
- ADR 0004 said enforcement should be "machine-checkable rather than a
  naming convention". This does not reverse that -- it makes the naming
  convention itself machine-checkable, from the same metadata.
- One grandfathered exception exists and is visible in the source. If
  `axiolid-graph-compile` is ever published, the entry is deleted and the
  gate enforces the rule with no exceptions.
