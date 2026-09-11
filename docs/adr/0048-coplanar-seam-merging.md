# 0048 — Coplanar seam merging across chained booleans

- **Status:** Proposed
- **Date:** 2026-09-11
- **Deciders:** Friedrich Schrödter
- **Supersedes:** —

## Context

The exactness suite scores `(A-B) u (A^B) = A`. It fails on topology
while volume is exact (axiolid/kernel#100):

```
lhs  vol=4.800000000  chi=4  comps=2  manifold=true
rhs  vol=4.800000000  chi=2  comps=1  manifold=true
```

Splitting a solid and re-uniting the pieces returns two closed shells
instead of one solid. Each shell is individually closed and
two-manifold, and the volume is exact, so every pre-existing check
passed.

### This is an upstream limitation, not absorption damage

The same case run through crates.io `boolmesh` 0.1.9 directly gives
identical output:

```
UPSTREAM  (A-B)u(A^B): tris=32  chi=4  comps=2
OURS      (A-B)u(A^B): tris=32  chi=4  comps=2
```

ADR 0047 absorbed the algorithm faithfully. `is_coplanar` is present
in our tree and semantically identical to upstream. Nothing was lost.

### Mechanism

Coplanar merging is keyed on `Tref { mid, pid }`:

- `mid` — which INPUT manifold a face came from (0 or 1)
- `pid` — coplanar group id WITHIN that input, from
  `compute_coplanar_idx`, which runs only in `Manifold::new_impl`

Two faces merge only when they agree on BOTH. In `(A-B) u (A^B)` the
seam faces of `A-B` and of `A^B` both originate from `A`, but they
arrive as operands of a SEPARATE boolean call, each having been
retriangulated independently. The `pid` numbering is local to the
call that produced it, so the two sides cannot be recognised as the
same plane.

Measured consequences on the benchmark operands:

```
area(A)       = 29.600000
sum of parts  = 34.669992
area(union)   = 30.353368   <- removed 4.3166
excess over A =  0.753368   <- should be 0
shared faces between A-B and A^B = 0
```

The union IS removing seam area and stops short. Zero exactly-shared
faces confirms the two sides are combinatorially distinct along a
geometrically identical surface.

Rotation-dependent: the axis-aligned equivalent reconstructs exactly
(`chi=2 comps=1 area=29.600000`). Rotated coordinates are not exactly
representable, so the independent retriangulations diverge.

## Decision

We will make coplanar identity GEOMETRIC rather than provenance-based,
so that faces lying in the same plane merge regardless of which
boolean call produced them.

Scope, in landing order, each independently gated:

1. **Fixture first.** Add the reconstruction identity as a kernel test
   asserting `chi` and component count, not only volume. It must fail
   before any fix lands, and mutation-testing must show it can fail.
2. **Recompute plane groups on boolean output.** Derive `pid` from the
   plane a face lies in (quantised normal plus signed offset), rather
   than inheriting an input-local group id.
3. **Key merging on the geometric plane id** in `is_coplanar`, keeping
   `mid` only where provenance genuinely matters (sharp-edge and
   normal-flip guards in `collapse.rs`).
4. **Re-gate the whole suite**, including the drift and exactness
   tables, to prove no regression in volume accuracy or determinism.

## Alternatives considered

| Option | Why not |
| --- | --- |
| Weld duplicate vertices after the boolean | Tried and measured. Seam vertices are bit-identical, so the weld joins the index graph but leaves interior faces in place: `chi=4 comps=2 nonmanifold=0` becomes `chi=3 comps=1 nonmanifold=1`. An odd Euler characteristic is impossible for a closed surface. It converts an honest two-shell result into a corrupt one-shell result that still passes a volume check. Strictly worse. |
| Loosen the vertex merge tolerance | The duplicated vertices are separated by exactly `0.000e0`. No epsilon changes anything. |
| Cancel coincident face pairs | Measured: `A-B` and `A^B` share zero exactly-coincident faces. There is no pair to cancel. |
| Fork or patch upstream boolmesh | ADR 0047 already absorbed this code in-repo; crates.io boolmesh is a dev-dependency only. Upstream has the same defect, so a fork buys nothing a local fix does not. |
| Accept and document | Volume stays exact and both shells are closed, so impact is currently low. Rejected as the default because chained booleans are the primary CAD workload and silent topology loss compounds. |

## Consequences

**Positive**

- Chained booleans stop silently shedding topology, which is the CGAL
  failure mode for inexact constructions over consecutive operations.
- Output triangle counts drop where seams currently survive.

**Negative / costs**

- Touches absorbed upstream CSG, so we carry the divergence.
- A geometric plane key needs a quantisation tolerance, which is a new
  tunable and a new way to merge faces that should NOT merge. This is
  the main risk and is why the fixture lands first.
- Determinism must be re-proven: the provider already documents
  nondeterministic vertex ordering from upstream hashing.

**Follow-ups / risks to watch**

- Over-merging near-coplanar faces would round off genuine shallow
  features. Guard with a fixture holding a small dihedral angle.
- `compute_coplanar_idx` currently runs only in `new_impl`; calling it
  on boolean output reintroduces cost that ADR 0047 removed. Measure
  before and after rather than assuming it is negligible.

## Relation to existing code

- `crates/providers/mesh/boolmesh/src/csg/common.rs` — `Tref { mid,
  fid, pid }`, the identity the merge is keyed on.
- `crates/providers/mesh/boolmesh/src/csg/simplification/collapse.rs`
  — `is_coplanar`, and the sharp-edge guards that legitimately need
  `mid`.
- `crates/providers/mesh/boolmesh/src/csg/triangulation.rs` —
  `update_reference`, which stamps `pid` from the input manifolds.
- `crates/providers/mesh/boolmesh/src/csg/manifold.rs` —
  `compute_coplanar_idx`, called only from `new_impl`.
- `crates/providers/mesh/boolmesh/src/csg.rs` — `compute_boolean`,
  which no longer rebuilds a `Manifold` (ADR 0047) and so never
  recomputes plane groups.
- Benchmarks: the `reconstruction` row of the exactness table, and
  `topology_report` in `drift.rs`, both of which must stay green.
