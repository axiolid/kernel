# 0048 — Coplanar seam merging across chained booleans

- **Status:** Proposed (steps 2-4 withdrawn, see amendment)
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


## Amendment 2026-09-11 — step 3 as written will not work

Probed before implementing. Two measurements invalidate the planned
approach.

**`is_coplanar` does not decide which faces survive.** Its only
callers are in `simplification/collapse.rs`, where it gates
EDGE-COLLAPSE eligibility during decimation. The retained face set
comes from the winding-number classification in
`boolean03/kernel03.rs`, which never consults `pid`. Re-keying
`is_coplanar` therefore cannot merge two shells into one.

**Making it maximally permissive destroys the mesh.** Forcing
`is_coplanar` to always return `true` and running the reconstruction
fixture gives:

```
InvalidInput("subject: mesh has no triangles")
```

Decimation collapses the solid away entirely. So the guard is load
bearing in the opposite direction from what step 3 assumed: it exists
to RESTRICT collapsing, and loosening it is destructive, not
corrective.

**The premise itself survives.** The seam IS coplanar: 62 triangle
pairs between `A-B` and `A^B` share a plane to 1e-9 in both normal
and offset. So a geometric plane key is still the right IDEA; it just
cannot be applied at `is_coplanar`, which runs too late and governs
the wrong decision.

### Revised direction

The union never treats the two operands as touching. With 62 coplanar
pairs present and zero EXACTLY-shared faces, the two shells are
coincident over a region but combinatorially disjoint, and the
classification stage retains both copies of the interface instead of
recognising it as interior.

That places the fix in coplanar-region handling during
INTERSECTION/CLASSIFICATION -- `intersect12` and `winding03` in
`boolean03/` -- not in decimation. This is materially deeper than
steps 2-4 assumed: it is the part of the algorithm that decides what
a boolean MEANS, and getting it wrong changes volumes, not just
triangle counts.

Step 1 (the fixture) stands and is unaffected. Steps 2-4 are
withdrawn pending a design for coplanar-region classification.


## Trace 2026-09-11 — winding number to emitted triangle

Requested before any further implementation. The path, with the file
and line where each step happens:

1. `boolean03/kernel03.rs::winding03` — produces `w03[v]`, a count
   per VERTEX of how many times that vertex is enclosed by the other
   operand. Computed by a planar-grid collision query, not by any
   face-pair test.
2. `boolean45.rs::Windings::new` — applies the operation
   coefficients: `i03[v] = c1 + c3*w03[v]`. For union `c1=1, c3=-1`,
   so `i03 = 1 - w03`.
3. `boolean45.rs::size_output` — accumulates
   `side_p[f] += |i03[tail]|` over each halfedge of `f`.
4. `boolean45.rs` — `keep_fs`: a face is emitted iff `side[f] > 0`.

### Why the reconstruction case degenerates

Measured on the benchmark operands:

```
A-B verts = 18, of which 18 lie ON the A^B surface
A^B verts = 10, of which 10 lie ON the A-B surface
```

Every vertex of both operands is a boundary case. Not one vertex is
strictly inside or strictly outside the other solid, so `w03` is
never driven to the value that would zero `i03` and drop a face.
With `side[f] > 0` for every face of both shells, both complete
surfaces are emitted — the doubled interface observed as `chi=4
comps=2`.

### `.abs()` is NOT the bug

The absolute value in step 3 looks like it prevents opposite-facing
coincident faces from cancelling, and removing it is the obvious
candidate fix. It is not: `side_pq` is also consumed at
`boolean45.rs:174-182`, where the per-face values are halfedge COUNTS
fed through `inclusive_scan` to produce `ih_per_f`, the output buffer
offsets. A negative or cancelled entry there corrupts allocation
rather than dropping a face. `side` is a retention count, and `|.|`
is correct for that role.

### Where the fix has to go

The gap is upstream of the keep decision: coincident boundary faces
never generate the cancelling windings that step 4 would act on. A
correct union of two solids meeting along a shared surface has to
classify that surface as INTERIOR and drop both copies, which means
recognising coplanar overlap during intersection/classification.
Confirmed such overlap exists here: 62 triangle pairs between `A-B`
and `A^B` share a plane to 1e-9 in both normal and offset, while zero
pairs are exactly equal as triangles.

`intersect12` (`boolean03/kernel12.rs`) finds edge-face crossings. A
face lying exactly IN another face crosses nothing, so it generates
no intersection record and the classifier never learns the two
surfaces are the same. That is the actual missing capability.

### Consequence for scope

This is a coplanar-overlap classification feature, not a repair of an
existing decision. It changes what the boolean considers interior,
so it moves volumes and not merely triangle counts, and it must be
gated on the full exactness and drift tables plus determinism before
it can be trusted. Estimated blast radius: `kernel12.rs` and
`kernel03.rs`, the two files that define boolean semantics.

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

## Measurement 2026-09-17 — the defect blocks chained booleans

Re-probed before attempting the fix, and one measurement changes the
priority argument recorded above.

The "Accept and document" alternative was rejected as the default but
described the impact as "currently low", on the grounds that volume
stays exact and both shells are closed. That understated it. The
doubled interface leaves edges with four incident faces, and a
half-edge mesh admits at most two, so the reconstructed solid cannot be
used as an operand at all:

```
provider.boolean(rebuilt, cutter, Difference)
  -> NotManifold("subject: edge (1, 3) has 4 incident faces;
                  a half-edge mesh admits at most two")
```

The axis-aligned control chains successfully (`chi=2 comps=1`, then a
further difference to `chi=0 comps=1`), which establishes the refusal
is caused by this defect and not by chaining itself.

So the failure mode is not a cosmetic topology count that a volume
check happens to miss. A chained boolean — the primary CAD workload —
hard-refuses on any solid that has been through a
split-and-reunite cycle with inexactly-representable coordinates.
Both measurements are now committed as `#[ignore]`d tests in
`tests/reconstruction.rs` alongside their passing axis-aligned
controls, so the fix has a gate and the controls guard the diagnosis.

### Re-confirmed, independently

The earlier findings were re-measured rather than trusted:

```
A        chi=2  comps=1        (one solid)
A-B      chi=0  comps=1        (torus: B punches clean through)
A^B      chi=2  comps=1        (the plug)
rebuilt  chi=4  comps=2        (two shells)

exact shared triangles between A-B and A^B : 0
coplanar pairs, same-facing                : 32
coplanar pairs, OPPOSED-facing             : 30
```

Thirty opposed-facing coplanar pairs is the interior-interface
signature: the torus tunnel wall and the plug's side wall face each
other across the same surface. Zero exactly-shared triangles confirms
no face-pair cancellation can act on it.

Instrumenting the winding computation on the union call shows the
classification is not wrong, only blind:

```
PROBE  nvP=18 nvQ=10  w03={0: 18}  w30={1: 10}  x12=20  x21=14
PROBE  nfP=36 keptP=36  nfQ=16 keptQ=8
```

`w03 = 0` for every `A-B` vertex and `w30 = 1` for every `A^B` vertex
are both CORRECT answers: torus vertices are not inside the plug, and
plug vertices are inside the torus region. With union coefficients
`i03 = 1 - w03`, every `A-B` face scores 1 and all 36 are retained,
including the tunnel wall that should have become interior.

That is the confirmation of the diagnosis already recorded: a face
lying exactly IN another face generates no edge-face crossing, so
`intersect12` never records the surfaces as coincident and the
classifier has nothing to act on. The windings do not need correcting;
the coincidence needs detecting in the first place.

### Refinement: the overlap is far smaller than "two doubled shells"

One measurement narrows the target usefully. By VERTEX INDEX the
rebuilt mesh is a clean two-manifold — every edge has exactly two
incident faces — which is why it passes validation in isolation:

```
rebuilt   verts=20  tris=32  edge incidence by index      = {2: 48}
                             edge incidence by COORDINATE = {2: 46, 4: 1}
                             distinct positions=20  distinct coords=18
```

Only two vertices are duplicated at identical coordinates, and only ONE
edge is four-incident. The two shells are not two fully independent
copies of a surface; they touch along a seam that is combinatorially
split at a single edge. That is what makes the next boolean refuse:
`Manifold` construction keys edges after welding by position, so the
pair collides there and nowhere else.

This also re-confirms why welding is not the fix. Welding by coordinate
and dropping the degenerate triangles gives:

```
welded  verts=18  tris=32  chi=3  comps=1
```

`chi = 3` is odd, which is impossible for a closed orientable surface —
exactly the result recorded in the "Alternatives considered" table. The
weld joins the index graph while leaving the interior faces in place,
converting an honest two-shell answer into a corrupt one-shell answer.
Re-measured here rather than carried over, and it reproduces.
