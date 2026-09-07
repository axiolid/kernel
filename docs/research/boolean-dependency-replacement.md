# Replacing `boolmesh`: dependency inventory and effort

Measured 2026-09-07 against `boolmesh 0.1.9`, the sole external dependency
of the mesh boolean path. Written because the boolean profile showed 0.0%
of that path is Axiolid code, which makes the dependency the whole story.

## What we depend on

The provider crate `axiolid-mesh-boolean-boolmesh` has exactly one external
dependency, and it has exactly one of its own:

```
axiolid-mesh-boolean-boolmesh
  └── boolmesh 0.1.9
        └── glam 0.30.10        (already ours transitively, via axiolid-core)
```

No C++, no system libraries, no build scripts. `rayon` is an optional
upstream feature we deliberately leave off (see `docs/architecture/
threading.md`). This is about as clean as a third-party dependency gets.

| property | value |
| --- | --- |
| license | **MPL-2.0** (file-level copyleft) |
| source | 3,248 code lines, 28 files |
| approach | from-scratch Rust port of Elalish's Manifold |
| arithmetic | `f32` or `f64`, selected by feature |
| robustness | epsilon/expansion (`expand: Real`), **not exact predicates** |


## Where the 3,248 lines sit

Measured per module, with a judgement on whether a replacement needs it:

| module | lines | do we need it? |
| --- | --- | --- |
| `triangulation/` (ear clip, flat tree) | 1,132 | **yes** — the core gap |
| `manifold/` (half-edge, collider, bounds) | 800 | mostly **no** — we own equivalents |
| `simplification/` (dedup, collapse, swap) | 664 | **yes**, eventually |
| `boolean03/` (intersection kernels) | 477 | **yes** — the core gap |
| `boolean45/` (assembly, classification) | 410 | **yes** |
| `tests.rs` | 442 | n/a — we have our own suite |
| `compose/` (sphere, torus, cube generators) | 362 | **no** — test fixtures, we have `axiolid-fixtures` |
| `common.rs`, `lib.rs` | 232 | trivial |

Subtracting what we already own or do not need (`manifold/`, `compose/`,
`tests.rs`) leaves roughly **2,700 lines of algorithm we would have to
write**, of which ~2,000 is the genuinely hard part: intersection kernels,
classification, and retriangulation.

## What we already own

Measured in this workspace:

| capability | lines | status |
| --- | --- | --- |
| exact predicates (`orient3d` etc, certified) | 1,436 | **stronger than boolmesh's epsilon approach** |
| BVH / spatial queries | 988 | done |
| mesh adjacency (`EdgeAdjacency`) | 1,245 | done |
| `ScalarBoolean` reference | 457 | partial — see below |


## The actual gap

`ScalarBoolean` is exact and total for disjoint, nested, and identical
operands. Its own doc names what it refuses:

> It reports `GeomError::Unsupported` when operand surfaces properly
> intersect, because resolving that requires retriangulating along the
> intersection curve -- the hard part of a real boolean, and the part an
> oracle must not fake.

So the replacement is not "3,248 lines". It is **one capability**:
resolve properly-intersecting surfaces. Everything else either exists or
is not needed.

That capability decomposes into:

1. **Triangle-triangle intersection segments** — we have
   `reference/triangle_triangle.rs` (140 lines) already.
2. **Stitching segments into intersection polylines** — new.
3. **Retriangulating each cut face against its polylines** — new, and the
   part boolmesh spends 1,132 lines on (constrained ear clipping).
4. **Classifying the resulting patches in/out** — `ScalarBoolean` already
   does this by exact ray parity for the non-intersecting case.
5. **Cleanup: dedup, degenerate collapse** — boolmesh spends 664 lines.

## Effort

An honest range, not a point estimate. Step 3 is where robust boolean
implementations historically fail, and where the epsilon-vs-exact choice
gets decided.

| scope | lines | calendar |
| --- | --- | --- |
| exact-predicate retriangulation, correct on our corpus | ~1,500-2,500 | weeks |
| plus performance parity with boolmesh | +unknown | months |
| plus the long tail of degenerate real-world IFC | +unknown | open-ended |


## Reasons to do it that are not about speed

The profile says replacing `boolmesh` is not a *performance* lever: its
cost is spread flat across intersection kernels, with 11.3% in the
allocator and no single hotspot to attack. A rewrite that merely matches
the algorithm would land in the same place.

Two non-speed arguments are stronger:

1. **Exactness.** `boolmesh` classifies with an epsilon expansion band;
   we own certified exact predicates. An exact-predicate boolean would be
   a genuinely different accuracy tier, not a faster version of the same
   thing -- and per the kernel's multi-path principle, that is a reason to
   ADD a path, not replace one.
2. **MPL-2.0.** File-level copyleft. It does not infect our crates, but it
   does constrain what a downstream consumer can do with modified
   `boolmesh` files, and it is the only non-permissive component in the
   dependency tree.

## Recommendation

Do not replace it. Extend `ScalarBoolean` toward the intersecting case as
a SECOND provider, exact-predicate based, and let dispatch choose:

- `boolmesh` where speed matters and epsilon classification is acceptable
- the exact provider where the answer must be certifiable, or where
  `boolmesh` refuses

That matches how the kernel already treats providers -- refusal is a
first-class outcome, and dispatch picks per request rather than globally.
It also means the work is incremental and useful at every step, instead
of a months-long rewrite that only pays off at the end.

Measure before believing any of the effort numbers above: they are
estimates, and the only measured figures in this document are the line
counts, the dependency tree, the license, and the profile shares.
