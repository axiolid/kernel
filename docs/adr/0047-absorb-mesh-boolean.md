# 0047 — Absorb the mesh boolean into the provider crate

- **Status:** Accepted
- **Date:** 2026-09-08
- **Deciders:** Friedrich, axiolid
- **Supersedes:** [0014](0014-adopt-boolmesh-mesh-boolean.md)

## Context

ADR 0014 adopted `boolmesh` as a normal cargo dependency and explicitly
rejected vendoring: "Triggers MPL-2.0 file-level copyleft on our
modifications, and forks us off upstream fixes." That reasoning has been
overtaken by measurement on two fronts.

### The licence half of the objection does not apply

Axiolid is MPL-2.0. `boolmesh` is MPL-2.0. MPL-2.0 copyleft is
*file-level*: modified files stay MPL-2.0 and must be published. Axiolid's
source is already MPL-2.0 and already public, so absorbing adds no
obligation the project has not already accepted. ADR 0014 counted a cost
that is zero for a project with this licence.

### The maintenance half is real, and now cuts both ways

"Forks us off upstream fixes" assumes upstream fixes arrive. Measured
against the current release (0.1.9, the newest on crates.io):

- A depth-2 Menger sponge panics inside `boolmesh` with an out-of-bounds
  index in its own half-edge construction. It is a panic, not an error
  return, so no caller can catch it. Recorded in the benchmarks repo with
  an ignored regression test.
- Repeated grid-aligned subtraction is the workload BIM produces, and it
  is exactly the workload that fails.

Waiting on a third party for a defect that breaks a target workload is
itself a cost, and it is the cost the project is currently paying.

### Performance is the trigger, not the reason

Sphere-sphere boolean scaling (`benchmarks/`, subdivision 1..8, measured
on 20 cores) puts axiolid at 1.45x-2.12x Manifold's wall clock, the ratio
flattening at the top of the ladder rather than diverging. A `perf record`
profile of the union at 81920 triangles per operand attributes:

```text
21.5%  find_collisions / MortonCollider  broad phase
20.6%  sort (all call sites)             face morton, tri_halfs
 7.8%  Hmesh::new                        half-edge assembly
 7.6%  kernel01/02/03/11/12              the intersection math
 5.2%  Manifold::new_impl / new
 0.4%  axiolid's own convert.rs
```

Over 99% of the time is inside `compute_boolean`, reachable only as one
opaque call. Axiolid's own glue is 0.4%. No optimisation is available to
us from outside the crate, including the clearest one: `compute_boolean`
builds a full `Manifold` for its *output* — Morton BVH, coplanar index,
and a whole-mesh manifold validation — every bit of which the provider
discards, because `from_manifold` reads only positions and indices and
axiolid re-validates orientation itself.

## Decision

We will absorb the `boolmesh` algorithm into
`crates/providers/mesh/boolmesh/` as internal modules, and drop the
crates.io dependency.

- The upstream files keep Saki Komikado's copyright headers and their
  MPL-2.0 notice. Provenance is recorded in the crate README and in this
  ADR, not merely in git history.
- `compose` (primitive generators: cube, sphere, torus, cone, cylinder,
  extrude, fractal) is **not** absorbed. Axiolid has its own primitives
  and profile extrusion; carrying a second set would create two answers to
  one question.
- The absorption lands **behaviour-identical first**. The existing
  conformance, differential, determinism, and corpus suites are the gate.
  Optimisation is separate work on top of a green tree.
- The `MeshBoolean` contract seam from ADR 0003 is unchanged, so no
  consumer moves. Four crates name the provider in their manifests
  (`facade/axiolid`, `execution/compile`, `algorithms/discrete/minkowski`,
  `algorithms/discrete/decompose`); all of them depend on the *adapter*,
  never on upstream, and none of them change.

## Alternatives considered

| Option | Why not |
| --- | --- |
| Keep depending on crates.io `boolmesh` | The status quo. Leaves a known panic on a target workload unfixable by us, and puts 99% of boolean runtime behind a call we cannot enter. |
| Fork `boolmesh` to a git dependency | Every cost of absorbing, plus a divergent version in the dependency graph, plus no clean way to apply axiolid's gate, docs, and predicate work to it. Strictly worse than absorbing. |
| Swap to `manifold-rust` | ADR 0014 recorded it as the primary alternative and it remains attractive, but it relocates the same problem: another third-party boolean we cannot fix. Worth re-testing as a comparison, not as an answer to ownership. |
| Write a boolean from scratch | ADR 0003 and 0014 both judged this the wrong use of the project's scarcest resource, and that judgement stands. Absorbing a working, tested implementation is not the same as rebuilding one. |

## Consequences

**Positive**

- The `hmesh.rs` panic and the coincident-plane robustness gap become
  fixable in-tree rather than reportable upstream.
- The measured hot paths (broad phase, sorting, redundant output
  construction) become reachable. The redundant output `Manifold` build is
  a contained first target.
- Axiolid's exact predicates become available to the boolean's
  orientation decisions, which is the standing direction of the kernel.
- One fewer external crate; `glam` was already a transitive dependency.

**Negative / costs**

- Axiolid now owns roughly 3,650 lines of geometry code it did not write.
  Bugs in it are ours, including ones that predate absorption.
- Upstream improvements no longer arrive by bumping a version. Tracking
  `komietty/boolmesh` becomes a deliberate, manual activity.
- The absorbed files are MPL-2.0 with a third-party copyright holder.
  Headers must survive refactoring; a careless rewrite that strips them is
  a licence problem, not a style problem.

**Follow-ups / risks to watch**

- Absorbed code does not yet meet axiolid's documentation conventions.
  Bringing it up to standard is follow-up work, not a precondition.
- The absorbed modules use `glam` 0.30 for `USizeVec3`; `axiolid-core`
  pins `glam` 0.29. Two versions coexist in the graph today via the
  crates.io dependency, so this is unchanged by absorption, but it is now
  visibly ours to resolve.
- Upstream is edition 2024 and the workspace is edition 2021. One
  let-chain in `collider.rs` needs rewriting; there are no other
  edition-2024 constructs.

## Relation to existing code

- `crates/providers/mesh/boolmesh/` — receives the absorbed modules
  alongside the existing `convert.rs`, `provider.rs`, `cellular.rs`.
  Layer, role, and domain are unchanged.
- `crates/providers/mesh/boolmesh/src/convert.rs` and `provider.rs` — the
  only two files that imported upstream, via three symbols (`Manifold`,
  `compute_boolean`, `OpType`).
- `docs/adr/0003-pure-rust-mesh-boolean.md` — the contract seam it
  established is what makes this a one-crate change.
- `docs/adr/0014-adopt-boolmesh-mesh-boolean.md` — superseded by this ADR.


## Addendum, 2026-09-08: the absorbed code is now ours to shape

The initial landing was deliberately a faithful port, so that the
differential test could prove absorption changed nothing. That property
has been demonstrated, so the code was then reworked to the workspace's
own standards. The differential test in
`tests/absorbed_differential.rs` still passes against upstream 0.1.9, so
every change below is behaviour-preserving by construction.

`crates/providers/mesh/boolmesh/src/lib.rs` carried nine
`#[allow(...)]` attributes scoped to `mod csg` when the port landed.
**All nine are gone**; the module builds warning-free under both the
default and `--all-features` configurations with no suppressions.

### Changes that were latent defects, not style

- Five `x.abs() as usize` casts became `unsigned_abs()`. The original
  panics in debug builds for `i32::MIN` (which has no positive
  counterpart) and wraps silently in release, so the behaviour depended
  on the build profile.
- `Manifold::translate`, `rotate`, and `scale` each rebuilt the mesh
  through `.unwrap()`. They had no callers, so removing them deleted
  three panic paths rather than merely three functions. Upstream removed
  the same methods after 0.1.9.
- `Vec4::default()` followed by four field writes became a single
  `Vec4::new(...)`, so the value is never observable half-initialised.

### Parameter bundles, where they name a real thing

Six functions exceeded the argument limit. Rather than raise the limit,
five small structs were introduced -- each one named an existing concept
that was previously passed as loose parallel slices:

| Struct | Replaces | In |
| --- | --- | --- |
| `ResultEdges` | `hs_r`, `rs_r`, `face_ptr_r` | the result mesh being filled |
| `SourceSide` | `i03`, `hs_p`, `vid_p2r`, `fid_p2r`, `fwd` | one operand plus its index maps |
| `Windings` | `i03`, `i30`, `i12`, `i21` | operation-adjusted winding numbers |
| `EdgePoints` | `pt_old`, `pt_new` | where new vertices accumulate |
| `ShadowOperands` | `ps_p`, `ps_q`, `hs_q`, `ns` | the two operands of a shadow test |
| `SwapWalk` | `tag`, `visit`, `stack`, `edges` | edge-swap traversal state |

`Windings` also absorbed the four parallel `let` bindings that derived
those arrays, so the operation's sign convention now lives in one
constructor instead of four adjacent expressions.

`SourceSide` makes the symmetry of the algorithm visible:
`append_partial_edges` and `append_whole_edges` are each called twice,
once per operand, and the call sites now differ only in which side is
passed.

### Dead code removed rather than allowed

`compute_orthogonal`, `query_two_d_tree`, `Rect::overlap`, `dir_r`,
`Mat3`, and the `bounding_box` / `original_idx` fields were unused.
`original_idx` was always constructed empty, which invites a caller to
trust a value that is never populated.

Two of these were traps for a blanket fix. `Rect::overlap` looks dead in
isolation but serves `query_two_d_tree`, so removing only the leaf
breaks the build; they had to go together. And `tri_halfs_single` is
reported dead only under `--all-features`: it is the serial path,
replaced by `tri_halfs_multi` when `parallel` is on. It is now marked
`#[cfg(not(feature = "parallel"))]`, which states the fact instead of
suppressing the question -- deleting it would have broken the default
build.

### Verification

`cargo clippy --all-targets` and `--all-targets --all-features` both
report zero warnings with no `allow` attributes. The full suite passes
in both feature configurations, and `absorbed_differential` continues to
assert bit-identical vertex positions and volumes within 1e-12 against
upstream `boolmesh` 0.1.9 for union, intersection, and difference.

