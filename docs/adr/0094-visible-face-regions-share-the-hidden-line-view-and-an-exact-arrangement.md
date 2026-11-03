# 0094 — Visible face regions of meshes share the hidden-line view and an exact arrangement

- **Status:** Accepted
- **Date:** 2026-10-30
- **Deciders:** Friedrich Schrödter
- **Supersedes:** —
- **Amends:** [0090](./0090-mesh-hidden-line-drawings.md), [0093](./0093-mesh-sections-merge-in-plane-faces-on-an-exact-arrangement.md)

## Context

Issue #383 (consumer: rheonda/cad). A drawing of an elevation or a plan fills
each visible face with a surface pattern and each cut solid with a cut
pattern. ADR 0090 returns visible and hidden edges only. A consumer that
wants fills has to rebuild which face is frontmost where, from rounded
coordinates, and cannot make its fill boundaries lie exactly on the lines the
kernel drew.

## Decision

A second operation trait, `MeshVisibleFaces`, with its own capability
`MESH_VISIBLE_FACES` (`org.axiolid.geometry.mesh-visible-faces.v1`) and
`Operation::MeshVisibleFaces`, lives in `axiolid-mesh-hidden-lines-contract`
beside `MeshHiddenLines`. It shares `OrthographicView`, `MeshId`, `MeshRole`,
`HiddenLineMesh`, `HiddenLineError` and `HiddenLineSource`, and has its own
request, limits, outcome, evidence and conformance suite.

- **Groups are caller-chosen ids per triangle.** `HiddenLineMesh` gains a
  `groups: FaceGroups` field (`Single`, the default, or `PerTriangle`).
  `MeshHiddenLines` ignores it. The kernel does not infer planar faces.
- **The result is a list of `VisibleRegion`s**: per mesh and `FaceSource`
  (`Group` or `Cap`), a closed region with an outer ring, holes and the
  frontmost triangles. Vertices are exact points rounded once to nearest, so
  they are bit-identical between regions and with hidden-line endpoints.
- **Equal depth over an area gives the area to every tied source**, as ADR
  0090 says equal depth never occludes. The evidence counts the tied sides.
  Regions of one source never overlap.
- **A face seen edge-on covers no area**: it makes no region, never occludes
  a region, and is counted in the evidence.
- **The near-plane cap is a source of its own** at depth `near`, reported
  apart from the solid's faces because a drawing fills it with the cut
  pattern. Sheets are never capped and the far plane has no cap.
- **Refusals** are the hidden-line refusals, a group slice whose length is
  not the triangle count, the request budgets (meshes, vertices, triangles,
  candidate pairs, regions, region vertices, memory), a region that collapses
  under rounding, and cancellation. No tolerance enters a decision.

A second registry, `MeshVisibleFacesRegistry`, in `axiolid-dispatch` checks
the outcome like the hidden-line registry does.

**The exact arrangement is a crate of its own.** The point store, the
segment arrangement and the ring tracing that ADR 0093 built inside
`axiolid-reference` move to `axiolid-exact-arrangement`
(`crates/algorithms/planar/exact-arrangement`, layer `algorithms`, role
`algorithm.planar`, domain `planar.exact-arrangement`). The section of
`axiolid-reference` and the visible faces of `axiolid-hidden-lines` both
depend on it, so a sign about points, crossings and containment has one
source of truth (ADR 0036 §4). What stays in `axiolid-reference` is what
knows the plane: the plane values, the interning of cut nodes, the coverage
of in-plane triangles and the checks of cut rings. The section's output is
unchanged bit for bit, which every existing section test pins.

- *Closure.* The crate adds `axiolid-core`, `axiolid-contracts`,
  `axiolid-exact`, `axiolid-predicates` and `num-bigint` (already in the
  closure through `axiolid-exact`) and no third-party crate. It joins
  `expected_internal` of the `rust-facade-application`, `c-abi-profile` and
  `full` profiles, the only ones that contain `axiolid-reference`, and
  `forbidden_internal` of every profile that forbids `axiolid-reference`.
  `axiolid-reference` and `axiolid-hidden-lines` list it in
  `allowed-internal-dependencies`. `round_ratio`, the single rounding of an
  exact ratio, moves into the crate too, since the hidden-line output and
  the new region vertices must round alike.
- *Orientation.* The store takes a frame sign: the section works in a frame
  that may mirror the projected axes and passes the sign, the visible faces
  pass the positive one.
- *Inconsistency errors.* The crate has no backend id, so a state that
  cannot occur (overlapping pieces leaving a vertex the same way, parallel
  segments asked to cross) is `GeomError::Degenerate`, no longer the
  section's `BackendContractViolation`. No input reaches it.

| Option | Why not |
| --- | --- |
| Inside `axiolid-exact` | The arithmetic crate would depend on `axiolid-contracts` and carry about 3,000 lines of topology every exact consumer pays for. |
| `axiolid-hidden-lines` depends on `axiolid-reference` | Pulls the umbrella (curve, evaluate, measure, surface, primitive) into a provider closure, against ADR 0090 and ADR 0036. |
| Copy the code | Two sources of truth for a sign. |
| `axiolid-overlay` or `axiolid-arrangement` | A tolerance-snapping overlay and an editable binary64 subdivision: neither decides on exact points. |

## Alternatives considered

| Option | Why not |
| --- | --- |
| An extra output of `MeshHiddenLines` | `HiddenLineOutcome` is not `#[non_exhaustive]`, so a field is breaking; every hidden-line provider would also have to refuse face requests. |
| A defaulted method on `MeshHiddenLines` | Every hidden-line provider would nominally implement faces, which needs the claimed/unclaimed machinery of ADR 0091. |
| A new contract crate | Needs a contract-to-contract edge to share the view and mesh types and a second release line, for nothing a second trait does not give. |
| One group per mesh | Splits closed solids into one mesh per material, which loses caps and shared edges. |
| Kernel-computed planar groups | A plane-equivalence pass the caller already has the answer to. |
| Refusing equal-depth ties | Drops duplicated and flush finishes. |

## Consequences

**Positive**

- Every change is additive (new items, a field on a `#[non_exhaustive]`
  struct, a new variant of a `#[non_exhaustive]` enum); patch releases and no
  `semver-exceptions` entry.
- A hidden-line provider without faces does not claim them.
- cad calls both operations with one view and gets fills that meet lines
  exactly.

**Negative / costs**

- A second registry and a second conformance suite to maintain.
- Cost. `scenario/visible-faces` (`tools/benchmark/benches/scenario.rs`, grids
  of boxes, three runs on a noisy machine, run-to-run spread up to 50%)
  times `ScalarVisibleFaces` (`bvh`) against a full pair scan (`scan`) and
  against the BVH without the depth prune of hidden edges (`unpruned`). All
  three return identical regions, which the bench asserts. The machine was
  loaded by other jobs, so read the ratios as ranges over the three runs.
  The BVH beats the scan on the cut plan by 1.4-2.1x at 400 boxes and by
  4.1-6.6x at 2,500 boxes (50x50: 1.4-2.0 s against 7.0-9.6 s). The prune
  saves 2.7-3.2x in the isometric view at 400 boxes (`bvh` 2.0-2.8 s,
  `unpruned` 6.3-7.8 s) and 1.4-2.7x at 100 boxes; on the 50x50 cut plan
  `unpruned` was slower than `bvh` in every paired run (1.3-2.0x), while at
  400 boxes on the cut plan the runs disagree and no gain is claimed. In the
  isometric view and on the interpenetrating 6x6 scenes the BVH is not
  consistently faster than the scan; no gain is claimed there. A request
  takes seconds at a few thousand triangles, because exact intersection
  curves and labelling dominate; the 50x50 isometric grid takes over ten
  seconds and is not in the bench.
- A known limit. A generic rotated scene can be refused as `Degenerate` when
  two exact ring vertices round to the same binary64 value, because the
  region would collapse under rounding. The refusal is typed and names the
  mesh; the remedy belongs to the snapping decision of #377 (`TODO(#410)`).

**Follow-ups / risks to watch**

- A `FaceGroups::Planar` variant and a `ties` option can be added later.

## Relation to existing code

- `crates/contracts/operations/mesh-hidden-lines`: `faces.rs`, the
  `groups` field of `HiddenLineMesh`, the conformance suite.
- `crates/contracts/common/base`: the capability id and `Operation`.
- `crates/execution/dispatch`: `visible_faces.rs`.
- `crates/algorithms/discrete/hidden-lines`: the shared preparation.
