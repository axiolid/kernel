# 0094 — Visible face regions of meshes share the hidden-line view and an exact arrangement

- **Status:** Proposed
- **Date:** 2026-10-30
- **Deciders:** Friedrich Schrödter
- **Supersedes:** —
- **Amends:** [0090](./0090-mesh-hidden-line-drawings.md)

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
the outcome like the hidden-line registry does. The provider and the exact
arrangement it needs arrive in later steps of the same issue; this ADR is
amended then.

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

**Follow-ups / risks to watch**

- A `FaceGroups::Planar` variant and a `ties` option can be added later.

## Relation to existing code

- `crates/contracts/operations/mesh-hidden-lines`: `faces.rs`, the
  `groups` field of `HiddenLineMesh`, the conformance suite.
- `crates/contracts/common/base`: the capability id and `Operation`.
- `crates/execution/dispatch`: `visible_faces.rs`.
- `crates/algorithms/discrete/hidden-lines`: the shared preparation.
