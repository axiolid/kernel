# 0091 — Mesh sections report in-plane faces as regions and cut curves as chains

- **Status:** Accepted
- **Date:** 2026-10-28
- **Deciders:** Friedrich Schrödter
- **Supersedes:** —
- **Amends:** [0033](./0033-mesh-plane-section-contract.md)

## Context

Issue #359 (consumer: rheonda/cad). A floor plan is a section at, say, 1.0 m.
A railing, a slab or a wall top that lies exactly in that plane is a
face in the cutting plane. ADR 0033's `MeshPlaneSection::section` refuses such
a mesh as degenerate, because a two-dimensional overlap is not a curve, and it
refuses an open mesh as not manifold. The consumer needs both answered: what
is in the plane, and on which side its material lies, kept apart from the
curves where the surface crosses the plane. The consumer's decision, taken
there: a face in the plane comes back as a region tagged with a side, no cut
loop is invented for it, and the kernel never decides "just below".

`SectionOutcome`, `SectionEvidence`, `SectionLimits` and the conformance
types are public structs with public fields and no `#[non_exhaustive]`, so
adding to them is a breaking change.

## Decision

The existing method keeps its behaviour and messages. A new opt-in method,
`MeshPlaneSection::section_request`, takes a `SectionRequest` and returns a
`PlaneSection`. It is a provided method whose default refuses with
`Unsupported`, so every existing provider still compiles and a registry falls
back past one that does not implement it. All new types are
`#[non_exhaustive]`; every change is additive.

- **Regions.** Faces lying in the plane come back as `InPlaneRegion`s: an outer
  ring, hole rings, and an `InPlaneSide`. `Below` means the face normal agrees
  with the frame's `+z`, so the material lies on the `-z` side; `Above` is the
  opposite; `Both` is reserved for coincident faces of both facings.
- **Chains.** `SectionChain`s are transversal crossings only. A mesh edge in
  the plane yields a cut segment only when both of its triangles are off the
  plane and their third vertices lie on opposite sides. An edge next to an
  in-plane triangle is a region boundary, and an edge whose neighbours all lie
  on one side only touches the plane: nothing is emitted and the touching
  edge is counted in the evidence. A chain may therefore be open, ending on a
  region (`ChainEnd::InPlaneRegion`); its end points and the region vertices
  are projected by one function and are bit-identical.
- **Exactness.** Every decision is an exact `orient3d` or `orient2d` on the
  original binary64 coordinates with one axis dropped, never on projected or
  interpolated values. Regions come from topology, without a planar boolean:
  same-facing in-plane triangles sharing an edge form a component, boundary
  edges are linked into rings by fan order at each vertex, and the outer ring
  is told from a hole by the exact turn at its lexicographically smallest
  vertex. A ring that visits that vertex twice (a footprint or hole touching
  the outer boundary at a corner) is decided at the next lexicographic
  extreme; such a region has one weakly simple outer ring and no hole.
- **Exact chains.** A closed chain is made counter-clockwise from the exact
  vertex signs of its triangles: the cut leaves a triangle where the winding
  goes from above the plane to below it. A node is dropped only when its two
  cut segments are exactly collinear: a vertex node between two vertex nodes
  with a zero `orient2d`, or an edge node whose two triangles are exactly
  coplanar (`orient3d` of their four vertices). The nodes of in-plane
  triangles and the ends of an open chain are never dropped. No tolerance and
  no rounded coordinate takes part in a decision.
- **Refusals by name.** Overlapping in-plane faces of one side, in-plane faces
  of opposite sides that overlap, and a saddle vertex (more than two cut
  segments at one vertex) are `UnsupportedInput` with the constants
  `REFUSED_OVERLAPPING_IN_PLANE_FACES`, `REFUSED_IN_PLANE_FACES_ON_BOTH_SIDES`
  and `REFUSED_SADDLE_VERTEX`. A face that collapses under projection is
  `Degenerate`. Budgets get their own names: in-plane regions, region
  vertices and in-plane pairs.
- **Registry.** `MeshPlaneSectionRegistry::section_request` validates the
  source as before, then checks the answer against the documented output
  invariants: frame, cardinalities, finiteness, adjacent duplicates, canonical
  ring starts and orders, budgets and evidence counts. It does not check ring
  orientation, since a floating-point area sign is unreliable far from the
  origin (#274); the conformance suite checks orientation on fixtures.
  `register_request_conformant` registers a provider after both suites and
  rejects a claimed failure; a provider that does not claim the method passes
  on the base suite.
- **Conformance.** `RegionConformanceSuite` runs hand-derived fixtures: a cube
  through its top and bottom, a square frame with a hole, a parapet with a
  post rising out of it (the railing case), edge and point contact, and the
  two overlap cases, where a provider may refuse by name or merge. A provider
  whose first case returns `Unsupported` does not claim the method and passes.

The first landing (this record's code) covers closed solids. A second landing
adds a `SectionRole { Solid, Sheet }` on the request, so an open mesh (a
terrain, a shell) is sectioned too: chains may then end on the mesh boundary
(`ChainEnd::MeshBoundary`), boundary edges are counted, and a sheet in the
plane is a region. Both additions are non-breaking.

## Alternatives considered

| Option | Why not |
| --- | --- |
| Change `SectionOutcome` to carry regions | Breaking: the type has public fields and no `#[non_exhaustive]`. |
| Add `axiolid-overlay` to `axiolid-reference` to merge overlapping faces and name `Both` | Contradicts the dependency rationale of ADR 0090 and ADR 0036 and needs a maintainer decision (also on overlay tolerance snapping, #222). The overlaps are refused by name instead. |
| A new capability id | The operation is the same; a provider without the method refuses with `Unsupported`. |
| Let the kernel decide "just below" for an in-plane face | The consumer owns that decision; the kernel reports the side the material lies on. |
| Emit a closed cut loop around an in-plane face | Invents a curve that is not a crossing, and the loop would hide the region. |

## Consequences

**Positive**

- A floor plan at a level that coincides with a slab or railing top gets the
  filled region and the crossing curves, instead of a refusal.
- Output is deterministic and every decision is exact.

**Negative / costs**

- Elements that are soups of overlapping bodies with flush tops are refused by
  name until faces can be merged.
- Open meshes are not covered by the first landing.
- The reference implementation is quadratic in the worst case in the overlap
  sweep and in ring simplification; it is an oracle, not a performance tier.

**Follow-ups / risks to watch**

- Landing 2: open meshes, `SectionRole`, `ChainEnd::MeshBoundary`.
- Landing 3 (maintainer decision): `Both` and the merge of overlapping
  same-side faces, through `axiolid-overlay` or a focused provider.
- Landing 4: pair cut segments at a saddle vertex by fan order.

## Relation to existing code

- `crates/contracts/operations/mesh-section/src/request.rs`: the request,
  result, evidence and refusal names; `fixtures.rs` and
  `region_conformance.rs`: the suite.
- `crates/execution/dispatch/src/section.rs` and `section_request.rs`: the
  registry method and the outcome checks.
- `crates/algorithms/reference/src/section/`: `request.rs` (orchestration),
  `regions.rs` (in-plane regions), `chains.rs` (chains); the legacy path is
  unchanged.
- `crates/facade/axiolid/src/application.rs`: `Application::section_mesh_request`.
- `scripts/probe_section_regions_mutants.py`: the mutation probe.
