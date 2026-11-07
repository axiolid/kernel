# 0106 — Planar face groups are derived exactly, and group boundaries are drawn

- **Status:** Accepted
- **Date:** 2026-11-07
- **Deciders:** Friedrich Schrödter
- **Supersedes:** —
- **Amends:** [0094](./0094-visible-face-regions-share-the-hidden-line-view-and-an-exact-arrangement.md), [0090](./0090-mesh-hidden-line-drawings.md)

## Context

Issues #384 and #385, both follow-ups to #383. ADR 0094 takes face groups from
the caller, one id per triangle, and says the kernel does not infer planes. A
caller with a tessellated B-rep or a CAD mesh has to compute the planar faces
itself, with its own tolerance, before it can ask for visible faces. ADR 0094
lists `FaceGroups::Planar` as a follow-up.

ADR 0090 draws boundary, silhouette and crease edges. An edge between two
triangles of different groups on one plane (a material change, or two planar
faces that meet without a crease) is drawn by none of them, while the visible
face regions of the same request have an outline along it. Region outlines and
drawn lines then disagree.

## Decision

- **`FaceGroups::Planar`** derives the groups from the mesh. Two triangles are
  joined across an edge that has exactly two triangles whose model-space
  normals `n1` and `n2` are parallel and point the same way, decided exactly as
  `n1 × n2 = 0` in all three components and `n1 · n2 > 0`: the interval
  filter first, dyadic arithmetic when it cannot decide, no tolerance. A
  group is an edge-connected component of the joined triangles. An edge with
  one, three or more triangles joins nothing.
- **A group's id is `FaceGroup(i)` for the smallest triangle index `i` in
  the component.** It is the same in every view and under reversing every
  triangle. Two coplanar components that touch only at a vertex, or not at
  all, are two groups.
- **The provider resolves the groups once per request** (`groups.rs` in
  `axiolid-hidden-lines`), polls cancellation every 1024 triangles and edge
  runs, and maps an allocation failure to the `memory` budget refusal. The
  result equals the `PerTriangle` request carrying the resolved ids, which the
  tests check against an independent reference.
- **`FaceGroups::of` panics for `Planar`**: the groups depend on the geometry,
  and providers resolve them first.
- **The registry checks structure only for `Planar`.** It does not recompute
  the groups. A region's group id must be a triangle of the mesh and no later
  than any triangle the region lists, which holds because the id is the
  smallest triangle of the group.
- **Conformance** gains the isometric cube with planar groups (each face is
  group `2k`), a flat sheet of four coplanar triangles that is one group, and a
  fold of two opposite-facing triangles over one edge that are two groups.
- **`EdgeKind::GroupBoundary`** (amends 0090) is an edge with two triangles of
  different resolved groups that is neither a silhouette nor a selected
  crease. The precedence is silhouette, then crease (only when
  `crease_cosine` is set and the edge is a crease), then group boundary.
- **`EdgeSelection::group_boundary`** selects it and defaults to `false`
  (`EdgeSelection::with_group_boundary`), so the output of every existing
  request is unchanged bit for bit and the groups of a hidden-line request
  stay ignored unless it is set. When it is set, the groups are resolved as
  for visible faces (`Single`, `PerTriangle`, `Planar`), and a `PerTriangle`
  slice of the wrong length is `InvalidInput` naming the mesh, from the
  provider and from the registry before it calls one.
- **Region outlines and drawn lines agree.** For `Planar` groups the drawn
  group boundaries are the edges along which two visible regions of one
  mesh meet inside a face of the mesh.

## Alternatives considered

| Option | Why not |
| --- | --- |
| A tolerance on the angle between normals | Every decision in these operations is exact; a tolerance would make the grouping depend on a number the caller must tune, and on the order of the triangles. A caller who wants a tolerance still passes `PerTriangle`. |
| Group by plane equation, not by edge-connectivity | Two disjoint coplanar patches would be one group, and a region could not be a face of the model. |
| The largest or a view-dependent triangle as the id | The smallest index is stable under reversing the mesh and under the view, and equals the first triangle of the face in the usual builders. |
| Draw group boundaries whenever groups are given | It changes the output of requests that carry groups for visible faces only, and `MeshHiddenLines` ignored them until now. |
| A silhouette between two groups is a group boundary | A silhouette is drawn or not by its own selection; two kinds on one edge would make the output depend on both. |
| Recompute `Planar` in the registry | It would run the exact pass twice per request; the provider's result is tested against a reference instead. |

## Consequences

**Positive**

- A mesh gets its planar faces without the caller choosing a tolerance.
- The change is additive: `FaceGroups` and `HiddenLineMesh` are non-exhaustive
  or built through constructors, and existing requests are unchanged.

**Negative / costs**

- A rotated box whose binary64 vertices are not exactly coplanar falls apart
  into more groups than it has faces. That is the exact answer; a caller who
  wants those merged passes its own groups.
- A third-party provider that matches on `FaceGroups` must handle `Planar`;
  the enum was already `#[non_exhaustive]`.

**Follow-ups / risks to watch**

- On a tessellated curved surface every facet is its own group (neighbouring
  facets are not parallel), so planar grouping does not help there; exact
  B-rep drawing (#371) does.

## Relation to existing code

- `crates/contracts/operations/mesh-hidden-lines`: `faces.rs`
  (`FaceGroups::Planar`), `contract.rs` (`with_planar_groups`),
  `face_fixtures.rs` and `conformance/faces.rs`.
- `crates/algorithms/discrete/hidden-lines`: `groups.rs`, `faces.rs`,
  `edges.rs` (`collect`), `provider.rs`.
- `crates/execution/dispatch`: `visible_faces.rs`, `hidden_lines.rs`,
  `view_checks.rs`.
