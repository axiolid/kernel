# 0090 — Mesh hidden-line drawings decide visibility exactly in the view's affine frame

- **Status:** Accepted
- **Date:** 2026-10-28
- **Deciders:** Friedrich Schrödter
- **Supersedes:** —

## Context

Issue #138 (consumer: rheonda/cad). A drawing of a model needs its edges as
seen from one direction, split into the parts a viewer can see and the parts
a solid hides. The ledger row C16 compares this with OCCT's `HLRBRep`
family. That family decides visibility with tolerances and polygonal
approximations, and its answer depends on them: a hidden edge that touches a
face, or sits in its plane, flips with the tolerance.

The kernel already draws on three pieces: exact predicates with an interval
filter and a dyadic fallback (`axiolid-exact`), the mesh plane-section
contract of ADR 0033 as the template for an operation contract with a
conformance suite, and ADR 0035 and ADR 0036 for where a provider lives.
Axiolid has no exact B-rep hidden-line work, and a smooth silhouette needs
curvature the triangles do not carry.

## Decision

Add an operation contract `MeshHiddenLines` for orthographic hidden-line
drawings of triangle meshes, a dispatch registry for it, and later a focused
scalar provider that decides visibility exactly. Nothing is rounded before a
decision; only the endpoints of output segments are rounded.

### Crates

- `axiolid-mesh-hidden-lines-contract` (`contract.operation`,
  `mesh.hidden-lines`, ADR 0064) in `crates/contracts/operations/`: request,
  view, edge selection, limits, segments, evidence, `HiddenLineError`, the
  `MeshHiddenLines` trait, and the conformance suite. It depends on
  `axiolid-contracts`, `axiolid-core` and `axiolid-mesh`.
- `axiolid-contracts`: capability `MESH_HIDDEN_LINES`
  (`org.axiolid.geometry.mesh-hidden-lines.v1`) and `Operation::MeshHiddenLines`.
  Both are additive: `Operation` is non-exhaustive, and the C ABI maps a
  variant it does not know to `Unknown`.
- `axiolid-dispatch`, feature `mesh-hidden-lines`: `MeshHiddenLinesRegistry`.
- `axiolid`: features `mesh-hidden-lines` and `dispatch-mesh-hidden-lines`,
  following the point-cloud pattern. They are not part of the `discrete`,
  `application` or `full` bundles: whether the application boundary and the
  C ABI carry hidden lines is a maintainer decision for when a provider exists.
- A provider, `axiolid-hidden-lines` (`algorithm.discrete`,
  `discrete.hidden-lines`), arrives in a later stage. It is its own package
  rather than part of `axiolid-reference`, which would pull `axiolid-exact`
  and its big-integer dependency into the umbrella crate that ADR 0036 keeps
  small (amended by ADR 0092, #377: `axiolid-reference` now depends on
  `axiolid-exact`).

### Semantics

The view is a frame `(o, x, y, z)`. A point `p` has image
`((p - o)·x, (p - o)·y)` and depth `W = (p - o)·z`; the viewer looks along
`-z`, so a larger `W` is nearer. The image and the depth are exact
polynomials in the binary64 inputs. Near and far clip by depth, view-parallel
and closed (`far <= W <= near`); an optional crop region is a box in the image.

- **Occlusion.** A point of an edge is hidden exactly when a kept triangle
  point lies strictly in front of it at the same image point. Triangles are
  closed sets, so touching an image boundary covers. Equal depth never
  occludes, so touching and coplanar geometry stays visible.
- **Faces seen edge-on.** A triangle of zero image area occludes through the
  depth range it spans above the image point (its fibre). A tall wall hides
  the edge of a low neighbour behind it.
- **Edge kinds.** One incident triangle: boundary. Two: silhouette when their
  image orientations differ, where "edge-on" is a third orientation of its
  own, so the test is symmetric and independent of winding; otherwise crease
  when the cosine of the angle between the triangle normals is strictly below
  the requested threshold. Flipping every triangle gives identical output.
- **Decisions.** Every visibility question along an edge is affine in the edge
  parameter. A split is a rational root; roots are ordered by exact sign
  evaluation, an interval filter first and a dyadic fallback after. No point
  is constructed to decide, and no tolerance flips a decision. The tolerance
  of `ExecutionOptions` applies to the mesh audit and frame validation only.
- **Output endpoints.** Every endpoint is the correctly rounded image of the
  exact rational point: each image coordinate (`U_A + t·(U_B − U_A)`, with
  `t` the exact split) is evaluated exactly and rounded once to nearest. An
  endpoint is a function of the exact point alone, independent of which
  constraint produced the split, so equal points are bit-identical and
  reversing every triangle changes neither the segments nor any evidence
  count except `exact_escalations`, which depends on the order the vertices
  are met in and is not part of the drawing. A piece whose rounded endpoints
  coincide is dropped and counted as collapsed.
- **Point edges.** A selected edge whose image is a single point (a vertical
  edge in plan) is counted in `point_edges` only, before clipping: it is not
  drawn and is never counted as clipped.
- **Near-plane caps** of closed solids occlude. A solid is cut when a vertex
  has `W > near` and another `W <= near`. A point of an edge with `W < near`
  (farther from the viewer than the plane) is hidden by the cap when the
  point of the plane above it, `(U, V, near)`, lies in the solid: a ray from
  it towards the viewer crosses the solid's triangles beyond the plane
  (`W > near`) an odd number of times. The crossings are counted with the image point perturbed by
  `δ = (ε, ε²)`, so a point on an image edge or on the diagonal two
  triangles share lies in exactly one of them: for an identically zero image
  side with `(dx, dy) = q' - p'` the perturbed side is `-sign(dy)`, or
  `sign(dx)` when `dy` is zero, times the triangle's orientation. Triangles
  seen edge-on cover no area and do not count. Sheets are never capped. A
  triangle whose plane is exactly at the near plane over the edge's image is
  not a crossing and needs no refusal: over its image every point with
  `W < near` is already hidden by that triangle itself (its point has depth
  `near`, which is kept), and elsewhere it contributes no toggle either way.
  A beam, lintel or sill with a face at the cut height is therefore drawn.
- **Evidence** counts selected, point, clipped and collapsed edges, candidate
  pairs and exact escalations, so a consumer can tell how hard a drawing was.

### Refusals carry the mesh as data

Consumer rheonda/cad drops and names elements one at a time. A refusal whose
only trace of the offending mesh is a sentence cannot be handled that way.
The trait and the registry therefore return `HiddenLineError { mesh:
Option<MeshId>, error: GeomError }` instead of a bare `GeomError`:

- `mesh` is `Some` when one input mesh is at fault: a duplicate id (the second
  occurrence), a mesh without triangles, an invalid index or non-finite
  position, a degenerate or exactly collinear triangle, a non-manifold edge,
  inconsistent winding, a `Solid` that is not closed.
- `mesh` is `None` for refusals of the whole request: limits, frame, near and
  far, crease threshold, region, an empty mesh list, memory, cancellation and
  a missing provider.
- `Display` prefixes `mesh <id>: `; the inner message does not repeat it. The
  type implements `std::error::Error` with the `GeomError` as its source, and
  `From<GeomError>` for a request-level refusal.
- Fallback to the next provider happens only when the inner error is
  `Unsupported` or `Unavailable`.
- The conformance suite names the mesh each refusal case must blame
  (`WrongRefusalMesh`).

The refusals are: `InvalidInput`, `NotManifold`, `Degenerate`,
`BudgetExceeded` (named resources "hidden-line meshes", "hidden-line source
vertices", "hidden-line source triangles", "hidden-line candidate pairs",
"hidden-line output segments", "memory"), `Cancelled`, `Unsupported`. A mesh
whose index count is not a multiple of three, or that has no triangles, is
`InvalidInput` naming the mesh (the section registry refuses the same).
A near plane that cuts a `Solid` is drawn with the solid's cap; the one
refusal is the tie described above. `UnsupportedInput` does not fall back to
another provider, by intent: a provider that supports the operation but not
an input is not asked to hand it on, and the caller should see that. A
provider refusal naming a mesh that is not in the request is a contract violation of
the registry, not a refusal passed on. There
is no undecided visibility in the first version; `Visibility` is
non-exhaustive so one can be added.

### Scope of the first version

Orthographic only; view-parallel clipping and a crop box; open and closed
meshes; boundary, silhouette and crease edges, visible and hidden. Left out
and refused by name or absent from the contract: perspective, oblique clip
half-spaces with caps (the request is non-exhaustive, so they can be added),
non-manifold edges, inconsistent winding, degenerate triangles, a `Solid` that
is not closed, exact B-rep and smooth-silhouette hidden-line removal. Ledger
row C16 stays narrow.

## Alternatives considered

| Option | Why not |
| --- | --- |
| Tolerance-based hidden-line removal as in OCCT `HLRBRep_PolyAlgo` | Its answers flip with the tolerance exactly where drawings are most delicate (touching boxes, coplanar faces). |
| Round the view transform first, then decide in floating point | The rounding moves points across the very boundaries the decision is about. |
| Rasterise or sample midpoints | Sampling misses slivers and cannot say which pieces are maximal. |
| Provider inside `axiolid-reference` | Pulls the exact tier and its big-integer dependency into the umbrella crate (ADR 0036). |
| Near-plane caps through `MeshPlaneSection` | Needs coplanar-face handling the section contract does not have yet (#359). |
| Return a bare `GeomError` and put the mesh id in the text | A consumer would have to parse it; the id is data. |
| One `MeshId` field on a new `GeomError` variant | Changes a shared, non-contract-owned enum for one operation. |

## Consequences

**Positive**

- A drawing's visibility is reproducible and independent of tolerance and
  winding; touching and coplanar cases have a defined answer.
- Consumers handle a refusal per mesh without parsing text.
- The contract and registry land before the provider, with a conformance
  suite that any provider can be held to.

**Negative / costs**

- Exact decisions cost more than a tolerance test where the interval filter
  fails; evidence reports how often.
- The brute-force candidate search is quadratic in the worst case, bounded by
  `max_candidate_pairs`.
- `HiddenLineError` is a second error type next to `GeomError`; only this
  operation uses it.

**Follow-ups / risks to watch**

- A spatial index for candidates.
- Oblique clip half-spaces with caps; coplanar faces of capped solids (#359).
- Hidden lines in the application bundle and the C ABI.
- Exact B-rep and smooth-silhouette hidden-line removal.

## Relation to existing code

- `crates/contracts/operations/mesh-hidden-lines/`
- `crates/contracts/common/base/src/capability_id.rs`, `capability.rs`
- `crates/execution/dispatch/src/hidden_lines.rs`, `frame.rs`
- `crates/facade/axiolid/` (`mesh_hidden_lines`)
- `docs/architecture/capability-ledger.toml` row C16
- Builds on ADR 0033, 0035, 0036 and 0064.
