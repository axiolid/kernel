# 0111 — See-through triangles in hidden lines and visible faces

- **Status:** Accepted
- **Date:** 2026-11-14
- **Deciders:** Friedrich Schrödter
- **Supersedes:** —
- **Amends:** [0090](./0090-mesh-hidden-line-drawings.md), [0094](./0094-visible-face-regions-share-the-hidden-line-view-and-an-exact-arrangement.md), [0106](./0106-planar-face-groups-are-derived-exactly-and-group-boundaries-are-drawn.md)

## Context

Issue #496. A drawing of a building shows what lies behind its glazing: the
edges behind a window are visible, and the face behind it is seen through
the pane. Under ADR 0090 and ADR 0094 every kept triangle occludes, so a
caller can only leave the glass out, which drops its own edges and its own
face, or keep it and lose everything behind it.

The constraints are those of ADR 0090 and ADR 0094: every decision is exact
and tolerance-free, the output is deterministic, and a provider written
before the change ignores any new field, so it must never be handed a
request it would draw wrongly. The maintainer decisions on the design are
recorded on the issue (comment 6867725203).

## Decision

A triangle carries an `Occlusion`: opaque, the default, or see-through. A
see-through triangle hides nothing. A hidden-line piece counts the
see-through triangles in front of it, and a visible-face region of a
see-through source lists the surfaces seen through it.

- **Input.** `Occlusion { Opaque, SeeThrough }` and
  `TriangleOcclusion { AllOpaque, PerTriangle(&[Occlusion]) }`, both
  `#[non_exhaustive]` and defaulting to opaque, in the style of
  `FaceGroups`. `HiddenLineMesh` gains the `occlusion` field and
  `with_occlusion`. A `PerTriangle` slice whose length differs from the
  triangle count is `InvalidInput` naming the mesh. A provider that meets a
  variant it does not know refuses with `InvalidInput` naming the mesh.
- **Claim.** Both traits gain a defaulted
  `supports_occlusion(Occlusion) -> bool`, `true` for `Opaque` only, the
  pattern of ADR 0107. The registries skip a provider that does not claim
  `SeeThrough` when the request holds a see-through triangle and answer
  `Unsupported` when none does. A separate `SeeThroughConformanceSuite`
  (`run`, `run_faces`) infers the claim from `supports_occlusion`: a provider
  that does not claim must refuse a see-through request, one that claims must
  draw the see-through fixtures. Registration verifies a claim with it, so a
  claim is never routed to unverified.
- **Neutrality.** A request with no see-through flag, `PerTriangle` with all
  flags opaque included, gives exactly the output and the evidence of a
  request without flags.
- **Hidden lines.** A point is hidden only by a kept opaque triangle point
  strictly in front of it (closed sets, as ADR 0090) or by a cap.
  `HiddenLineSegment::seen_through` counts, for a visible piece, the kept
  see-through triangles of nonzero image area whose plane is strictly in
  front and whose image contains the point shifted by `(ε, ε²)`, the
  perturbation the near-plane cap already uses. So two panes that share an
  edge above a piece count once, an edge under a pane's border is decided by
  the shifted point, a pane's own edges do not count the pane (equal depth),
  and a pane seen edge-on is never counted. A piece is maximal with one
  `(visibility, seen_through)`; a hidden piece has the count 0 and is not
  split. The evidence gains `see_through_triangles` and
  `seen_through_segments`.
- **Visible faces.** The kept surfaces covering an image point are grouped
  into depth layers by exactly equal depth, layer 0 frontmost. The seen stack
  is every layer up to and including the first that holds an occluding
  surface (an opaque triangle or a cap), or all layers if none does. Every
  source in layer 0 gets a region, as before. `VisibleRegion` gains
  `behind: Vec<BehindSurface>`, where `BehindSurface { mesh, source, layer }`:
  empty for an opaque source or a cap, and for a see-through source the seen
  stack without its own layer-0 entry, so a source tied with it appears at
  layer 0. Within a layer there is one entry per `(mesh, source)`, ordered by
  the mesh's position in the request and then by source; a source crossed at
  two depths is listed once per depth. A region is maximal with one
  `(mesh, source, behind)`, and the output order breaks a tie of the old key
  by `behind`. `tied_sides` also counts a side whose stack has a deeper layer
  with two or more sources.
- **Uniform groups.** For visible faces every triangle of one face group has
  the same occlusion, else `InvalidInput` naming the mesh. `FaceGroups::Planar`
  joins two triangles only when their occlusion is equal (amends ADR 0106).
  With mixed groups a consumer could not tell from `(mesh, source)` whether a
  region occludes; allowing them later is not a breaking change.
- **Budget.** `VisibleFaceLimits::max_behind_surfaces` bounds the entries
  summed over all regions. `VisibleFaceLimits::new` keeps its six arguments
  and sets it to 0, so a see-through request fails with a named
  `BudgetExceeded` until the caller sets the limit, as for every explicit
  limit. The exact depth comparisons made only for stacks count against
  `max_candidate_pairs`; a request without see-through triangles makes none.
- **Caps.** The near-plane cap of a `Solid` is decided as before, from all of
  its triangles whatever their flags, and it occludes: a cut pane shows its
  cut face and stays see-through elsewhere. A see-through cap can later be a
  per-mesh opt-in without breaking this.
- **Projection.** The hidden-line counting is written on the shared
  constraint machinery, so it serves orthographic and perspective views.
  Visible faces stay orthographic until #387 P4 and inherit the stacks then.
- **The C ABI** is a separate step with its own symbols (ADR 0040 forbids
  layout changes); until then a C caller cannot mark a triangle.

## Alternatives considered

| Option | Why not |
| --- | --- |
| A `&[bool]` per mesh | Cannot grow a third kind of occlusion without a break. |
| A flag per face group | `FaceGroups::Planar` ids are derived by the provider, so the caller cannot name them. |
| A new `Visibility::SeenThrough` variant | A consumer's `_` arm would drop those edges, and the count would still need a field. |
| One region per frontmost source with a separate list of stacked sub-areas | Two geometries to snap-round consistently; one stack per region is what a consumer composites. |
| Mixed face groups | A consumer could no longer tell from `(mesh, source)` whether a region occludes. |
| Closed sets for the count under a pane's border | Two panes sharing a diagonal above an edge would count it twice. |
| Refuse a cut solid with see-through triangles | Every 3D view whose section box cuts a window would lose the glazing for that mesh. |
| An unlimited default for `max_behind_surfaces` | Breaks the rule that every output is bounded by an explicit limit. |

## Consequences

**Positive**

- The change is additive: every touched struct is `#[non_exhaustive]`, the
  constructors keep their arguments and the claim is a defaulted method.
- A provider that predates the field is never given a see-through triangle,
  and a request no provider can honour is `Unsupported`, never an opaque
  drawing.

**Negative / costs**

- A visible-face request with see-through triangles prunes fewer edge pieces,
  so the arrangement grows; an all-glass scene approaches the unpruned
  provider.
- A hidden-line edge behind glass loses the early-out of full cover by one
  triangle.

**Follow-ups / risks to watch**

- The C ABI symbols, sequenced with #387 P5.
- Perspective visible faces (#387 P4) must carry the stacks.

## Relation to existing code

- `crates/contracts/operations/mesh-hidden-lines`: `contract.rs`
  (`Occlusion`, `TriangleOcclusion`, `HiddenLineSegment::seen_through`,
  `supports_occlusion`), `faces.rs` (`BehindSurface`, `VisibleRegion::behind`,
  `VisibleFaceLimits::max_behind_surfaces`), `conformance/see_through.rs`,
  `fixtures.rs`, `face_fixtures.rs`.
- `crates/algorithms/discrete/hidden-lines`: `occlusion.rs` (the count),
  `prepare.rs`, `provider.rs`, `groups.rs` (planar groups split by
  occlusion), `faces/label.rs` (the stacks), `faces/regions.rs`,
  `faces/provider.rs`.
- `crates/execution/dispatch`: `hidden_lines.rs`, `visible_faces.rs`
  (routing, outcome validation), `view_checks.rs` (the length check).
- `crates/facade/axiolid`: `application/drawing.rs` (the claim verdict).
- `scripts/probe_see_through_mutants.py`.
