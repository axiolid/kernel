# 0107 — Perspective views share the hidden-line and visible-face traits

- **Status:** Accepted
- **Date:** 2026-11-07
- **Deciders:** Friedrich Schrödter
- **Supersedes:** —
- **Amends:** [0090](./0090-mesh-hidden-line-drawings.md), [0094](./0094-visible-face-regions-share-the-hidden-line-view-and-an-exact-arrangement.md)

## Context

Issue #387, a follow-up to #383 and #138. `MeshHiddenLines` and
`MeshVisibleFaces` take an `OrthographicView` only (renamed `View` below). A perspective drawing
(an interior seen from an eye point) needs projective exact points, which
divide by depth, and a different near-plane cap. The issue asks for a design:
which view type, and whether one trait serves both projections.

The constraints come from ADR 0090 and ADR 0094. Every visibility decision is
exact and uses no tolerance; nothing is rounded before it is decided, and the
endpoints of the output are the exact points rounded once. A provider written
before this decision ignores any new view field, so an unchanged provider
must never be handed a view it would draw wrongly.

## Decision

One trait per operation serves both projections. The view gains a
`projection` field, the registries route by it, and perspective is
implemented in later steps behind the same contract.

- **The view.** `ViewProjection { Orthographic, Perspective }`
  (`#[non_exhaustive]`, default `Orthographic`) is a field of the view, which
  is already `#[non_exhaustive]`, so adding it is additive. The type is
  renamed `View`, since a type that can be perspective should not be called
  orthographic; `OrthographicView` stays as a `#[deprecated]` type alias, so
  the rename is semver-compatible. `View::new`
  builds an orthographic view as before; `View::perspective(frame,
  near)` builds a perspective one and refuses a `near` that is not finite and
  below zero with `InvalidInput`. The registries validate the same rule for a
  view built field by field.
- **The projection.** A perspective view is the central projection from the
  eye at `frame.origin` onto the plane `W = -1`: the image of a point with
  view coordinates `(U, V, W)`, `W < 0`, is `(U / -W, V / -W)`. `near` is
  required and below zero, so every kept point lies strictly in front of the
  eye and the division is defined.
- **Claim by `supports_projection`.** Both traits get a defaulted
  `fn supports_projection(&self, ViewProjection) -> bool` that is `true` for
  `Orthographic` only. The registries skip a provider whose answer is `false`
  and answer `Unsupported` when no provider supports the view's projection,
  before the topology audit, as for any request no provider can take. This
  protects an out-of-tree provider that predates the field and would ignore
  it. The in-tree providers keep the default and also refuse a perspective
  view by name (`UnsupportedInput`, input `perspective view`) when called
  directly.
- **Conformance.** `ConformanceReport` and `FaceConformanceReport` are
  exhaustive, so they cannot gain a field. A separate
  `PerspectiveConformanceSuite` (`run` for hidden lines, `run_faces` for
  visible faces) follows ADR 0091's `RegionConformanceSuite`: the claim is
  inferred from `supports_projection`, not probed. A provider that does not
  claim perspective must refuse a valid perspective request with
  `Unsupported` or `UnsupportedInput`. A provider that claims it fails with
  `ClaimUnverifiable` until the perspective fixtures exist (P3, P4), so a
  claim is never accepted unverified. `register_perspective_conformant` on
  both registries runs the base suite and this one. `register_conformant`
  also runs it whenever the provider claims perspective, so a claim is never
  routed to unverified: today such a provider is refused, as a
  `ReturnedError` of the case `perspective claim` (the base reports are
  exhaustive and take no new variant).
- **The exact method (P2 to P4).** Inverse depth `1 / (-W)` takes the place of
  `W`; it is affine on every triangle and along every edge. Every predicate
  becomes a homogeneous determinant of view-space vectors from the eye,
  multiplied by positive weights, so it stays exact and affine in the 3D edge
  parameter. The near-plane cap uses ray parity along the eye ray instead of
  the vertical ray. Endpoints are exact rationals rounded once.
  `ExactPoints` is already homogeneous, so snap rounding does not change.
- **Steps.** P1 is this contract: the type, the constructor, the claim, the
  routing, the refusals and the suite. P2 makes the provider's `Vert`
  homogeneous (`[Dyadic; 4]`) as a pure refactor: the segment hash on every
  scenario case is unchanged and `scenario/hidden-lines` shows no wall-clock
  regression. P3 draws perspective hidden lines and extends the binary64
  oracle to cast rays from the eye. P4 does the same for visible faces. P5
  adds the C ABI (`AXIOLID_PROJECTION_PERSPECTIVE = 2`). The capability
  ledger keeps "perspective" among the refusals until P4.

## Amendment: perspective hidden lines (P3)

`ScalarHiddenLines` and `ScanHiddenLines` draw perspective views and claim
them through `supports_projection`. The visible-face providers still refuse.
The exact method of the Decision is made concrete as follows.

- **The vertex.** A perspective vertex keeps its exact view coordinates
  `[U, V, W]` and the weight `h = -W`, so `[U, V, W, h]` is the point on the
  plane `W = -1` the vertex projects to; the image is `U / h, V / h`, and the
  endpoint of a piece is `(a·q - b·p) / (a·h_q - b·h_p)` rounded once, the
  formula the orthographic path already uses with `h = 1`. `U, V, W` stay
  affine in the edge parameter, so every constraint of the sign machinery
  (`recipe`, `param`) is unchanged in kind. The enclosure of a vertex holds
  the image coordinates `U / h`, `V / h` (the whole line when `h` is not
  certainly positive) and `W`; a box over vertices that straddle the eye
  plane is therefore unbounded in the image and the broad phase keeps it.
- **The predicates.** All are evaluated on vectors from the eye, with no
  division, so a triangle behind the eye needs no case. The side of an image
  line is the determinant of the rows `(h, U, V)`. An eye ray through `P`
  meets a triangle exactly when its three side determinants have the
  triangle's sign, and then reaches its plane at `P / s`. "In front of `P`"
  is `s > 1`, which is the sign of `det(t1 - t0, t2 - t0, P - t0)` against
  the triangle's orientation, the orthographic plane determinant read in view
  space. "Not past the near plane" is `near·s - W_P >= 0`, the determinant
  with `near·P - W_P·t0` in place of `P - t0`. A crop bound `b` is
  `b·h - U` or `U - b·h`. A depth bound compares `W` as before.
- **Edge-on triangles.** A triangle whose plane contains the eye has image
  orientation zero. Over an edge in the same plane through the eye, the
  extent of a triangle edge `CD` is `(C × P)·(C × D) >= 0` and
  `(P × D)·(C × D) >= 0`; the edge is in front when
  `(P × (D - C))·(C × D) > |C × D|²` and not past the near plane when
  `near·(P × (D - C))·(C × D) - |C × D|²·W_P >= 0`.
- **The cap.** The point of the near plane on the eye ray through `P` is in
  the solid exactly when an odd number of the solid's triangles are crossed
  by the eye line on the eye's side of the plane. The provider counts the
  other side instead, the triangles the ray meets between that point and
  infinity away from the eye (strictly beyond the plane, `near·s - W_P > 0`),
  which needs no triangle behind the eye and uses the same perturbed
  image-side clips as the orthographic cap. The independent oracle counts
  the eye line's own side, including behind the eye.
- **Conformance.** The perspective suite draws five fixtures (a box in
  front, the same box cut by the near plane, an eye inside a box, an exact
  non-dyadic boundary, an edge along the eye ray) and requires the refusal
  of a near depth that is missing, zero, positive or not finite. A provider
  that claims hidden-line perspective and ignores the projection fails them.
  The binary64 oracle of the provider tests casts rays from the eye.
- **Broad phase.** The hierarchy splits along the longest extent, and a
  perspective image is in units of one over the depth, tiny beside `W`. Its
  depth coordinate is `W` times a positive constant that brings the depth
  extent to that of the image. A positive factor keeps every order, so the
  probe still encloses what the exact test admits, and the pairs and the
  drawing do not change.
- **Compatibility.** Orthographic output is bit for bit unchanged: a
  constraint carries a `perspective` flag and the orthographic arm of every
  evaluation is the code it was.

## Alternatives considered

| Option | Why not |
| --- | --- |
| A second pair of traits for perspective | Doubles the registries, the conformance suites, the bundles and the C surface for one extra field of the same request. |
| The caller pre-transforms the mesh to a projective space and uses the orthographic provider | Rounds the transformed vertices before any decision is made, which ADR 0090 rejects; the result would not be the drawing of the input. |
| A new `PerspectiveView` type and a second request | The request, outcome and evidence carry the view unchanged, and a second type duplicates all of them. |
| A required trait method for the claim | Breaks every existing provider; a default of `Orthographic` only keeps them compiling and correct. |
| Probe the claim by running a perspective request | A provider that ignores the field would look like it draws perspective. The claim is a declared fact, and the suite checks the declaration against behaviour. |
| A field on `ConformanceReport` | The report is exhaustive with public fields; adding one is a breaking change. |

## Consequences

**Positive**

- The change is additive: the view type (under its old name too) and the
  traits keep their shape, and an existing provider is routed past a
  perspective view by the default.
- A perspective request that nothing can draw is `Unsupported`, never an
  orthographic drawing of a perspective view.

**Negative / costs**

- Code that names `OrthographicView` gets a deprecation warning until it
  moves to `View`.
- Until P4 a registered visible-face provider can only decline perspective.

**Follow-ups / risks to watch**

- A provider that claims perspective must also decide the near-plane cap by
  ray parity along the eye ray; the hidden-line fixtures arrived with P3, the
  visible-face ones come with P4.
- The C ABI (P5) needs its own enum value and a view struct field, and must
  keep the 0.4 call signatures.

## Relation to existing code

- `crates/contracts/operations/mesh-hidden-lines`: `contract.rs`
  (`ViewProjection`, `View`, `View::perspective`,
  `MeshHiddenLines::supports_projection`), `faces.rs`
  (`MeshVisibleFaces::supports_projection`), `conformance/perspective.rs`.
- `crates/execution/dispatch`: `hidden_lines.rs`, `visible_faces.rs` (routing,
  `register_perspective_conformant`), `view_checks.rs` (`validate_view`).
- `crates/algorithms/discrete/hidden-lines`: `view.rs`
  (`require_orthographic` for the face providers, `Vert::perspective`),
  `recipe.rs` (the eye-ray reading of each constraint), `caps.rs`,
  `occlusion.rs`, `provider.rs`, `faces/provider.rs`.
- `scripts/probe_perspective_routing_mutants.py`,
  `scripts/probe_perspective_lines_mutants.py`.
