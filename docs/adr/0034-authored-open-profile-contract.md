# ADR 0034: Authored open profiles are graph declarations, not areas

- **Status:** Accepted
- **Date:** 2026-09-01
- **Deciders:** Axiolid maintainers
- **Supersedes:** —

## Context

A source model can declare an exact bounded open two-dimensional path as a profile. Such a declaration carries neither an enclosed area nor a width. Closing it would fabricate an edge and an area; treating it as a centre line would fabricate width semantics.

Axiolid already has two relevant but intentionally different vocabularies:

- `axiolid-profile::Profile` describes sweep sections with area semantics;
- `axiolid-model` can preserve atomic and relational exact 2D curves, including trimmed, composite, offset, and parameter curves.

Putting a graph `NodeId` in `axiolid-profile` would reverse crate layering. Copying a relational curve into an owned atomic profile value would lose exact composition and trim intent. Reusing mesh-plane section contours would replace authored exact geometry with evaluated approximate linework and would cross the IFC/plan boundary forbidden by ADR 0033.

## Decision

We will represent an authored bounded open profile as `axiolid_model::OpenProfile`, a graph payload whose `path` references a prior exact 2D curve-family node.

- `GeometryNode::OpenProfile` is distinct from `GeometryNode::Profile`.
- It carries no area, width, closure edge, tessellation, or manufacturing semantics.
- Graph construction accepts only conservatively bounded-open 2D paths:
  - source-open polylines with at least two finite controls and distinct authored endpoints;
  - source-open B-splines with finite controls and knots, coherent degree/knot/multiplicity cardinalities, and optional strictly positive finite weights;
  - trimmed 2D curves whose start and end selectors are finite, 2D-compatible, present in the declared preferred representation, and not exactly equal in a shared representation;
  - non-empty composites, finite-distance 2D offsets without a 3D reference direction, and parameter curves whose referenced paths satisfy the same contract;
  - instances preserve the referenced curve family and require finite transforms.
- Curve-family and bounded-open validation visit each reachable graph node at most once per check, so shared append-only DAGs do not expand into exponential work.
- Infinite lines, atomic circles/ellipses, source-closed or exactly self-closing polylines, source-closed or structurally invalid B-splines, non-finite curve/offset/instance data, 3D offset reference directions, empty composites, missing or non-finite trim endpoints, 3D trim selectors, exactly equal parameter/point trim endpoints, 3D curves, and non-curves are rejected.
- Endpoint coincidence and relational closure that cannot be established from identical authored selectors remain evaluator validation responsibilities. The payload preserves the source declaration; it does not claim a numerical proof that arbitrary relational endpoints differ.
- Existing solid operations continue to accept only area `Profile` nodes. Passing an `OpenProfile` to extrusion, revolution, or any other area-profile edge is a graph type error.
- The scalar mesh compiler reports standalone open-profile compilation as missing `CurveEvaluation`; it never triangulates or closes the path.

## Alternatives considered

| Option | Why not |
| --- | --- |
| Add an open variant to `axiolid_profile::Profile` | Makes area-only generator inputs admit a no-area state and pushes runtime rejection into every consumer. |
| Store `Curve2` directly in `axiolid-profile` | Cannot preserve graph-level trimmed/composite/offset curve relations without copying or flattening. |
| Store `NodeId` in `axiolid-profile` | Reverses the L1 representation → L2 graph dependency direction. |
| Reuse `CenterLineProfile` | A centre line implies a width and area; an authored open profile does not. |
| Close the path | Fabricates an unstated edge and enclosed area. |
| Use mesh-plane section contours | Substitutes evaluated approximate linework for authored exact geometry and confuses representation with plan derivation. |

## Consequences

**Positive**

- Importers can preserve exact open-profile intent without format-specific types.
- Invalid use as an area profile is rejected at graph construction.
- Existing area-profile and generation APIs remain semver-compatible.
- Composite and trimmed source curves remain graph references rather than flattened approximations.

**Negative / costs**

- Consumers that need points must explicitly request curve evaluation.
- Exact endpoint openness for arbitrary relational curves cannot be proven structurally and requires an evaluator.
- Importers must lower the path curve first and then append the declaration node.

**Follow-ups / risks to watch**

- Add a result-bearing bounded evaluator if consumers need a certified open-path endpoint/closure report.
- Do not broaden solid-operation profile edges to accept `OpenProfile`.
- Do not infer plans, regions, offsets, or manufacturing layers from this declaration.

## Relation to existing code

- `crates/representations/modeling/graph/src/node.rs` owns `OpenProfile` and `GeometryNode::OpenProfile`.
- `crates/representations/modeling/graph/src/validation.rs` enforces graph-reference and conservative bounded/open semantics.
- `crates/representations/modeling/graph/src/value.rs` provides sealed builder insertion.
- `crates/execution/compile/src/compiler.rs` classifies evaluation without manufacturing a mesh.
- ADR 0009 defines the graph layering; ADR 0033 separates authored geometry from mesh-derived section contours.

## Amendment 2026-10-14: open profiles compile to polylines (#311)

A consumer that needs an open profile's points no longer has to evaluate
the curve itself. `axiolid-mesh-compile-contract` gains a separate
operation, `PolylineCompiler::compile_polyline`, that compiles a curve node
to a 3D polyline and reports a certified chord bound the way a mesh's
deviation is reported: the requested budget and an upper bound on the
distance from every point of the exact curve to the polyline. The bound is
never absent; a curve that cannot be bounded is refused by name.

- It is a second operation, not an output enum on `MeshCompiler`: changing
  `compile_mesh`'s return type would break every consumer to serve one node
  family. `compile_mesh` keeps refusing curve nodes, open profiles
  included, with `CurveEvaluation`; it still never triangulates or closes
  the path.
- The reference compiler compiles `Curve2` (in `z = 0`), `Curve3` and
  `OpenProfile` over atomic curves, and their placements, by bisecting each
  smooth piece until every span's certified chord bound fits the budget.
  Curve relations (trims, composites, offsets) are refused by name for
  now, so an open profile whose path is a relation is too.
- A polyline carries no area or width, as the declaration does not.

Relation to code: `crates/contracts/operations/compile/src/polyline.rs`
owns the contract, `crates/execution/compile/src/polyline.rs` the
reference implementation.
