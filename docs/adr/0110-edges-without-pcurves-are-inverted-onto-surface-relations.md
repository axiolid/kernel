# 0110 — Edges without pcurves are inverted onto surface relations

- **Status:** Accepted
- **Date:** 2026-11-14
- **Deciders:** Friedrich Schrödter
- **Supersedes:** —
- **Amends:** [0103](./0103-b-rep-faces-on-swept-and-ruled-surface-relations.md) (D6: an edge use without a pcurve on a relation face)

## Context

Issue #454. ADR 0103 meshed B-rep faces over `SurfaceRelation::Revolution`,
`LinearExtrusion` and `Ruled`, but only when every edge use carried a pcurve:
a use without one was refused before any meshing. IFC advanced faces on
swept surfaces and most Revit edges carry their 3D curves only, so these
faces could not be meshed.

Faces over atomic surfaces already handle such uses (ADR 0086): the edge's
3D samples, shared by every face using it, are inverted onto the surface,
unwrapped along periodic axes, and bounded per span by the chord bound `d`,
the parameter chord's image `e` and the lens between the edge's preimage and
the chord. What was missing for a relation is the inversion itself: a
certified, unique `(u, v)` for a point, or a refusal by name.

## Decision

We invert a sample onto a relation along its swept parameter, with a
certified one-dimensional search, and hand the parameters to the inverter
the atomic faces use; nothing downstream of the inversion changes.

**D1. Generators.** Each relation is a one-parameter family of generators
along a swept parameter `t`: the circle about the axis through `c(v)`
(revolution, `t = v`), the line `c(u) + s d` (linear extrusion, `t = u`) and
the ruling from `A(u)` to `B(u)` (ruled, `t = u`). The distance `g(t)` from a
sample `P` to the generator at `t` changes no faster than the swept curve
moves: `|g(s) - g(t)| <= L |s - t|`, `L` the curve's certified speed bound
(`curve_derivative_bounds3`; the larger of the sides' for a ruled surface).

**D2. The certificate.** A branch-and-bound over the bracket of `t` (a swept
B-spline's domain, `[0, 2 pi]` for a conic, `[0, 1]` for a ruled surface,
cut into 64 cells up front) bounds a cell `[a, b]` with midpoint `m` from
below by `g(m) - L (b - a) / 2`. It drops a cell whose bound exceeds
`min(tau, best + slack)`, `tau` the linear tolerance, `best` the least `g`
evaluated and `slack = 1e-3 tau`, and keeps a cell as a leaf once `L (b - a)
/ 2 <= slack`. A cell split from a coarse one takes, where smaller, the
generator's speed at its middle plus the coarse cell's acceleration bound
over half its width, so the leaves follow the curve's local speed rather
than its fastest (an ellipse's flat side, a ruled side's slow rim). The
leaves cover every `t` with `g(t) <= min g + slack` and lie within `min g +
3 slack`. A swept line is unbounded and is handled exactly:
extruded, it sweeps a plane, inverted in closed form; revolved, its bracket
is cut per sample to where both the height and the distance from the axis
can be within `tau` of the sample's.

**D3. Uniqueness.** The leaves must merge into one interval, or on a conic's
periodic bracket into two touching both ends, which meet across the seam.
Gauss-Newton then sharpens `t` inside the interval to `t*`, and the other
parameter is read in closed form (`u = atan2` about the axis; `v` by
projection onto the line or the ruling, clamped to the ruling). Every
generator over the interval passes within the slack of the sample. With
`o*` the other parameter read at `t*`, the points `G(t) = S(t, o*)` over the
interval must stay within `tau` of `S(t*, o*)`, or the surface runs along a
fold there and the sample does not fix `t`. That spread is bounded locally,
at `t*`: with `r` the interval's farthest end from `t*`, by `r |G'(t*)| + r^2
sup|G''| / 2`, where `G` moves as the swept curve for a revolution (turned by
a fixed `u*`) and an extrusion (shifted by a fixed `v* d`) and as `(1 - v*) A
+ v* B` for a ruled surface, `sup|G''|` certified over the interval. Where a
tangent break of the swept curve (a knot of full multiplicity, a ruled
side's corner) lies inside the interval, the bound is the certified speed
over the interval times `r` instead: the speed jumps there. For the same
reason the cover's coarse cells are cut at every tangent break, so no cell's
speed bound, nor its Taylor step from the middle, spans one. A ruling pivoting about a slow rim,
or an ellipse's flat side, is not a fold: the generator barely moves near
the sample, however fast its far end or its vertex moves. The residual
`|S(u, v) - P|`, measured by the relation's own evaluator, must be within
`tau`. The parameters returned depend on the relation, the tolerance and the
point only.

An edge whose swept parameter moves more than half the domain between two
samples of a closed swept curve without a period (a closed B-spline, ruled
sides that both close) crosses its seam, where its parameters jump by the
domain; it is refused before its parameter ring is triangulated.

**D4. Refusals.** Each is named:

| Case | Error |
| --- | --- |
| a sample within `tau` of the axis of revolution | `Degenerate`: runs into the axis, where every `u` names the point |
| a sample within `tau` of a ruled apex | `Degenerate`: runs into the apex |
| a sample on a ruling shorter than `tau` | `Degenerate`: runs into a collapsed ruling, where every `v` names the point |
| no generator within `tau`, or a residual above it | `Degenerate`: leaves its surface by at least a certified floor |
| a cover in two parts, a spread led by the first-order term (a fold), or an edge across a closed swept curve's seam | `Degenerate`: no unique parameters (a closed swept curve's seam, a surface covered twice, a fold or a surface meeting itself) |
| a spread led by the second-order term | `Degenerate`: no unique parameters, its generators turn about it (near an apex or a slow rim of the swept curve, where it stops and turns back) |
| more than 100,000 distance evaluations | `BudgetExceeded`: certified inversion of an edge sample onto a surface relation |
| a line extruded along its own direction | `Degenerate`: it is a line, not a surface |
| a non-finite sample, or a swept line with a zero or non-finite direction | `InvalidInput` |
| finite input whose arithmetic comes to no finite number (a distance, a bracket, a residual) | `InvalidInput`: met a non-finite value |

A profile that runs along the extrusion direction over a stretch exhausts
the work cap. The quarter-period guard, the edge-sample budget and the
relation evaluator's own domain refusals stay as they are.

**D5. The bound.** No new deviation term. A mesh vertex is the edge's 3D
sample itself, and its distance from `S(u, v)` is the vertex offset the
triangle bound adds already. The lens reads a relation's closed-form
tangents and open domain (a swept B-spline's domain, the ruled unit square)
through the face's support, and reports the face unbounded by name where it
cannot show the surface invertible (near a fold or the axis).

**D6. Scope.** Relation poles (the axis, a ruled apex, a collapsed ruling)
stay refused: #463. A surface covered twice (a full circle revolved about its
diameter, a line perpendicular to the axis, a profile extruded in its own
plane) and the seam of a closed swept B-spline or of closed ruled sides are
refused, never resolved by a guess, as ADR 0086 refuses a clamped closed
B-spline surface's seam; ADR 0103 gives neither a period. The inverter stays
crate-private, as the relation evaluators do.

## Alternatives considered

| Option | Why not |
| --- | --- |
| Convert the relation to a NURBS surface and use `invert_surface_certified` | A two-dimensional search; a conic's rational parameter is not the angle ADR 0103 fixes; a full turn's seam reads ambiguous; an extrusion needs a finite `v` range it does not have; the ruled sides would need unifying. |
| Reduce each kind to `project_curve3_certified` | Only B-spline swept curves and profiles in a meridian plane; nothing for ruled surfaces; three mechanisms instead of one. |
| A two-dimensional search over the ruled unit square (#454's sketch) | `v` is a projection onto the ruling in closed form, so the search is one-dimensional in `u` like the others. |
| Lift the vertex locator's cover (`brep/trim.rs`) into a shared module now | Its mutation probe anchors in that body; refactoring a certified path inside a capability landing mixes risks. |
| Resolve a double cover by continuity along the edge | A guess (ADR 0086); a trimmed relation's range (#455) is the right restriction. |

## Consequences

**Positive**

- Faces over revolutions, linear extrusions and ruled surfaces whose edges
  carry 3D curves only mesh with a certified bound, sharing samples with
  their neighbours: a revolved wall closes with its planar caps.
- Where the relation's derivative bounds are an atomic surface's (a line
  revolved about a parallel axis, a circle extruded along its normal), the
  bound is the atomic cylinder's to rounding.

**Negative / costs**

- As on an atomic surface, an inverted edge's terms are sampled to the
  budget and the face adds its triangles' term, so an exact face without
  pcurves reports above the budget: on the tested fixtures at 1e-3, a
  cylinder sheet 1.09e-3 (on its atomic cylinder as on its relation), a
  flat ellipse beside its vertex 1.06e-3, and a frustum patch, an arc ruled
  against a line and quarter arcs ruled 1:0.03 1.6e-3, so up to 1.6 times
  the budget; the tests hold every such face within twice it, and its
  bound against the exact surface.
- The edge sampler refines a piece of an edge uniformly, until the lens can
  show the surface invertible along every span. Where the surface's speed
  across the edge varies by orders of magnitude (a ruling pivoting about a
  slow rim), the slowest stretch sets the spacing of the whole edge, and the
  interior is refined around those samples: quarter arcs ruled 1:0.03 mesh
  in about 49,000 triangles without pcurves against about 1,800 with them,
  1:0.01 in about 170,000, and 10:0.1 and 10:0.03 run into the curved-face
  vertex budget (`BudgetExceeded`, after 38 s and 220 s) although every
  sample inverts. Sampling graded by the local speed would lift this; it
  changes the edge sampler that atomic faces share, so it is left to #497.
  Measured on this branch (optimised test build, 1e-3, no pcurves):

  | Face | Before (fold rule by the coarse speed) | After (local spread, local arc, local cover speed) |
  | --- | --- | --- |
  | ellipse 3 x 0.1 over its vertex, edges 0.1 tau off | refused as a fold | 350 triangles, 0.17 s, 8.5e-4 |
  | ellipse 1.2 x 0.1 over its vertex, edges 0.7 tau off | refused as a fold | 770 triangles, 0.12 s, 1.0e-3 |
  | ellipse 1.2 x 0.1 over its vertex, exact | 32,768 triangles, 1.3 s, 6.2e-7 | 4,096 triangles, 0.22 s, 4.0e-5 |
  | ellipse 1.2 x 0.1 beside its vertex, exact | 3,626 triangles, 3.1 s, 1.0e-3 | 1,126 triangles, 0.28 s, 1.1e-3 |
  | quarter arcs ruled 10:0.03, exact | refused as a fold | `BudgetExceeded` (curved-face vertices), 220 s |
  | quarter arcs ruled 10:0.1, exact | `BudgetExceeded` (curved-face vertices), 54 s | the same, 38 s |
  | quarter arcs ruled 1:0.03, exact | 49,366 triangles, 5.6 s, 1.6e-3 | the same |

  The ellipse rows owe their change to the span's arc bound
  (`brep/trim.rs`), now `|C'(m)| h + sup|C''| h^2 / 4` where that is below
  the curve's speed bound times `h`: a conic's speed bound is its largest
  semi-axis over any span, which made the lens's box far too wide along an
  ellipse's flat side. That bound is shared with atomic faces, which it only
  tightens.
- Double covers are common in authored data: a dome revolved from an arc
  centred on the axis is refused when the arc is stored as a full circle,
  because a relation carries no swept range. Adapters must lower trimmed
  arcs to B-spline arcs until #455.
- A kinked or degree-1 profile is not smooth across the lens's box, so such
  a face meshes but reports unbounded by name.
- A full turn bounded by its two rims with no seam edge is not triangulated,
  on an atomic cylinder as on a relation.

**Follow-ups / risks to watch**

- #463: faces reaching a relation's pole.
- #497: grade B-rep edge sampling by local speed (`TODO(#497)` at
  `sample_edge`).
- #455: a trimmed relation's range should restrict the bracket of `t`.
- A period for closed swept B-splines and closed ruled sides would remove
  their seam refusal; it amends ADR 0103.
- One shared one-dimensional cover for `brep/trim.rs` and this inverter.

## Relation to existing code

`crates/execution/compile/src/surface_relation/invert.rs`
(`RelationInverter`), `src/surface_relation.rs` (`RelationSurface::tangents`,
`open_domain`), `src/brep/support.rs` (`Support::partials`, `open_domain`),
`src/brep/inversion.rs` (`SurfaceInverter` over an atomic or a relation
target, its seam check, `FaceBound::lens`) and `src/brep/trim.rs` (the
span's local arc bound). Tests: the module's unit tests,
`tests/relation_edge_inversion.rs` (fixtures shared with
`tests/relation_faces.rs` in `tests/common/mod.rs`) and
`tests/ruled_surface.rs`. Mutation probe:
`scripts/probe_relation_inversion_mutants.py`. The ledger rows B5 and C8
record the level.
