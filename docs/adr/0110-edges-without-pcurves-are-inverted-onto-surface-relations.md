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
/ 2 <= slack`. The leaves cover every `t` with `g(t) <= min g + slack` and lie
within `min g + 3 slack`. A swept line is unbounded and is handled exactly:
extruded, it sweeps a plane, inverted in closed form; revolved, its bracket
is cut per sample to where both the height and the distance from the axis
can be within `tau` of the sample's.

**D3. Uniqueness.** The leaves must merge into one interval, or on a conic's
periodic bracket into two touching both ends, which meet across the seam.
The interval must also hold no generators more than `tau` apart (its width
times the largest `L` over it): generators that far apart all passing within
the slack of the sample mean a fold of the surface, where the sample does
not fix `t`. Gauss-Newton then sharpens `t` inside the interval, the other
parameter is read in closed form (`u = atan2` about the axis; `v` by
projection onto the line or the ruling, clamped to the ruling), and the
residual `|S(u, v) - P|`, measured by the relation's own evaluator, must be
within `tau`. The parameters returned depend on the relation, the tolerance
and the point only.

**D4. Refusals.** Each is named:

| Case | Error |
| --- | --- |
| a sample within `tau` of the axis of revolution | `Degenerate`: runs into the axis, where every `u` names the point |
| a sample within `tau` of a ruled apex | `Degenerate`: runs into the apex |
| a sample on a ruling shorter than `tau` | `Degenerate`: runs into a collapsed ruling, where every `v` names the point |
| no generator within `tau`, or a residual above it | `Degenerate`: leaves its surface by at least a certified floor |
| a cover in two parts, or along a fold | `Degenerate`: no unique parameters (a closed swept curve's seam, a surface covered twice, a fold or a surface meeting itself) |
| more than 100,000 distance evaluations | `BudgetExceeded`: certified inversion of an edge sample onto a surface relation |
| a line extruded along its own direction | `Degenerate`: it is a line, not a surface |
| a swept line with a zero direction | `InvalidInput` |

A sample on a fold at a fine tolerance runs into the work cap before the
fold is seen. The quarter-period guard, the edge-sample budget and the
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
  budget and the face adds its triangles' term: a cylinder sheet without
  pcurves reports 1.09e-3 at a 1e-3 budget, on its atomic cylinder as on its
  relation.
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
- #455: a trimmed relation's range should restrict the bracket of `t`.
- A period for closed swept B-splines and closed ruled sides would remove
  their seam refusal; it amends ADR 0103.
- One shared one-dimensional cover for `brep/trim.rs` and this inverter.

## Relation to existing code

`crates/execution/compile/src/surface_relation/invert.rs`
(`RelationInverter`), `src/surface_relation.rs` (`RelationSurface::tangents`,
`open_domain`), `src/brep/support.rs` (`Support::partials`, `open_domain`)
and `src/brep/inversion.rs` (`SurfaceInverter` over an atomic or a relation
target, `FaceBound::lens`). Tests: the module's unit tests,
`tests/relation_edge_inversion.rs` (fixtures shared with
`tests/relation_faces.rs` in `tests/common/mod.rs`) and
`tests/ruled_surface.rs`. Mutation probe:
`scripts/probe_relation_inversion_mutants.py`. The ledger rows B5 and C8
record the level.
