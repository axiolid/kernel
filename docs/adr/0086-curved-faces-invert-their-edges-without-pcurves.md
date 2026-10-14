# 0086 — Curved faces without pcurves invert their shared edge samples

- **Status:** Accepted
- **Date:** 2026-10-14
- **Deciders:** Friedrich Schrödter
- **Supersedes:** —

## Context

Issue #310. `axiolid-mesh-compile` tessellated a curved B-rep face only
when every edge use in its loops carried a pcurve, and refused it with
`Unsupported { operation: Tessellation }` otherwise. Exchange formats
usually give an edge only its 3D curve: IFC's `IfcEdgeCurve` without an
`IfcPcurve` is one example. Cylindrical and B-spline faces from such
files lowered exactly into the graph and then could not compile.

The triage on #310 binds this record: don't fit a pcurve; invert the
edge's existing shared 3D samples onto the face's surface, closed form for
analytic surfaces and by the certified global search for B-spline
surfaces; unwrap periodic parameters along the ring by continuity and
refuse an edge that runs into a singularity; report a deviation that adds
a bound on how far the surface image of each parameter chord strays from
the edge on top of the edge's chord bound. Edges with pcurves keep the
current path.

## Decision

We will build a curved face's parameter boundary, where an edge use has
no pcurve, by inverting the 3D samples every face sharing the edge uses
(`crates/execution/compile/src/brep/inversion.rs`).

- The first face to reach the edge samples its 3D curve (or the segment
  between its vertices when it has none) at uniform parameters, doubling
  until each span's midpoint meets the chord budget, consecutive samples
  turn less than a quarter period on a periodic axis, and the span's
  certified terms meet the budget where that is cheap. Later faces invert
  the same points, so shared vertices stay shared.
- Analytic surfaces invert with `axiolid_evaluate::surface::invert`
  (reached through `axiolid-reference`'s re-export); a sample within the
  linear tolerance of a cone apex or sphere pole is refused by name.
  B-spline surfaces invert with `axiolid_nurbs::invert_surface_certified`
  (the periodic form when the closed axes validate as periodic), at a
  thousandth of the linear tolerance first and at the tolerance itself for
  a point the first pass finds just off the surface; ambiguous,
  off-surface and unresolved samples are refused. Gauss-Newton sharpens
  the parameters inside the proven enclosure.
- The reported deviation adds per span `d + e + M W^2 / 2`: the edge's
  chord bound, the surface image of the parameter chord, and the lens
  between the edge's preimage and that chord, whose width `W` is bounded
  through a lower bound on the surface derivative's singular values over a
  box the preimage provably stays in. Where that bound cannot be had the
  face is reported unbounded by name.

A face on a declared plane with an edge that follows a curve takes the
same path with the plane as its surface. The planar path triangulates a
face's vertices only, so a cap bounded by one closed circle had no
polygon at all; through the shared samples it meets the wall across its
rim by index, and its bound is the edges' chord bounds plus each
sample's distance from the plane (the plane bends nowhere, so the image
and lens terms are zero). A face with no declared surface keeps the
planar path.

`axiolid-mesh-compile` gains a direct dependency on `axiolid-nurbs`
(`allowed-internal-dependencies`). `axiolid-evaluate` was already in its
closure through `axiolid-reference`. No closure profile in
`docs/architecture/closure-profiles.toml` contains `axiolid-mesh-compile`,
and every profile that may contain it already contains `axiolid-nurbs`
and its dependencies, so no `expected_internal` list changes.

## Alternatives considered

| Option | Why not |
| --- | --- |
| Fit a 2D pcurve to the inverted samples | A fit is a second approximation with its own error to certify, and the shared samples would no longer be the boundary's vertices. |
| Invert onto the surface per face, sampling the 3D curve independently | Two faces would place different vertices on one edge, and the mesh would open along it. |
| Keep refusing | IFC exports would keep failing for every curved face. |

## Consequences

**Positive**

- Curved faces whose edges carry only 3D curves compile, watertight
  along shared edges, with a certified deviation.

**Negative / costs**

- `axiolid-mesh-compile` depends on `axiolid-nurbs` (and through it on
  `axiolid-exact` and `axiolid-numeric`).
- A B-spline sample where four Bézier patches meet exhausts the certified
  search's work budget and is refused (`BudgetExceeded`).

**Follow-ups / risks to watch**

- An edge shared with a face that sampled it through a pcurve is inverted
  at those samples; its 3D chord bound is then known only for a straight
  edge, and a curved one reports the face unbounded by name.
- A periodic B-spline face whose boundary crosses the seam outside the
  surface's domain is refused.

## Relation to existing code

- `crates/execution/compile/src/brep.rs` (`curved_boundary`,
  `SharedEdge`, `edge_samples`)
- `crates/execution/compile/src/brep/inversion.rs`
- `crates/execution/compile/tests/edge_inversion.rs`
