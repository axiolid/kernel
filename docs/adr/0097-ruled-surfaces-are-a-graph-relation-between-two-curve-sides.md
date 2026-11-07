# 0097 — Ruled surfaces are a graph relation between two curve sides

- **Status:** Accepted
- **Date:** 2026-11-03
- **Deciders:** Friedrich Schrödter
- **Supersedes:** —
- **Amends:** [0085](./0085-geometry-graph-wire-format.md) (format 1.5)
- **Amended by:** [0103](./0103-b-rep-faces-on-swept-and-ruled-surface-relations.md) (a B-rep face on a ruled surface)

## Context

Issue #402 (consumer: openbimrs/revit#108). Revit's `RuledSurface` joins two
profile sides: each is a curve over an interval (`m_endParams`, possibly
reversed) or a point (an apex), and the surface is `S(u, v) = (1 - v) A(u) +
v B(u)`. The consumer holds lines, arcs, ellipses and Hermite splines for the
sides, and no way to hand the surface over.

Two obvious routes do not work. Converting the surface to a NURBS fails for a
conic: a circle's angle parameter is not rational, so a ruled surface between
two arcs at equal angle is not a rational B-spline surface, and the equal-`u`
pairing is part of the surface's meaning. A new `axiolid_surface::Surface`
variant would make a two-curve construction an atomic value that owns its
curves inline, against ADR 0009 and the surface crate's README, which keep
relations between nodes in the graph.

## Decision

`SurfaceRelation::Ruled { first, second }` in `axiolid-model`, each side a
`RuledSide`: `Curve { curve, interval }`, an atomic `Curve3` node (no
instance, no curve relation) over an oriented interval whose `u` reads
`interval.start + u (interval.end - interval.start)`, or `Apex(point)`. The
builder refuses a non-finite or empty interval, a non-finite apex, two apexes
and one curve over one interval (`GraphError::InvalidRuledSurface`), and
reads the sides' references first, then second. Wire format 1.5 carries it
(ADR 0085 note).

`ReferenceMeshCompiler` meshes it as an open sheet, `MeshClosure::Surface`,
with a certified `DeviationBound` under the new `DeviationPath::RuledSurface`;
the evaluation stays `pub(crate)` in `axiolid-mesh-compile` (`ruled.rs`). A
side is a line, a circle, an ellipse or a B-spline (rational too). A conic
may run at most one turn, a B-spline stays in its domain; any other family is
`UnsupportedInput`. The exact compiler and a B-rep face on any surface
relation refuse with typed errors (`UnsupportedInput`), no longer an
`InvalidInput` for the face.

**The bound.** The grid has columns `u_i` and `m` rows; every vertex
`P_ij = (1 - v_j) a_i + v_j b_i` is an exact surface point. A cell is two
triangles split along the `(i, j)`-`(i+1, j+1)` diagonal. With `l_A`, `l_B`
the sides' linear interpolants over the cell and `L = (1 - v) l_A + v l_B`:

1. *Curve term.* `S - L = (1 - v)(A - l_A) + v (B - l_B)` and `|A - l_A| <=
   h_u^2 / 8 sup|A''|`, with `A'' = (t1 - t0)^2 c''`. The cell is within
   `h_u^2 / 8 max_sides (t1 - t0)^2 sup|c''|`, `sup|c''|` from
   `curve_derivative_bounds3`. This is the Taylor form: a circle's sagitta
   bounds the distance to the chord's line, not the gap between points of
   equal parameter, so `chord_bound3` does not apply.
2. *Twist term.* `L` is bilinear over a sub-cell with twist `T = h_v T_i`,
   `T_i = (b_{i+1} - b_i) - (a_{i+1} - a_i)`, first order in the cell's
   width. The triangulation differs from
   `L` by `w (1 - s) T` under the diagonal and `s (1 - w) T` over it, each at
   most `|T| / 4`. The cell adds `h_v |T_i| / 4`.
3. *Apex.* `L` is the planar triangle `(a_i, a_{i+1}, apex)`, which the mesh
   covers: one row, twist term 0.

The report is the worst cell times `1 + 1e-12`. Columns: the sides' tangent
breaks (`continuity_breaks3`) are forced, each piece is cut so the curve term
fits half the chord, every cell is split by the corners' twist estimate, and
the rows are chosen so the twist term fits the other half. A cell whose
two end rulings both vanish would be a segment hiding the lens between the
sides, so it is halved until a midpoint ruling does not vanish. If the
midpoint vanishes too and the certified gap `h^2 / 8 (K_A + K_B)` is at most
half the chord (or both sides are straight), the sides are the same curve
there and the sheet is refused as `Degenerate`, wholly or in part. A column
whose ruling has no length is one vertex, as an apex is; a triangle of
three distinct vertices and exactly zero area is a collapse, also refused
as `Degenerate`, so no mesh carries a degenerate triangle. Refusal of
coincident sides is guaranteed where the rulings vanish bitwise or a
triangle collapses; otherwise a near-coincident sheet meshes with a sound
bound. Caps: 4096 cells per direction, 2^20 in all, `BudgetExceeded`.

## Alternatives considered

| Option | Why not |
| --- | --- |
| Convert to a NURBS surface | A conic's angle parameter is not rational; the equal-`u` pairing would change. |
| `axiolid_surface::Surface::Ruled` | An atomic value owning two curves inline; ADR 0009 keeps relations in the graph. |
| The issue's flat shape with a `Point3` node as apex | A point node has no interval, so the side would carry a meaningless one. |
| A public ruled evaluator in `evaluate` or `construct` | No consumer needs one yet; a public API is a commitment (ADR 0078). |
| `chord_bound3` for the curve term | It bounds the distance to the chord's line, not the same-parameter gap `L` needs. |

## Consequences

**Positive**

- A Revit ruled surface lowers to one node and meshes with a bound the
  report checks against a closed-form oracle.
- The model stays the input seam: no closure change, no new crate.

**Negative / costs**

- The twist term is conservative for strongly twisted strips (sound, within
  the chord budget, with more rows).
- Only line, circle, ellipse and B-spline sides mesh; others are refused
  until a bound exists for them.

**Follow-ups / risks to watch**

- Swept and revolved relations are still carried but refused by the
  compilers as sheets; ADR 0103 reads them as a B-rep face's support.
- A ruled surface as a B-rep face support is read by ADR 0103 when every
  edge use carries a pcurve; it is refused otherwise.

## Relation to existing code

`crates/representations/modeling/graph/src/surface_relation.rs`
(`RuledSide`, `SurfaceRelation::Ruled`), `validation.rs`,
`crates/execution/compile/src/ruled.rs`, `tests/ruled_surface.rs` in both
crates and `scripts/probe_ruled_surface_mutants.py`. The ledger row B5
records the level.
