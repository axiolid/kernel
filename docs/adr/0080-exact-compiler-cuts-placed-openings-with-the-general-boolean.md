# 0080 — The exact compiler cuts placed openings with the general boolean

- **Status:** Accepted
- **Date:** 2026-10-03
- **Deciders:** Friedrich Schrödter
- **Supersedes:** —

## Context

Issue #228: `ReferenceExactCompiler` compiled `SolidOperation::Boolean`
only through `boolean_prisms_exact`, on unplaced sharp rectangles extruded
along `+z`. Building models place their operands: a consumer scan of 40
IFC models found 1243 opening bodies, about 98% of them extrusions, and
door and window openings usually run perpendicular to the wall's own
extrusion, so the coaxial prism reduction (one plan arrangement crossed
with height intervals) cannot express them. Slab shafts run parallel.

The general exact boolean (`axiolid-brep-boolean`, ADR 0075) builds the
difference of any two exact B-reps with analytic faces, and since #223
the compiler places an exact B-rep with `ExactBRep::transformed`.
`axiolid-mesh-compile` did not depend on `axiolid-brep-boolean`; taking
the edge adds it, and only it, to the crate's production closure
(`axiolid-measure`, `axiolid-nurbs` and `axiolid-evaluate` were already
there through `axiolid-construct`). No declared closure profile
(`docs/architecture/closure-profiles.toml`) carries
`axiolid-mesh-compile`, so none changes.

## Decision

The exact compiler cuts a placed extrusion with placed extrusions through
the general exact boolean, and `axiolid-mesh-compile` depends on
`axiolid-brep-boolean` (an allowlisted internal edge in its
`[package.metadata.axiolid]`).

- Dispatch: two unplaced sharp rectangles along `+z` keep the prism path
  (integer-exact plan decisions, every operator). Otherwise a difference
  whose subject is a placed extrusion or an earlier difference, and whose
  tool is a placed extrusion, compiles both operands exactly (placements
  through `transformed`) and runs `axiolid_brep_boolean::boolean`. Nested
  differences compose, one opening at a time.
- The coaxial arc-prism path is not used for placed operands: it takes a
  single-ring subject in a shared frame, so it would serve a slab's first
  shaft and no opening after it, and a second route would be a second set
  of results to keep in agreement.
- Refused by name: unions and intersections of placed operands, operands
  that are not extrusions, a tool that is itself a boolean, scaled or
  sheared placements, and every `BooleanError` the general boolean
  raises. An emptied subject is `GeomError::Degenerate`.
- Independently placed operands meet faces that agree only up to
  rounding. The general boolean reads them within tolerance where its
  exact predicates would otherwise see the residue: a section along an
  edge, a tangent crossing, near cuts and boundary splits, and a plane
  parallel or perpendicular to, or touching, a cylinder.

### What a result guarantees

- **Source of the tolerance.** The caller's `ExecutionOptions` tolerance,
  passed through unchanged: its linear part `eps` for every distance, its
  angular part `alpha` for the coincidence and contact of normals and
  axes. No decision uses a built-in constant; bookkeeping in face
  parameters (welding one vertex's two evaluations, ordering pieces at a
  vertex) uses slacks far below any tolerance and refuses what it cannot
  order.
- **Exact decisions.** Which section two supports have (except the
  plane/cylinder reading below), where a section crosses the surface next
  to an edge, whether a point lies in a face (certified) or in a solid
  (ray parity over exact intersections).
- **Within-tolerance decisions.** Coincident supports, contact, a point on
  an edge, a section along an edge (proved over the whole edge for lines
  and circles, never read for other pairs), a tangent crossing recovered
  from the edge, cuts and splits within `eps` of each other, and a plane
  parallel (`|n . a| * extent <= eps` over the faces' common box),
  perpendicular (`r sin(theta) <= eps`) or touching (`eps`) a cylinder.
- **Guarantee.** With no within-tolerance decision the result is the exact
  boolean of the operands. Otherwise it is the exact boolean of operands
  whose faces moved by at most `eps` (turned by at most `alpha` for a
  direction decision), every surface and curve exact for them. A reading
  that no single such perturbation explains is refused by name (cuts
  chained over more than `eps`: `BooleanError::NearCoincidence`). Features
  further apart than `eps` go through the exact predicates: a skin or a
  sliver ten tolerances thick is kept (tests in
  `crates/execution/compile/tests/exact_placed_boolean.rs`).

## Alternatives considered

| Option | Why not |
| --- | --- |
| Bring both operands into the subject's frame and keep the 2D arc path | Only parallel extrusions; perpendicular openings, the common case, stay refused. |
| Arc path for parallel cases, general boolean otherwise | Covers the first shaft of a slab only; two result shapes for one operation. |
| Snap operands to a shared lattice before the boolean | Moves geometry the caller gave; the tolerance belongs where decisions are made. |

## Consequences

**Positive**

- Walls and slabs with placed openings compile exactly: perpendicular and
  parallel, rectangle, circle and line-and-arc profiles, through, blind,
  flush and touching, several per body; volumes and certified distances
  match closed forms.

**Negative / costs**

- The general boolean is slower than the prism path: tens of milliseconds
  per opening in release. Exact section and classification over every
  face pair made the cost grow with the openings already cut; sound
  bounding boxes (enlarged by the tolerance, so a pair that touches is
  never skipped) now skip face pairs, edges and ray tests that cannot
  matter. `tools/benchmark/benches/exact_openings.rs`, a 6 m wall under a
  general placement losing `n` windows, validated against the exact
  volume, on a shared 20-core machine under load: 1 window 42 ms to 18 ms,
  3 windows 224 ms to 85 ms, 10 windows 1.67 s to 0.28 s.

**Follow-ups / risks to watch**

- Half-space clipping (`IfcHalfSpaceSolid`,
  `IfcPolygonalBoundedHalfSpace`) is the next step.

## Relation to existing code

- `crates/execution/compile/src/exact/boolean.rs`
- `crates/algorithms/construction/brep-boolean/src/section.rs`,
  `crates/algorithms/construction/brep-boolean/src/split.rs`
- Tests: `crates/execution/compile/tests/exact_placed_boolean.rs`,
  `crates/algorithms/construction/brep-boolean/tests/openings.rs`; probe
  `scripts/probe_placed_boolean_mutants.py`.

## Amendment 2026-10-03: half-space clipping (#234)

A boolean whose tool is a half-space is dispatched first, to
`crates/execution/compile/src/exact/clip.rs`: a difference or intersection
whose right operand is a `HalfSpace` (`IfcHalfSpaceSolid`) or a
`BoundedHalfSpace` (`IfcPolygonalBoundedHalfSpace`, polyline boundary),
read through rigid placements. The subject is anything the rest of the
dispatch compiles from placed extrusions, so clips compose with each other
and with openings, in either order.

- **Semantics.** The mesh compiler's: `agreement` selects the normal side;
  a bounded half-space is that side within the prism of its boundary
  swept along the plane normal, the boundary framed by its placement's
  axes projected into the plane and its origin projected onto it.
- **Finite tool, general boolean.** The half-space becomes a prism standing
  on the plane over a sound envelope of the subject (the box of a placed
  extrusion's edges, carried through placements; a difference or clip
  keeps its subject's), plus a margin `m` (a quarter of the envelope's
  diagonal plus four tolerances). Its footprint is the envelope's
  projection widened by `m` (unbounded) or the boundary (bounded); it
  reaches `m` past the envelope. Since the subject lies in the envelope,
  `S - P = S - H` and `S ∩ P = S ∩ H` for every `m > 0`, and the faces `m`
  moves stay more than the tolerance from the subject, so no decision of
  the general boolean involves them: the result does not depend on `m`.
  A subject whose envelope lies wholly off the half-space (or, unbounded,
  wholly in it) is decided without the boolean.
- **Not a direct plane clip.** A dedicated clip of an exact B-rep by a
  plane would be a second section, split and classification to keep in
  agreement with the general boolean, for a tool of six faces whose other
  five the boolean's bounding boxes already skip.
- **Refused by name:** a union with a half-space (unbounded), a
  half-space as the subject, a subject not built from placed extrusions,
  a bounded half-space whose boundary is not a polyline or whose plane is
  placed by an instance, scaled or sheared placements, and the general
  boolean's refusals. A curved base surface cannot be expressed: the
  model's half-space carries a plane. A clip removing the whole subject
  is `GeomError::Degenerate`.
- **Known limit.** A plane within the tolerance inside a cylinder face
  (crossing it, read as touching) is refused by the general boolean as a
  point too close to classify; a plane tangent to it, or outside, keeps
  the whole solid.

Tests: `crates/execution/compile/tests/exact_half_space_clip.rs`; probe
`scripts/probe_half_space_clip_mutants.py`.
