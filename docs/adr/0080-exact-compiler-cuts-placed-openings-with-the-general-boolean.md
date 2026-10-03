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
  exact predicates would otherwise see the residue (see the
  `section` module docs): a section along an edge, a tangent crossing, near
  cuts, and a plane parallel or perpendicular to a cylinder's axis.

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
  per opening in release, growing with the faces earlier openings add
  (exact section and classification over every face pair).

**Follow-ups / risks to watch**

- A face-pair bounding prefilter in `section_edges` and classification
  would make many openings per wall cheaper.
- Half-space clipping (`IfcHalfSpaceSolid`,
  `IfcPolygonalBoundedHalfSpace`) is the next step.

## Relation to existing code

- `crates/execution/compile/src/exact/boolean.rs`
- `crates/algorithms/construction/brep-boolean/src/section.rs`,
  `crates/algorithms/construction/brep-boolean/src/split.rs`
- Tests: `crates/execution/compile/tests/exact_placed_boolean.rs`,
  `crates/algorithms/construction/brep-boolean/tests/openings.rs`; probe
  `scripts/probe_placed_boolean_mutants.py`.
