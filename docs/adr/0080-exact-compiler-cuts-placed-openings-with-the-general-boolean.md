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
  the whole solid. (Lifted by the #243 amendment below.)

Tests: `crates/execution/compile/tests/exact_half_space_clip.rs`; probe
`scripts/probe_half_space_clip_mutants.py`.

## Amendment 2026-10-03: exact first, and a report (#236)

A consumer could not tell an exact difference from one built within
tolerance, and `Tolerance::ZERO` failed for every placed difference, even
for openings placed by matrices with entries `0` and `+-1`. So every result
had to be treated as perturbed.

- **Why zero failed.** Naming a point's parameters on the curve or surface
  it was evaluated from (`axiolid_evaluate::surface::locate`, `locate3`)
  checks the `f64` round trip against the linear tolerance; at zero any
  rounding residue fails it (`BooleanError::Evaluation`). That is
  bookkeeping, not a decision. It, and welding two evaluations of one
  vertex, now allow the rounding of the operands' coordinates: `2^-40` of
  their extent. A residue that small is one exact point evaluated twice;
  only a larger one is a decision within tolerance.
- **Exact first.** Readings about the operands' own surfaces (coincident
  supports, a plane parallel or perpendicular to or touching a cylinder)
  ask an exact predicate on the operands' numbers first (dyadic arithmetic,
  `axiolid-exact`, already in the crate's closure through `axiolid-nurbs`;
  the direct edge is allowlisted). Exactly coplanar, parallel and
  perpendicular faces are decided without the tolerance at any tolerance.
- **Report.** `axiolid_brep_boolean::boolean_with_report` returns a
  `BooleanReport`: per kind of within-tolerance reading that fired, the
  furthest it moved (`linear`) and turned (`angular`) the operands.
  `ReferenceExactCompiler::compile_exact_with_report` (and a batch form)
  merges the reports of every general boolean and clip beneath a body.
- **Guarantee, amended.** An empty report means the result is the exact
  boolean of the operands as given; at `Tolerance::ZERO` the report is
  always empty, and operands that miss coincidence by rounding are refused
  by name rather than guessed. A non-empty report bounds the perturbation
  of the existing guarantee. Placed with exact axis matrices, through,
  blind and flush openings succeed at zero and report nothing, at any
  tolerance; built from sines and cosines under a general placement, a
  flush opening's caps are read as coplanar and reported, and at zero it
  is refused. A through opening under a general rotation meets the wall
  only transversally and reports nothing: it needs no reading.
- **Not changed.** The within-tolerance readings themselves, their bounds
  and their refusals.
- **Floor exposed (#244).** Points closer than the rounding floor are one
  point and unreported, so a consumer had to mirror the private factor and
  extent. `BooleanReport::rounding_floor` now gives the floor that applied
  (`ROUNDING_FACTOR` × `BooleanReport::extent`, which states the extent's
  definition), on exact reports too. Merged reports, and the compiler's
  per-body reports, keep the largest floor beneath the body; a rigid
  instance keeps it unchanged.

Tests: `crates/algorithms/construction/brep-boolean/tests/report.rs`,
`crates/execution/compile/tests/exact_boolean_report.rs`; probe
`scripts/probe_placed_boolean_mutants.py` (report mutants).

## Amendment 2026-10-03: one plane/cylinder contact for every face pair (#243)

A placed difference whose round hole touches a planar face (an I-beam's
web hole tangent to the flange, consumer-reported) was refused ("split
face pieces do not close", "a point too close to classify", "do not
sew"), and so was the #234 known limit above. Both are a plane touching a
cylinder along a ruling, which has to be sewn into every face it reaches.

- **Why it failed.** The plane/cylinder face pair read the contact within
  tolerance as one ruling (`PlaneTouchesCylinder`). Every other face pair
  the two surfaces meet in kept its exact roots: the plane's chord across
  each cap disk, the circle a face across the axis cuts from the cylinder,
  and the plane's own edges cross the cylinder twice, `2 sqrt(2 r d)`
  apart (`d` how far the plane reaches in: a rounding residue, or a
  fraction of the tolerance). That is far more than the tolerance, so the
  split faces disagreed and did not sew. The readings were also per face
  pair, over each pair's own common box. Two bookkeeping faults made it
  worse: cuts merged across a circle's parameter origin averaged to the
  opposite side, and an imprinted edge (a flush hole's cap against the
  web) kept the two roots apart, leaving a sliver too short to classify.
- **Decision: one reading, placed on the ruling.** A plane/cylinder pair of
  supports is read once (parallel, perpendicular, touching), over the
  largest common box of any pair of their faces. A pair read as touching
  within tolerance stands for the plane moved by `d` onto the cylinder;
  wherever a curve on one of the two is cut by the other, the cut is where
  the moved plane puts it, on the contact ruling: a conic of the cylinder
  where the ruling pierces its plane (on the cylinder exactly), a line of
  the plane where it passes closest to the ruling (at most `d` off,
  `TangentCrossing`). Where the given and the moved plane disagree on a
  side (cylinder points beyond the plane, plane points inside the
  cylinder) no point decides whether a piece lies in a face or a region
  in a solid. The result is the exact boolean of operands whose plane
  moved by at most the reported distance.
- **Refused by name.** A curve whose meeting with the ruling this reading
  cannot place (a curve of the plane that is no line, a curve of the
  cylinder that is no ruling or conic) is
  `BooleanError::UnsupportedContact`, never sewn against roots that
  disagree with the reading.
- **Exact first, unchanged.** An exactly tangent pair (the exact predicate
  holds) is decided by exact double roots: with exact placements and
  dyadic sizes the difference succeeds at `Tolerance::ZERO` with an empty
  report. A section's window now reaches `2^-20` of its diagonal past the
  faces' common box, so a tool resting exactly on a face no longer leaves
  the window next to the touching point at zero tolerance.
- **Not changed.** The bound of the touching reading (`eps`), the other
  readings, and the refusal of contacts beyond the tolerance (a groove ten
  tolerances deep is cut exactly). A hole tangent to an I-beam with root
  fillets is still refused: the hole is tangent to each fillet cylinder
  where the fillet meets the flange, and that quartic with a double point
  is a section the general boolean does not build.

Tests: `crates/execution/compile/tests/exact_tangent_hole.rs`,
`crates/algorithms/construction/brep-boolean/tests/tangent_contact.rs`,
the I-beam case of `crates/execution/compile/tests/boolean_deviation.rs`;
probe `scripts/probe_placed_boolean_mutants.py` (contact mutants).

## Amendment 2026-10-03: nothing read at zero tolerance, enforced (#251)

The #236 guarantee (an empty report at `Tolerance::ZERO`) was broken by a
support reading whose `f64` measure was exactly `0` while the exact
predicate rejected it: the gate `measure <= tolerance` held as `0 <= 0`,
and the reading was taken and reported with `linear` 0. A round hole
across a wall under a general rotation hit it: its axis and the wall's
normals are the same `f64` vector, but the numbers given are not exactly
perpendicular (the frames are orthonormal only to rounding).

- **Decision.** A reading the exact predicate rejects moves or turns the
  operands by a positive amount, whatever its `f64` measure. A tolerance
  with a zero part admits none, so it never takes one: the exact answer
  stands (at zero, the general closed form cuts the hole). Both parts
  positive, nothing changes. Point readings were already sound: they
  record only residues above the rounding floor and within the tolerance.
- **Enforced structurally.** The report's session carries the caller's
  tolerance. Recording a decision at `Tolerance::ZERO`, or beyond the
  tolerance, is a debug assertion, and otherwise the entry point refuses
  it by name (`BooleanError::ToleranceExceeded`), so no future reading can
  break the guarantee silently.

Tests: `crates/algorithms/construction/brep-boolean/tests/report.rs` (a
general-placement round hole, and random general placements by every
operator, at zero); probe `scripts/probe_placed_boolean_mutants.py` (#251
mutants).

## Amendment 2026-10-04: a hole touching an I-beam's root fillets (#249)

The #243 amendment left one case refused: a round web hole touching the
flange of an I-beam with root fillets, which most rolled sections have.
The hole is tangent to each top fillet where the fillet runs into the
flange, and the two cylinders (axes perpendicular, skew by the difference
of their radii) meet in a quartic with a double point there: two loops
round the fillet cylinder, crossing at the fillet/flange edge.

- **Why it failed.** The section is traced (ADR 0077) on the fillet over
  its face's parameter box, and that face ends at the double point: the
  crossing sat on the box's edge, where the trace cannot certify the
  pieces around it, and refused. Behind that, the section touches the
  flange and web planes along the fillet's smooth edges, so cutting it
  against the adjacent surfaces meets double roots, and the face's
  pcurve would have been traced over the same box again.
- **Decision: trace the whole turn, cut smooth edges by the edge.**
  - A trace that cannot be decided in a face's box of a periodic carrier
    is taken over the carrier's whole turn. The crossing is then inside,
    the trace ends both loops at it, and the loops are cut there. On the
    carrier's own face a section piece keeps the traced curve itself as
    its pcurve.
  - A section is cut at a smooth edge of its face (the adjacent surface
    tangent to the face along it) where the edge crosses the section's
    other surface: a line against a cylinder in closed form, transversal
    at the web and an exact double root at the double point. Smoothness
    only chooses between two exact formulations of the same crossing.
  - Under a general placement the flange is read as touching the hole
    (#243), and the fillet/flange edge, a line in the flange plane, is
    cut on the contact ruling; the traced section must pass there within
    the tolerance, which it does when the residue is rounding. A traced
    section on a cylinder in a contact is cut where it crosses the plane
    through the contact ruling and the axis.
- **Refused by name.** A hole a fraction of the tolerance into or short of
  the flange meets each fillet in two arcs up to `4 sqrt(2 r eps)` from
  the contact point. Read as touching, the flange would have to join them
  there; no single move of one operand does (moving the hole onto the
  flange makes it touch the fillets as well, moving the flange alone
  leaves the fillets crossing it), so it is
  `BooleanError::UnsupportedContact`. Beyond the tolerance the exact
  predicates decide (a groove ten tolerances deep is cut exactly).
- **Exact at zero tolerance.** With exact axes and dyadic sizes the double
  point is the exact double root of the edge against the hole, and the
  report is empty at `Tolerance::ZERO`. The fillet's radius comes from
  the profile arc's bulge and is an ulp off, a residue the trace's
  crossing absorbs below its field's rounding.
- **Certified mesh.** `axiolid_evaluate::bound::chord_bound2` now bounds
  implicit pcurves cell by cell (implicit function theorem with interval
  bounds of the field's partials; a bridge cell by its cubic's control
  points), so the boolean's deviation is `Certified` where its exact
  result has traced pcurves.

Not covered by these tests: a hole flush with the web faces of a
filleted beam (its cap lies inside the fillet material and leaves a
pocket), and fillets running into a tapered flange.
Tests: `crates/algorithms/construction/brep-boolean/tests/fillet_tangent.rs`,
the filleted cases of `crates/execution/compile/tests/exact_tangent_hole.rs`
and `crates/execution/compile/tests/boolean_deviation.rs`; probe
`scripts/probe_placed_boolean_mutants.py` (fillet mutants).
