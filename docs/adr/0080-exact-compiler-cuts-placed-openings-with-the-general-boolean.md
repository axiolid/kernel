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

## Amendment 2026-10-09: a skin thinner than the tolerance (#276)

A door exported `4.5e-15` m short of its wall's face left a skin that
thin. Subtracted in the wall's own frame (as `ifc-geometry` 0.13 lowers
openings) and then placed at georeferenced coordinates (about `5.6e6`, a
coordinate rounding to about `1e-9`), the mesh's two skin faces crossed,
and the volume kernel refused the solid as self-intersecting: 32 of 40
walls on one consumer model. The exact boolean read the door's end face as
the wall's (`CoincidentSupports`), but refused a door standing on the
wall's floor face that reached a rounding error past the face.

- **Reading unchanged.** The guarantee above already allows the result to
  be the boolean of operands perturbed within the tolerance, reported. A
  skin or sliver thinner than the tolerance is such a perturbation away
  from none; one thicker is kept, by both compilers.
- **Exact boolean.** Where coincident faces imprint each other's edges, a
  cut within `eps` of the imprinted edge's own end is that end
  (`MergedPoints` above the rounding floor), and cuts within `eps` of each
  other merge as on a section. The cut was the other face's boundary
  crossing the edge where the door's end face is already read as the
  wall's; kept apart, it left a piece shorter than the tolerance lying on
  the other face's boundary, which no point classifies.
- **Mesh compiler.** The mesh boolean decides on the `f64` numbers it is
  given and has no tolerance, so the reading is made before it: each
  vertex of the tool within `eps` of a subject triangle moves onto that
  triangle's plane (onto the common line or point of two or three such
  planes, by the smallest move), then each subject vertex within `eps` of
  the moved tool onto the tool's planes. A vertex moves only by at most
  `eps`, never between two parallel planes within `eps` of it, and no move
  is made that would turn or flatten an operand's triangle. Every operand
  point then moved by at most `eps`. Onto an axis-aligned plane the vertex
  lands exactly, so a host subtracted in its own frame meets the tool's
  face exactly and the mesh boolean's coplanar rule removes the skin.
- **Reported, and nothing at zero.** The largest move is a `Certified`
  contribution of the boolean path under `deviation::SNAPPED_OPERANDS`,
  carried through instances and enclosing booleans, including those
  measured against their exact result. At a zero linear tolerance nothing
  moves.
- **Placement is not checked differently.** A placed mesh is still
  refused when it intersects itself: a skin thicker than the tolerance but
  thinner than the placement's rounding is a real near-degeneracy the
  caller's tolerance does not cover.
- **Known limit.** Under a general rotation an `f64` point lands only
  within rounding of a plane, so a boolean computed far from the origin in
  a rotated frame (operands placed in world coordinates rather than in the
  host's frame) can still leave a skin at the rounding of those
  coordinates.

Tests: `crates/execution/compile/tests/thin_skin.rs` (the issue's wall,
doors short and long by `4.5e-15` to `1e-9`, at the origin and placed at
`(6e5, 5.6e6)` turned by 2.3 degrees and at `(4e5, 4e5)`; ten tolerances
kept; zero tolerance), the snap's unit tests in
`crates/execution/compile/src/snap.rs`, and
`crates/algorithms/construction/brep-boolean/tests/openings.rs`; probe
`scripts/probe_thin_skin_mutants.py`.

## Amendment 2026-10-09: the snap closes rounding only, and never refuses (#291)

The #276 snap reached the whole linear tolerance, and consumers pass one
of a millimetre. On a real model it took authored skins of 0.1-1 mm
(seven hosts lost their whole body; walls lost their skin at near-flush
slab openings), and a slab 0.3 m thick extruded downward from its top,
with 22 openings flush with both faces, that mesh-compile 0.3.15 meshed
was refused as "a void tangent to its host's face". The snap had moved
the openings' bottoms onto planes within rounding of the slab's bottom
face: a normal normalised to `0.9999999999999999`, a frame turned by a
rounded half turn (`sin pi = 1.2e-16`), a profile edge's end evaluated as
`a + (b - a)`. Some vertices landed above the face, some below, and only
those within the distance of a face's triangle moved, so the openings'
faces bent across the slab's.

- **Decision: the reach is the rounding scale.** A coordinate moves by at
  most `reach = min(eps / sqrt 3, 16 f64::EPSILON S)`, `S` the largest
  coordinate magnitude of the two operands, so a vertex moves by at most
  `eps`. Sixteen relative epsilons cover a few composed placements (each
  about `3 f64::EPSILON S`) and a 15-digit export (`4.5 f64::EPSILON`
  relative); the #276 wall's `4.5e-15` is `5.5 f64::EPSILON S` there
  (`S = 3.67`, reach `1.3e-14`). At georeferenced coordinates (`S` about
  `6e6`) the reach is about `2e-8`: a residue of an ulp or a few closes,
  a micrometre skin is fifty reaches and stays. Anything above the reach
  is authored and kept, whatever the tolerance. A skin of `1e-12`-`1e-9`
  computed near the origin is above it, kept, and placed at
  georeferenced coordinates may cross itself again: the residue is the
  rounding of the coordinates it was made at, which the snap does not
  see.
- **Decision: an exact landing only.** Only faces whose three corners
  share a coordinate exactly (an axis-aligned plane `x_k = c`) are landed
  on, and the coordinate becomes `c` itself. It moves by value: every
  vertex of the moving operand with that coordinate value gets `c`, so
  an operand's own axis-aligned face moves whole and stays planar even
  where most of it lies past the other operand (a door a micrometre below
  its wall's floor and a rounding error short of its face). A value near
  two planes, or on one and near another, stays. A face axis-aligned only
  within rounding is cut as given: no `f64` point lies on it, and a
  rounding-scale move onto it lands on either side.
- **Decision: never worse.** Where the snapped operands' mesh boolean is
  refused, or its result touches itself (#194), the boolean is cut from
  the operands as given; that result is returned with no snap reported,
  and when both are refused, the refusal of the operands as given. The
  guard is per boolean: a chain whose earlier snapped result makes a
  later boolean fail is not compared with the whole chain unsnapped.
- **The exact boolean keeps its reading, and gets the same guard.** It
  reads coincident faces within the caller's tolerance by design (above)
  and reports it; the #276 end merge is that reading's consequence on an
  imprinted edge, and the door a rounding error past the face
  (`1e-12`-`1e-9`) still needs it, so it is not narrowed to rounding.
  But it left a door `1e-6`-`5e-4` past the face under a millimetre
  tolerance unsewn, refused where before #276 it was cut. A boolean
  refused after an end merge above the rounding floor is now cut again
  with ends merged only within it, its report holding only that run's
  decisions.

Tests: `crates/execution/compile/tests/flush_through_cut.rs` (the slab in
exact half-turn and reflected frames, openings exact, reflected or under a
rounded half turn, with storey-rounded tops, in its frame and placed; the
reproduction under a rounded half turn with openings round-tripped through
world coordinates; a property over 160 random walls at `(6.1e5, 5.6e6)`
with faces flush, a few ulps off, or `1e-6`/`1e-4` off: none compiled as
given is refused after the snap, rounding residues close, authored skins
stay), `crates/execution/compile/tests/thin_skin.rs` (skins of
`1e-6`-`5e-4` kept under a millimetre; a nanometre at georeferenced
coordinates closed; the fallback with a mesh boolean that refuses the
snapped operands), the snap's unit tests, and
`crates/algorithms/construction/brep-boolean/tests/openings.rs`; probe
`scripts/probe_thin_skin_mutants.py`.

## Amendment 2026-10-13: a warped operand's slab width in a boolean's bound (#266)

Since #261/#257 a warped authored polygon face, or a warped B-rep face
without a surface, compiles with a certified bound: the width `w` of the
slab its corners span about its fit plane, within which its triangles
and every reading of the face (either diagonal, the flattened face, a
bilinear patch) lie of each other, both ways. Under a boolean the same
operand made the report unbounded: the exact compiler refuses meshes, so
nothing was measured, and a consumer compiled every authored operand
alone to find its warp.

- **Reading.** This ADR's guarantee already states a boolean's result as
  the exact boolean of operands perturbed within a reported amount. A
  warped operand's triangles are such a perturbation of it, by `w`, and
  they are exact polyhedra. So the deviation report measures a boolean
  with a mesh operand against the exact boolean of the operands
  flattened to their triangles, and reports the flattening.
- **Decision: flatten for the report only.** A crate-internal mode of the
  exact compilation (`compile_flattened_with_report`) takes an authored
  polygon or triangle mesh, and a B-rep whose faces are all planes with
  straight edges or warped faces without a surface, as the exact
  polyhedron of the triangles the mesh compiler makes of it (one planar
  face per triangle, line edges shared by index; a mesh that is not a
  closed, consistently wound two-manifold or is wound inside out is
  refused by name), as subject or tool, through rigid placements, and as
  a clip's subject. Where an operand is such a mesh it takes union and
  intersection too, since meshes are unioned as often as cut; unions of
  placed extrusions stay refused. It records the largest `w` flattened
  beneath each node, through instances and booleans. The public
  `ReferenceExactCompiler` is unchanged and refuses meshes by name: the
  triangles are one reading of a warped face, not its surface.
- **Decision: the composition.** The boolean's contribution is
  `Certified(d + w)`, `d` the certified distance from that exact result to
  the mesh, `w` the largest slab width flattened beneath the boolean, the
  sum rounded up, under the detail "measured against the exact boolean of
  operands flattened within their slab width, plus it". The search
  settles within the budget less `w`, so slab widths that fit the budget
  leave the measurement the rest and `meets_requested` follows; widths
  that do not fit settle within the budget. Each warped operand's own
  contribution travels with the boolean (as the #276 snap does), so the
  report names what was flattened. An operand's one-sided bound still
  never composes into a boolean's (#235); the slab width does because it
  is two-sided.
- **What it covers.** Relative to the flattened operands, the bound is
  exact surface to mesh, as for every measured boolean. The `+ w` covers
  the boolean of the unflattened operands where they cross transversally
  (the tests sample the boolean of the bilinear and both diagonal
  readings densely); near tangency the boolean of the unperturbed
  operands can lie further, and the guarantee is this ADR's: the exact
  boolean of operands perturbed within the reported deviation.
- **Refused by name.** Where the exact boolean cannot take the flattened
  operand (an open mesh, a B-rep with a curved face, a tool that is
  itself a boolean, a configuration the general boolean refuses, such as
  a point too close to a fine fan of triangles to classify), the boolean
  is `Unbounded` with that refusal's name, and the warp is still named;
  never unbounded without a reason.
- **Unwarped meshes.** `w = 0`: the polyhedron is the mesh's exact
  surface, and an unwarped authored box cut, unioned or intersected with
  a placed box is measured against the exact boolean of the operands as
  given (an empty report), within the budget.

Tests: `crates/execution/compile/tests/warped_boolean.rs` (a unit box with
one corner lifted 5 cm and a saddle of +-5 cm, authored and as a B-rep,
unioned with and cut by a planar box, and as the tool; at budgets of
1 mm, 5 cm and 20 cm; an unwarped pair by every operator; a chain under an
instance; a refusal by name), the polyhedron's unit tests in
`crates/execution/compile/src/exact/polyhedral.rs`; probe
`scripts/probe_polygon_faces_mutants.py` (#266 mutants).

## Amendment 2026-10-14: the snap reaches the placement's rounding (#312)

The #291 reach is the rounding of the operands' own coordinates. A host
cut in its own frame and placed by an `Instance` at georeferenced
coordinates is cut where a coordinate rounds to `1e-15` and measured
where it rounds to `1e-9`: a wall `3.65 x 0.25 x 3.67` whose doors stop
`2e-14` to `1e-10` short of one face kept a skin that thin (the #291
amendment's "residue the snap does not see"), and once placed 5600 km
out and turned to grid north its faces crossed and the volume was
refused as self-intersecting. At the origin every case measured.

- **Decision: the reach covers the scale the result is placed at.**
  `S = max(S_operands, S_placed)`, with `S_placed = P + S_operands`: `P`
  is the world coordinate magnitude the instances above the boolean add,
  in the operands' units. `Instance` and `InstanceAtStation` carry it
  down to their source with the tolerance and the chord budget
  (`instance_local_budget`): a source's `P` is its instance's `P` plus
  the largest component of the instance's translation, divided by the
  smallest column length of its linear part (a world residue is never
  read smaller than it is), rounded up to three significant bits. A
  turn's growth of a box's largest coordinate (at most `sqrt 3`) is left
  to the sixteen epsilons. Nothing placed, `P = 0` and the reach is the
  #291 one.
- **Unchanged.** The reach stays capped at `eps / sqrt 3`; only exactly
  axis-aligned faces are landed on, exactly; a snapped boolean that is
  refused is cut as given (never worse); the move is reported under
  `deviation::SNAPPED_OPERANDS`. At `5.6e6` the reach is about `2e-8`,
  so every rounding residue of the host's frame closes there and a
  micrometre skin (fifty reaches) stays authored.
- **Decision: the cache keys on `P`.** A source placed near and far is
  compiled once per scale, and the three-bit rounding lets instances at
  similar distances share one entry; a nested placement grows `P` by at
  most a quarter per level.
- **Rejected: collapsing slivers after placement.** That is mesh surgery
  after the boolean, with nothing to certify it against, where the snap
  already has a proven exact landing.

This supersedes the #291 amendment's note that a skin computed near the
origin and placed far may cross itself again, wherever an instance does
the placing. The #276 amendment's known limit stands: a boolean computed
in a rotated frame far from the origin meets faces axis-aligned only
within rounding, which no snap lands on.

Tests: `crates/execution/compile/tests/placed_skin.rs` (the reproduction:
residues `0` to `1e-10` placed and at the origin, volumes agreeing within
the residue, the reported snaps and the placement's rounding; a
micrometre skin kept placed; the placement split over nested instances;
one source compiled unplaced and placed in one batch), the snap's unit
tests (`placed_scale`, the placed reach); probe
`scripts/probe_thin_skin_mutants.py` (#312 mutants).

## Amendment 2026-10-14: a placement's rounding is part of every bound (#295)

An instance's bound was its source's grown by the transform's stretch
(#258), a relative factor. Applying the `f64` transform rounds each
vertex by a few ulps of its coordinates: an absolute error, about `1e-9`
at georeferenced coordinates of `6e6`, that no relative factor covers.
A mesh placed far from the origin could lie further from the exact
surface than its `Proven` or `Certified` bound, by that rounding. The
rule this ADR states for a boolean's guarantee (the mesh against the
exact result of operands perturbed within the reported deviation) needs
the reported deviation to cover how the operands were placed.

- **Decision: every bound carries the rounding of each transform
  applied.** An `Instance` or `InstanceAtStation` reports `s d + r` for
  each source bound `d` (kind kept), `s` the stretch bound, and `r =
  4 eps M`: `M` the largest `sum_j |a_ij| |p_j| + |t_i|` over the
  source's vertices, the magnitude a placed coordinate sums (its own
  where nothing cancels). The constant is derived: a coordinate is four
  terms each rounded at most four times in any order, so it is off by at
  most `gamma_4 S_i < 2.000_001 eps S_i`; over three coordinates `sqrt 3`
  of that, below `3.47 eps M`, with `M` itself computed low by at most
  `1 + gamma_4`. A row that copies a coordinate (one coefficient `+-1`,
  no translation) adds nothing, so an identity or a mirror is exact.
  Nested instances grow the inner `r` by the outer stretch and add their
  own. The term is named under `DeviationPath::Placement`
  (`deviation::PLACEMENT_ROUNDING`).
- **Decision: the source's budget leaves room for it.** The source is
  meshed to the budget less a room, shrunk by the stretch so that
  `s b + room` rounds to within the budget. The room is the floor at the
  translation's scale (`r` at the origin) plus a `2^-20` share of what is
  left, for the source's own coordinates; where the first mesh rounds by
  more and its bound would carry the instance over the budget, the source
  is meshed once more with room for an eighth more than it showed.
- **Decision: a budget at or below the rounding is refused by name.**
  Below the floor at the translation's scale before meshing, or below the
  rounding of the placed vertices after, the compilation is refused with
  `UnsupportedInput` `deviation::BUDGET_BELOW_PLACEMENT_ROUNDING`, in the
  reported and the plain compile paths alike, rather than reported as
  met. A torus at `(6e6, 6e6, 100)` with a `1e-9` budget is refused (the
  floor is about `5.3e-9`).
- **Decision: polylines follow the same rule.** `compile_polyline`
  (#311) through an `Instance` or `InstanceAtStation` scales its chord
  bound by the same stretch and moves its vertices by the same `f64`
  transform, so its bound adds the same term, its source leaves the same
  room, and a budget at or below the rounding is refused by the same name.
- **Rejected: exact-arithmetic bounds that exclude vertex rounding,
  left to consumers.** The rounding of placing is known to the compiler
  and is the dominant error at georeferenced scale; a bound that leaves
  it out is not a bound on the mesh delivered.
- **Unchanged.** The rounding of computing a construction's own vertices
  is not part of any bound; the exact compiler places exactly
  (`ExactBRep::transformed`) and is unaffected.

Tests: `crates/execution/compile/tests/placement_rounding.rs` (the
torus at `(6e6, 6e6, 100)`: a `2e-8` budget met, about a quarter of it
the rounding, and dense samples of the exact torus against the `f64` mesh
moved back exactly within the bound; `1e-9` and the floor refused by
name in both paths; a source far out turned about the origin, meshed
once more; an identity and a mirror adding nothing; a torus at a station
of a far line), `crates/execution/compile/tests/polyline.rs` (a polyline
and a circle placed at the site); the unit tests of `deviation::placement_rounding`; probe
`scripts/probe_placement_rounding_mutants.py` (#295 mutants).

## Amendment 2026-10-16: a cutter face on a curved face of its host (#324)

A quarter-ring wall less a niche whose inner arc lies on the wall's inner
cylinder, or one unit in the last place inside it, was refused ("exact
boolean whose kept faces do not sew into a solid"), so the boolean's
deviation was unbounded. The general boolean already read the two
cylinders as one surface (`same_support`: exact predicate first, else
axis and radius within `eps`, reported as `CoincidentSupports`) and split
each face along the other's boundary. The defect was the side: a region
lying on the other operand's face takes that face's outward normal to
decide whether the operation keeps it, and that normal ignored the
face's winding. An extrusion's face on its profile's inner arc winds
clockwise in its parameters and so faces against its flag, which the
region's own normal accounted for and the other face's did not, so the
difference kept the wall's patch under the niche.

- **Decision: coincident curved faces are decided like coplanar ones.**
  The same reading, split and keep-or-drop rule, with the other face's
  outward normal read through its flag and its winding
  (`split::face_against`).
- **Decision: the reading's move is a deviation contribution.** Where
  faces of the two operands were read as one surface within tolerance,
  the furthest that moved an operand is a `Certified` contribution of the
  boolean path under `deviation::COINCIDENT_OPERANDS`, next to the
  measured one, as the snap's move is (#276). It names the perturbation
  the exact result stands for; it does not bound the mesh's distance to
  the boolean of the unperturbed operands. The 1-ulp niche is within
  tolerance, so it is decided as flush, its move (about `4e-15`)
  reported; the flush niche reports its rounding (`2.2e-15`), since the
  two cylinders are each built from rounded profile points.
- **Unchanged.** A box cutter tangent to the inner face (touching it
  along one ruling) is a different configuration, refused by the general
  boolean as a point too close to a face boundary, and by the mesh path as
  a result touching itself; this amendment does not cover it.

Tests: `crates/algorithms/construction/brep-boolean/tests/coincident_cylinders.rs`
(the three niches, audited, closed, at the closed-form volume, the
reading reported for the two flush ones only),
`crates/execution/compile/tests/curved_flush_niche.rs` (the issue's
reproducer: all three certified at `0.99973 mm` under a `1 mm` request);
probe `scripts/probe_coincident_curved_faces_mutants.py`.

## Amendment 2026-10-18: abutting openings in a curved wall, a box tangent to its face (#332, #327)

Curved standard-case walls with openings were refused by the general
boolean ("a point lies too close to a face boundary to classify"), so the
boolean's deviation was unbounded. The undecided point lay on the side
face of one opening, a few micrometres from the side plane of the
opening next to it: delivered models give an opening in a curved wall
radial sides, and the side of one reaches past its neighbour's by
rounding of its authoring (tens of micrometres in the case reported,
just under a millimetre in another pair). The two planes cross on the
wall's axis, turned by far more than the angular tolerance, so they are
not read as one surface; in the wall they lie within the linear one.
Every point of such a region was refused, and where it was not, the
sections' cuts at the two sides merged at their mean, on neither face's
boundary, and the faces did not sew. The box of #327 was refused the same
way: its back face touches the inner cylinder along a ruling and lies
within the tolerance of it near the ruling.

- **Decision: the start of a classifying ray is read exactly where the
  tolerance decides nothing.** A point within `eps` of a face's support
  is first not trusted, as before. Where no point of a region decides
  so, its points are classified again with a point on a support only
  within the rounding floor, but for the faces read as the point's own
  support, which keep `eps`. The decision is exact for the operands as
  given; the region lies within the tolerance of that support wherever
  its points were taken. The tangent box is decided this way at the
  caller's tolerance.
- **Decision: a boolean refused at the caller's tolerance is cut again
  at `FINER_FACTOR` (`2^-10`) times it.** Features closer than `eps` but
  further apart than the finer tolerance are then apart, decided by the
  exact predicates; the readings still cover what they exist for, the
  rounding of separately placed operands (orders of magnitude below a
  micrometre under a millimetre, and above the `2^-40` floor and the
  `sqrt(f64::EPSILON)` split of a double root). The result is the exact
  boolean of operands moved by at most the finer tolerance;
  `BooleanReport::finer` says it was cut so, and its readings lie within
  it. A boolean that succeeds at the caller's tolerance is never cut
  again, and one refused at both returns the refusal at the caller's.
  This amends "Source of the tolerance" above: the caller's tolerance is
  passed through unchanged, and a fixed fraction of it is the only other
  tolerance any decision uses.
- **Decision: a split's vertices weld only within the tolerance on the
  surface.** Piece ends of a split face weld in parameters within a
  relative slack (`1e-7`), bookkeeping for one vertex's evaluations. On
  an arc of radius `24` it spans `7 um`, more than the finer tolerance,
  and welded two vertices the split had kept apart; they now also have to
  lie within the tolerance in space.
- **Not a new deviation contribution.** Both decisions are exact or
  within the finer tolerance, reported by the boolean as before; the
  measured contribution's detail still says when its result is the
  boolean of operands within tolerance, and a coincidence the retry reads
  is named under `deviation::COINCIDENT_OPERANDS` as before.
- **Changed refusal.** A web hole half a tolerance off a filleted flange
  (#249) was refused as `UnsupportedContact`: read as touching, the
  contact could not join the fillet arcs. At the finer tolerance the gap
  is a gap, and the boolean is cut at the reference volume.
- **Known refusal.** The tangent box's exact result is right, but its
  material is as thin as nothing along the contact, on the tangent plane
  and a unit in the last place past it: the mesh compiler refuses that
  mesh as a result touching itself (#194), so #327's mesh stays
  uncertified there. A unit short of the tangent plane the box touches
  nothing but material and is certified.

Tests: `crates/algorithms/construction/brep-boolean/tests/curved_wall_openings.rs`
(a wall of radii `6` to `6.25` less five wedge openings, abutting by
`15 um` and `0.6 mm`, and one of radii `24` to `24.1` abutting by `7 um`,
each step audited, closed and at the closed-form volume; the tangent box at the closed-form volume, decided at the
caller's tolerance), `crates/execution/compile/tests/curved_wall_openings.rs`
(the wall certified within the millimetre; the box exact, its mesh
refused where it touches itself), `fillet_tangent.rs` and
`exact_tangent_hole.rs` (the half-tolerance hole, cut finer); probe
`scripts/probe_finer_boolean_mutants.py`.

## Amendment 2026-10-20: a faceted face's `0` means planar (#343)

The #257 rule above gave a B-rep face without a surface `0` whenever its
corners lay within the linear tolerance of their fit plane, and the slab
width `w` only beyond it. So `0` meant "within a millimetre", not
"planar": the faceted IFC2X3 walls #343 measured at `(1.5e6, 5.0e6)`
have faces off their plane by up to `4.1e-10` m in exact arithmetic, and
each reported a proven `0`.

- **Decision (maintainer).** A face without a surface whose corners are
  not exactly coplanar reports `Certified(w)` under
  `"non-planar face without a surface"`, the width #261 already computes,
  however small. Coplanarity is decided on the stored corners by the
  certified `orient3d` (a plane through two distinct corners and a third
  off their line, every corner exactly on it), so `0` is now a proof. A
  tilted face whose corners are exactly coplanar still reports `0`.
- **Unchanged.** A face warped beyond the tolerance reports the same `w`
  as before, the triangulation is unchanged, and authored polygon faces
  keep their tolerance rule (#261). A faceted face rounded off its plane
  by a few ulps now reports about 64 ulps of its reach rather than `0`,
  well within any budget, so whether a mesh meets its requested chord
  error does not change for it.

Tests: `crates/execution/compile/tests/brep_warped_faces.rs` (a corner
lifted `9e-4`, `1e-6` and `4e-10` m reports half its lift; an exactly
coplanar tilted top reports `0`), the `planar` unit tests (exact
coplanarity, a diagonal wall on a national grid); probe
`scripts/probe_faceted_planarity_mutants.py`.
