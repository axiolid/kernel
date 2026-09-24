# Plan: arc speedup, #119, publish axiolid-exact, #120 (round 2)

User (2026-09-24): speed up arc booleans, then #119, then #120; publish
axiolid-exact before or after, whichever makes sense.

## S. Arc overlay broad phase — done, 4e9e2f9
Padded f64 edge boxes (chord box + sagitta), box-filtered crossing and
shared-edge tests, sorted start index for linking. 256x256 wavy rings
779 ms -> 7 ms. 12/12 mutants.

## #119 tranche 1: exact analytic intersections (this session)
Full #119 (NURBS x NURBS marching, tangent/overlap in general) is OCCT
IntPatch scale (~32k lines). Tranche 1 is the analytic core every curved
boolean needs, built on axiolid-exact so decisions are exact:

1. Curve/surface (B9): Line3, Circle3, Ellipse3 x Plane, Cylinder, Cone,
   Sphere, Torus. Substitute the curve's polynomial / half-angle rational
   parametrization into the surface's implicit equation, written for the
   GIVEN f64 frame (non-unit axes kept as |z|^2 factors, not normalised),
   so coefficients are exact dyadics; isolate real roots with IntPoly /
   RealRoot. Report per hit: parameter (exact root + f64), point (f64),
   multiplicity (tangent iff p' also vanishes). Identically zero poly =>
   curve lies in the surface (Contained). Half-angle misses theta = pi:
   test that point exactly. Cone slope = tan(semi_angle) rounded once
   (documented: the cone is the one with that f64 slope).
2. Curve/curve 2D (B8): Line2/Circle2/Ellipse2 pairs via axiolid-exact
   conic (M3) — gives conic.rs its first user.
3. Surface/surface (B10): exact_surface_intersection's degeneracy tests
   (parallel, tangent, coincident: raw `== 0.0` on rounded f64) become
   exact sign decisions. Unrepresentable quartics stay refused by name.
Ledger: B8/B9/B10 stay narrow (NURBS generality open) with the analytic
subset named; #119 stays open for NURBS.

## Publish axiolid-exact 0.1.0 — after #119 tranche 1
#119 is the last planned change likely to add exact API (poly
arithmetic helpers); publishing after it avoids a 0.1 -> 0.2 break days
later. Overlay/construct next versions need it on crates.io anyway.

## #120 tranche 2
Candidates (pick by value/effort after #119):
- arc prism cut by a plane (column under a sloped roof): cylinder x plane
  = ellipse, vertical edges x plane = points (#119 pieces).
- disconnected coaxial results returned as several solids.
- stepped coaxial spans (see boolean_stepped.rs for the polygon path).

## Progress

- Broad phase: done, 4e9e2f9 (gate + CI green). 2-116x on the scaling bench.
- #119 tranche 1 (exact analytic curve/surface and curve/curve): done
  locally; 12 tests, 6/6 mutants caught (scripts/probe_exact_curve_mutants.py).
  B8/B9 ledger rows extended, still narrow (B-spline operands).
- #119 tranche 2: done, f99e346 (gate green). Sphere/plane tangency,
  cylinder/plane parallel/perpendicular/tangent and cone/plane
  perpendicularity decided exactly; 37 tests, 5/5 mutants caught
  (scripts/probe_exact_surface_mutants.py). Three of the old float
  decisions were proven wrong on exact inputs, each now a regression test.
- axiolid-exact 0.1.0: published to crates.io from f99e346, checksum
  8584a50f..., tag axiolid-exact-v0.1.0; a fresh registry consumer builds.
- #119 stays open: B-spline operands, general quadric/quadric curves.
- Next: #120 tranche 2.


## #120 A: arc prism cut by a sloped plane (option 1, user-approved 2026-09-24)

Problem: a sloped plane cuts a cylinder wall in an ellipse. Its pcurve on the
cylinder is v(u) = mean + a cos u + b sin u, which no Curve2 variant held
exactly (B-spline only approximates, u is an angle).

Steps (each gated before the next):
1. axiolid-curve: `Curve2::Sinusoid(Sinusoid2 { mean, cosine, sine })`,
   point(t) = (t, mean + cosine cos t + sine sin t). Additive
   (#[non_exhaustive]). ADR 0071.
2. axiolid-evaluate: domain2 (full turn), evaluate2, derivative2,
   second_derivative2, closed-form inversion (t = p.x). Graph validation
   accepts finite coefficients.
3. construct: arc-prism builder takes lower/upper levels (Flat | Sloped
   plane). Straight walls: exact Line3/Line2. Curved walls: Ellipse3 in its
   principal frame (t = cylinder angle - const), Sinusoid2 pcurve on the
   cylinder, Ellipse2 pcurve on the sloped cap. Flat/flat output unchanged.
4. Public: clip an ArcPrism by a HalfSpace. Plane must clear both caps over
   the whole section (arc interior extremes included); plane above/below the
   whole prism -> unchanged / Degenerate; crossing a cap or vertical plane ->
   refused by name.
5. Tests: B-rep volume from sampled cap faces (walls vertical, so
   V = sum of int z n_z dA over the caps) against the closed form; geometric
   audit (pcurves lifted onto 3D curves); refusals; mutation probe.
6. Ledger C9 + B1, changelogs, gate, push, CI, #120 comment.

Status: done, 1f4f1b5 (CI green). 14 clip tests, 11/11 mutants.


## #120 remaining (user: "do the remaining 120 work", 2026-09-24)

Issue scope: close C9/D3 to implemented, or to narrow with the refused subset
named. Remaining refusals: stepped spans, planes crossing a cap, non-coaxial
curved solids, cones/spheres/tori/NURBS.

Key observation: every remaining VERTICAL-column case (stepped coaxial
booleans, a plane crossing a cap, a stepped result under a roof) is one shape:
a planar arrangement of cells, each cell carrying a stack of z-intervals
bounded by planes (flat or sloped). Cell data depends only on which operands
contain the cell, so it is a function of a membership bitmask. One builder
for that shape replaces two diverging special cases.

Steps (each gated before the next):
T1. axiolid-overlay: exact N-operand arc arrangement (`arc_arrangement`):
    split all operand edges at all crossings (exact), label each piece's left
    and right membership masks, dedupe shared pieces, snap f64 vertices;
    `region(pred)` links the pieces bounding {mask : pred(mask)} into
    outer/hole rings of piece uses. Refactor `exact_arc::boolean` to share
    the ring assembly. One vertex table => every face built from it agrees on
    every vertex bit-for-bit (the reason not to call arc_overlay per band).
T2. construct: column-solid builder. Input: arrangement, planes, and
    mask -> intervals [lo plane, hi plane]. Walls per piece per maximal
    symmetric-difference span; caps per (plane, facing) via region();
    vertical edges split at every height met at a vertex (heights snapped
    within tolerance, zero-length sides dropped -> triangle walls); connected
    components -> one ExactBRep each; void shells attached to their outer.
    Refuse by name: bands touching only along an edge (non-manifold).
T3. Stepped coaxial booleans (arc and polygon prisms) through T2: the
    `_solids` variants and single-solid variants stop refusing stepped spans.
    Provenance names: walls from (operand, ring, edge), caps from the operand
    whose bottom/top lies at that level.
T4. clip_arc_prism_exact: a plane crossing a cap through T2 (cells split by
    the plane's intersection lines with the top/bottom levels).
T5. Non-coaxial curved / cones / spheres / tori / NURBS: general
    surface-surface B-rep boolean (OCCT BOPAlgo scale). Not attempted here;
    ledger rows stay narrow with this subset named, per the issue's scope
    rule. Propose a follow-up issue.
T6. Ledger C9/D3, changelogs, ADR 0072, gate, push, CI, #120 comment.

Status: T1-T4 done (484c96d arrangement, ab708eb measure fix, column
commit on top). T5 not attempted: non-coaxial curved / cones / spheres /
tori need a general surface-surface B-rep boolean; C9 stays narrow with that
subset named; follow-up issue proposed to the user, not filed. Cavities are
refused until tessellation/measure read void shells.

### #120 tranche 3 design notes (column builder)

- Result solid = union of "column cells": plan region R_k x height interval
  [lo_k(p), hi_k(p)] where lo/hi are flat or sloped planes (z = a + gx x + gy y).
- Input: arrangement of all section rings (ArcArrangement), and per arrangement
  face a list of disjoint z-intervals (bottom plane, top plane). Adjacent
  faces with identical interval lists merge (arrangement.regions predicate).
- Faces emitted:
  * caps: for each distinct (plane, side) group, regions of the arrangement
    where that plane is an interval end -> planar face with arc/line pcurves
    (sloped plane: Ellipse2 in plane frame as in sloped_cap_loop).
  * walls: for each arrangement edge, the two adjacent cells' interval lists
    differ -> the vertical strip set difference (symmetric in z) along the
    edge is a wall. Each wall strip on one side between two planes becomes
    a face on the carrier (plane for line, cylinder for arc) bounded by
    bottom/top curves (line/circle/ellipse) and vertical lines.
- Vertical edges at arrangement vertices must be split at every height
  where any incident wall or cap boundary meets that vertex -> build per
  vertex a sorted list of distinct z values and emit vertical edges between
  consecutive ones, shared by all walls at that vertex.
- Mesh compiler tessellates only solids()[0].outer: build ONE shell per
  connected component; results with an enclosed void are refused by name
  (a stepped column with an internal cavity cannot occur for coaxial
  prism booleans of two operands anyway -- verify).
- Planes compared exactly? Levels come from input heights / half-space
  planes directly (no derived values), so equality is by value equality of
  the defining coefficients -> merge caps only when coefficients are equal.

### T2 implementation decisions (column.rs)
- Two phases: (A) abstract faces over edge keys Rim(piece, class rep plane) /
  Vert(arr vertex, height cluster k); (B) union-find faces by edge keys ->
  shells; each key must be used exactly twice (>2 = touching along an edge,
  refused by name); emit one ExactBRep per outer shell.
- Height classes per piece: planes equal (tol) at start/mid/end of the piece.
  Vertex heights: cluster rim endpoint heights per arrangement vertex (tol).
- Walls: per piece, symmetric difference of the two side stacks in class
  gaps; maximal same-side runs = one face. Built on the DIRECTED piece with
  solid on the left, always face Forward (the tested prism convention);
  rim uses Reversed when the directed piece runs against the piece.
- Caps: ArcArrangement::regions(pred: some block ends on plane P from that
  side); up = Forward, down = Reversed (build_arc_rings convention).
- Outer vs void shell: signed volume from caps only (walls vertical =>
  n_z = 0): sum over caps of +-(h A + gx Mx + gy My), Green closed forms for
  segments and arcs. Voids attach to the single outer shell; several outer
  shells plus a void is refused (unreachable for 2-operand prism booleans:
  a void needs a difference whose tool is enclosed, so one component).
- Mesh compiler must tessellate void shells too (currently outer only).
