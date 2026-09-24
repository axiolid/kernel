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

Status: step 1 in progress.
