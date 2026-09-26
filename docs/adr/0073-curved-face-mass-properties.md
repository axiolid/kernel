# 0073 — Exact mass properties over curved faces by Green's theorem

- **Status:** Accepted
- **Date:** 2026-09-26
- **Deciders:** Friedrich Schrödter
- **Supersedes:** —

Part of axiolid/kernel#125 (ledger row C17).

## Context

`exact_properties` measured an `ExactBRep` without tessellating it, but only
when every face was planar. A cylinder, cone, sphere, torus, elliptical
cylinder or B-spline face was refused by name, and a planar face bounded by an
arc could not be measured either: the planar path fans the loop's *vertices*,
which drops the segment between a chord and its arc. Every curved construction
since v0.6 (revolutions, arc prisms, column booleans, sloped clips) was
therefore checked against hand-rolled cap sums or sampled volumes.

OCCT's `BRepGProp` integrates each face over its parameter domain with Gauss
quadrature, walking the domain's boundary (the pcurves) and integrating the
other parameter from a reference line. Every exact face here already carries
pcurves with oriented intervals (ADR 0024), so the same reduction applies
without inventing any data.

## Decision

We integrate every non-polygonal face by **Green's theorem over its own
parameter domain**, with the same cone-from-origin fields the planar fan uses,
so the two paths sum consistently on one solid.

- The face contributes `int int_D g du dv` with
  `g = (S . N) * [1/3, S/4, S^2/5]` and `N = S_u x S_v`, plus `|N|` for area.
  Green's theorem turns it into `- oint H(u, v) du` round the pcurves, with
  `H` the integral of `g` in `v` from a reference line.
- Both integrals are adaptive Gauss-Kronrod (G7/K15), accepted per panel at a
  relative error of `1e-13` of `int |g|`, with a noise floor of `1e-15 L^d`
  for components that vanish by symmetry. A face that does not converge
  within 2048 panels is refused (`NotConverged`), never returned approximate.
- **Seams.** On a surface periodic in `u`, a loop may end whole periods from
  where it began (a cylinder wall bounded by two circles). Along a seam `u` is
  constant, so it contributes nothing to `- H du`; open loops are integrated
  as they are, and junctions are unwrapped by whole periods.
- **Poles.** If the loops' net winding in `u` is non-zero, the domain reaches
  a pole (sphere) or apex (cone) on the side the winding keeps on its left.
  The reference line is placed on that pole, where the missing boundary
  contributes exactly zero. A net winding with no pole on that side (a lone
  circle on a cylinder) bounds nothing and is refused (`ParameterDomain`).
- **Tube winding.** On a torus face whose loops wind in `v` (bounded by
  meridians), the roles swap: `+ oint G dv` with `G` integrated in `u`.
- **Orientation** is the loops' winding about `N`, after the face, shell-use
  and bound flips — the sense the planar fan, the tessellator and edge-use
  pairing already read. A loop anticlockwise in `(u, v)` runs anticlockwise
  about `N` whatever the frame's handedness, so no declared normal is
  consulted.
- A planar face whose edges are all straight lines keeps the vertex fan. Its
  area is now the length of the summed *vector* area, so a hole subtracts.

### The revolution frames were left-handed

Measuring revolutions exactly exposed that every exact revolution was built
inside out: its surface frames were `(x, y, z) = (X, Z, Y)`, which is
left-handed, and every loop is built anticlockwise in its parameters, so each
`Forward` face pointed into the solid. The topological and geometric audits
compare faces with each other and passed it; `exact_properties` reported
`-2 pi R A` for every Pappus fixture. The frames are now `(X, -Z, Y)`,
right-handed, which flips every face together and changes nothing else
(`revolve_contour::frame_at`).

## Alternatives considered

| Option | Why not |
| --- | --- |
| Closed forms per surface family and pcurve type | Exact to rounding, but a separate derivation for every (surface, pcurve) pair: a sinusoid rim on a cylinder, a circle rim on a B-spline disc, and each future pcurve type. K15 already integrates every elementary integrand exactly on one panel (a low-order trigonometric polynomial); closed forms buy nothing measurable. |
| Keep refusing curved faces | Leaves C17 narrow and every curved construction unmeasured, which is how the revolution orientation bug went unseen. |
| Integrate the declared surface normal, with loop direction read from the domain's signed parameter area | Makes every face measure positive regardless of how it is wound, which hides exactly the defect found here. The winding is the orientation data; it must be read, not normalised. |
| Tessellate and measure the mesh | The error `exact_properties` exists to avoid. |

## Consequences

**Positive**

- Every exact solid the kernel builds is now measurable exactly, and every
  curved constructor test can use a signed closed-form oracle.
- Planar faces bounded by arcs or ellipses measure correctly.
- The revolution orientation defect is fixed and pinned by a geometric test
  that does not depend on this module.

**Negative / costs**

- `ExactMeasureError` gains `ParameterDomain`, `Evaluation` and
  `NotConverged` and becomes `#[non_exhaustive]`: a breaking change for
  `axiolid-measure` 0.x.
- `axiolid-measure[exact]` now depends on `axiolid-evaluate` and
  `axiolid-curve`.
- The result is quadrature to a stated bound, not symbolic. For elementary
  faces the bound is met on the first panel; B-spline faces with many knots
  subdivide.

**Follow-ups / risks to watch**

- A torus face bounded in both directions without seams (a whole torus with
  no edges) is refused by name.
- A loop that runs along a pole by half a turn or more relies on the
  pole-jump branch in `join`, which no fixture exercises yet.
- The mesh `revolve` of a profile on the negative side of its axis also
  produces an inside-out mesh; it is a separate path and is not changed here.

## Relation to existing code

- `crates/algorithms/query/measure/src/exact.rs`: routing, orientation, errors.
- `crates/algorithms/query/measure/src/exact_face.rs`: the face integral,
  seams, poles, quadrature, and hand-built pole/tube/B-spline fixtures.
- `crates/algorithms/construction/construct/src/revolve_contour.rs`
  (`frame_at`) and `revolve_exact.rs`: right-handed revolution frames.
- `crates/algorithms/construction/construct/tests/curved_measure.rs`:
  closed-form oracles per family.
- `scripts/probe_curved_measure_mutants.py`: 18 faults, all caught.
