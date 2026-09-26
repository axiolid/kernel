# 0074 — Certified boundary distance by branch and bound

- **Status:** Accepted
- **Date:** 2026-09-26
- **Deciders:** Friedrich Schrödter
- **Supersedes:** —

Part of axiolid/kernel#125 (ledger row C18).

## Context

Distance between shapes existed only as floating-point closest-point
witnesses between primitives (`ClosestPoints3`) and meshes. A rule checker
that compares a clearance with a limit needs more than a value: a
circular column 1.2 mm from a wall that must be 1.2 mm away cannot be
passed or failed on a number within rounding of the limit. It needs an
interval certain to contain the true distance, or an explicit
"indeterminate".

OCCT's `BRepExtrema_DistShapeShape` finds extrema numerically; its answer is
a value with a tolerance, not a certificate.

## Decision

`axiolid-measure` (feature `exact`) computes the distance between the
**boundaries** of two exact B-reps as a certified interval, by branch and
bound over pairs of elements: face patches (rectangles of a face's
parameters) and edge spans.

- **Lower bound.** Each element has an enclosing sphere from Lipschitz
  bounds on the exact surface or curve, and an exact range of `d . x` along
  a direction `d` for planes, cylinders, elliptical cylinders, cones,
  spheres, tori and line/circle/ellipse edges. `d` is the line between
  centres and each patch's normal, which makes the bound second order where
  surfaces face each other. Rounding margins are added everywhere.
- **Upper bound.** Only points on the boundaries: edge points, and surface
  points at parameters a face's domain is certified to contain. The domain
  classifier splits each pcurve where either parameter turns, so every
  piece is monotone and a ray crossing is decided from its ends; it answers
  "don't know" near the boundary rather than guess, and understands seams,
  poles and tube-wound torus faces (shared with ADR 0073's integrator).
- **Pruning.** A patch certified outside its face is dropped. A face patch
  whose normals cannot point at the other element is dropped from that
  pair: a closest pair of separated boundaries is critical on each face it
  lies inside or lies on an edge, and edges are elements of their own. A
  cone patch reaching its apex is never dropped (the apex is not smooth),
  nor is a patch of a face with an edge the query cannot bound.
- **API.** `boundary_distance(a, b, accuracy, tolerance)` refines to an
  accuracy; `boundary_clearance(a, b, limit, tolerance)` refines only until
  the interval clears the limit and otherwise returns
  `Clearance::Indeterminate`.

## Alternatives considered

| Option | Why not |
| --- | --- |
| Tessellate and take the mesh distance with the chord error as a margin | The chord error is a target of the tessellator, not a certified Hausdorff bound. |
| Numerical extrema (Newton on the distance function) | Fast and precise, but a local method: it certifies nothing about the global minimum. |
| Sphere bounds only | First order: an isolated nearest pair needs `1/eps` refinement per dimension; the tests did not close at `1e-9` within a 400k-step budget. |
| Scoring every possible split by the bound it yields | Tried: bounds do not rise monotonically under refinement, and the scored choice stalled on cases the simple longer-side rule closes. |

## Consequences

**Positive**

- Rule checks get a certified verdict or an explicit indeterminate.
- Isolated nearest pairs (pole to pole, apex to sphere, wall to block
  face, block over a disc cap) close to `1e-9` in well under a second.

**Negative / costs**

- Where the nearest points form a whole line (two parallel columns), every
  slice along it is a near-minimal pair; refinement slows to about `1e-6`
  within the step budget. The interval stays sound.
- It is the distance between boundaries: a solid wholly inside another
  measures the gap between their surfaces, not zero. Containment is a
  separate classification.
- B-spline faces are bounded by their whole control net, which is sound but
  never refined; edges on B-spline or intrinsic curves are not elements, so
  their faces are never pruned by the normal test.

**Follow-ups / risks to watch**

- A second-order bound along lines of contact (e.g. splitting across the
  direction the normal turns) would close parallel columns faster, but it
  must not reintroduce the stalls the scoring experiment hit.

## Relation to existing code

- `crates/algorithms/query/measure/src/exact_distance.rs`: bounds, pruning,
  search, API.
- `crates/algorithms/query/measure/src/exact_domain.rs`: pcurve domain
  assembly shared with `exact_face.rs`, and the certified classifier.
- `crates/algorithms/construction/construct/tests/boundary_distance.rs`:
  closed-form fixtures.
- `scripts/probe_boundary_distance_mutants.py`: 11 faults, all caught.
