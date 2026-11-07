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
- `crates/algorithms/query/measure/src/exact_distance/witness.rs`: the certified
  foot witnesses (#318).
- `scripts/probe_distance_rounding_mutants.py`: 9 faults in the derived margins,
  all caught (#318).

## Amendment 2026-11-07: rounding margins in ulps of the coordinates, a floor, and foot witnesses (#318)

**Problem.** The margins of the lower bounds were `1e-12` of the points'
magnitudes. At georeferenced coordinates (450000, 5600000) a point cannot be
computed to better than a few `1e-9`, so the search could not close below an
interval about `1e-5` wide, and the Hausdorff search shared the margins.

**Decision.**

- Every margin is derived from `M`, the summed magnitude of the terms a point
  of an element is computed from, taken at its farthest parameters (never
  growing under a split). A side's scale `S` is the largest `M` over its root
  elements, and every element of the side is padded by `K eps S` with the same
  `S`, so a pad cannot reorder the queue: a pair's lower bound is its
  geometric bound less a constant of the search. The constants are `POINT`
  16, `PROJECTION` 8, `RATIO` 8 and `FLOOR` 128, derived in the module docs
  of `exact_distance.rs` (`FLOOR` covers the 48 `eps` of the bounds at zero
  size and the 32 `eps` of the witnesses' evaluation).
- A search settles at the floor `128 eps (S_a + S_b)`, about `3.2e-7` at
  (450000, 5600000) and `6e-14` at unit scale. An accuracy below it is never
  reported as met; the interval returned shows its width. The floor only
  stops refinement, so failing closed is unchanged. The touch margin of #273
  is the same 128 ulps of the points' magnitudes.
- Where a lower bound does not depend on patch size (along a planar patch's
  own normal, against a parallel planar patch or a line edge) the minimum is
  reached on a whole region and only the upper bound never closes. Certified
  foot witnesses close it at once: the foot of a witness on a planar face,
  the closest points of two line spans, and for point queries the foot of the
  point or its closest point on a line span. A foot counts only when its
  parameters lie in a patch already certified inside, or its face's domain
  certifies them.
- Witnesses are counted when their pair or node is queued, and the point
  queries of the Hausdorff search pop equal bounds first in, first out, as
  the key's own documentation said they did.
- The sphere and torus branches of a projected range carried no pad; they do
  now. The Hausdorff matched bound, a translate's re-charted surface and cone
  pole, and the displacement gate use the derived margins, and a plane,
  cylinder or elliptical cylinder is Lipschitz by its axes, not `1 + |v|`
  times as steep.

**Alternatives rejected.**

- Coarsest-first or grid tie-breaks on the queue (tried in a first
  implementation): the upper bound still improves only by patch centres
  lining up. A block over a disc cap ended at `1.19e-9` after the step budget
  and the `meets_early` tests of #273 went to 9 and 12 steps.
- Exact recentring of each B-rep, as for #418 and #419: it leaves
  world-charted parameters untouched (the `u x` terms of caps trimmed in
  world coordinates stay near `5.6e6`, so the derived margins are needed
  anyway), costs a copy of each B-rep per query, and saves a floor of about
  `3e-7` m, far below building tolerance.
- A second phase of the search: with uniform pads it has no ordering
  criterion the first lacks.

**`cut.rs` audit.** The `1e-12` pads at `exact_hausdorff/cut.rs` (L321, L359)
stay. They are sound over-enclosures, and for a translate the coefficient
differences are zero, so they do not widen the bound.

**Follow-ups.** The upper bound excludes its witnesses' evaluation rounding,
about `32 eps S`, so far out `Clearance::Below` can be claimed within about
`1e-8` of a limit (#439). The `1e-9`-relative gates of `plan_patches_overlap`,
`exact_bodies::separate` and `exact_domain::slack` reach millimetres far out
and lead to refusals or no verdict, never to an unsound result (#440).
