# Plan: exact curves (nested roots, conics) and their first users

Status: in progress. Delete this file when every milestone is landed.

## Goal

1. `axiolid-exact` gains values with nested square roots and conics.
2. #155: an exact arc/segment polygon boolean replaces `cavalier_contours`.
3. #120: curved-face exact booleans get their first real slice.

## Milestones (each ends at a green crate test run and a commit)

- M1 `tower.rs`: `Tower<T>` = `a + b*sqrt(r)` with `a`, `b`, `r` towers.
  Radicals carry creation ids; a radicand only contains older radicals, so
  the recursive sign (square with case analysis) terminates. Interval tier
  evaluates numerically with an outward-rounded sqrt (`Arith::sqrt_enclosure`).
  Oracle: fixed-point bigint with nested isqrt; exact zeros from denesting
  identities (sqrt(3+2sqrt2) = 1+sqrt2, sqrt(5+2sqrt6) = sqrt2+sqrt3).
- M2 `poly.rs` + `algebraic.rs`: integer polynomials (num-bigint), primitive
  PRS gcd, Yun squarefree split, Sturm counts, root isolation by bisection
  on dyadic endpoints (split points never roots), `AlgebraicReal` with
  exact compare (gcd test in the overlap) and sign of g(alpha).
- M3 `conic.rs`: `Conic` (six exact coefficients), line/conic hits (reuse
  the quadratic along a line), conic/conic intersection via the resultant
  in y after an integer shear x' = x + k*y chosen so both y^2 coefficients
  are non-zero and no two points share x'. Each point: its x' root,
  y = -s0/s1, multiplicity from the squarefree split; crossing iff odd.
  Common component refused by name.
- M4 #155 in `axiolid-overlay`: edges = segments and arcs (bulge is exact:
  centre = (2b(P+Q) + (1-b^2) perp(Q-P)) / 4b). Intersections exact
  (points in Q(sqrt D)), split edges, classify each sub-edge inside/outside/
  shared against the other operand by the germ order at on-boundary
  vertices (tangent then curvature), else one ray cast. Chain, orient by
  the leftmost point, nest holes by ray casting (nested roots appear here).
  No tolerance in any decision; output coordinates rounded once.
  `cavalier_contours` dependency removed; ADR 0069 superseded in part.
- M5 #120: whatever the exact planar layer unlocks in `boolean_exact`
  (arc prism results with holes, disconnected parts) and the next curved
  family, scoped after M4 with the ledger rows C9/D3 named precisely.
- M6 ledger, ADR, changelogs, issues closed/updated, gate, push, CI.

## Invariants

- Every topological decision is an exact sign (filter, then exact).
- Nothing reads a tolerance except input validation.
- Each milestone: tests with an independent oracle + mutation probe entries.
