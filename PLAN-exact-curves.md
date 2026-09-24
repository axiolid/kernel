# Exact curves: nested radicals, conics, #155, #120

Standing goal (user, 2026-09): take on the "not done" list of `axiolid-exact`
(values with more than one nested square root, conics), then #155 (exact arc
booleans) and #120 (curved B-rep booleans) so the crate has real users.

## Milestones

| # | What | State |
|---|------|-------|
| M1 | `Tower`: nested square roots, any depth (capped at 6) | done, a6636bb |
| M2 | `IntPoly` + Sturm + `RealRoot`: exact real roots of integer polys | done, 9659680 |
| M3 | `Conic`: exact line/conic and conic/conic points | done, 17961b2 |
| M4 | #155: exact arc overlay replaces `cavalier_contours` | done, ADR 0070 |
| M5 | #120: coaxial curved booleans with holes and raised bases | done (C9 stays narrow; #120 open on #119) |
| M6 | ledger rows, ADR, changelogs, close issues, gate, push, CI | pending |

M1-M3 pushed (a85fd09), gate + CI green.

## M4 design: exact arc overlay (#155)

Public API unchanged: `arc_overlay(subject, clip, op, tolerance)`.
`tolerance` only feeds `validate_arc_ring`; no decision uses it.

### Arc edge, exactly
Edge P0 -> P1 (f64), bulge b, chord D = P1 - P0, perp(x, y) = (-y, x).
- Circle (conic form, all coefficients exact dyadic): with k = 4b,
  kC = k*M + (1 - b^2) * perp(D)  (M = chord midpoint),
  k^2 |X|^2 - 2 X.(k*kC) + |kC|^2 - |D|^2 (1 + b^2)^2 = 0.
- On-arc test for X on the circle: X is on the arc iff
  orient(P0, P1, X) == -sign(b), or X is P0 or P1. (The chord line splits
  the circle into two arcs; b > 0 runs CCW through the right side.)
- Rational parametrization by s in [0, b] (s = b at P0, s = 0 at P1):
  w = b (1 + s^2)^2, g = (b - s)(1 + b s),
  X(s) = P0 + (g / w) * ((1 - s^2) D - 2 s perp(D)).
  Dyadic s gives an exact dyadic homogeneous point on the arc. This is the
  half-angle substitution s = tan(phi/2); the bulge is the parameter bound.
- Order along the arc: A before B iff sign(b) * orient(P0, A, B) > 0.

### Exact points
x = (ax + bx sqrt(d)) / w, y = (ay + by sqrt(d)) / w, all dyadic, w != 0.
Covers f64 vertices, segment/segment (rational), line/circle and
circle/circle (one radical per point). Predicates on up to three points
embed them in one `Tower` (depth <= 3) and run filter-then-exact.

### Algorithm
1. Validate; orient both operands CCW.
2. Split points per edge: all subject x clip edge pairs (closed edges, so
   touching vertices count). Circle/circle via the radical line.
   Co-circular and collinear overlaps: split at endpoints lying on the other.
3. Sort split points along each edge exactly, dedup by exact equality.
4. Sample each sub-edge at an exact dyadic interior point (bisection on the
   edge parameter, exact comparisons; seeded by f64).
5. Classify the sample against the other operand:
   - on its boundary: shared edge; compare tangent directions exactly.
   - else ray parity with a dyadic ray direction, re-chosen whenever the ray
     hits a vertex or is tangent to a circle (finitely many bad directions).
6. Select (A = subject, B = clip):
   - union: A out of B, B out of A, shared same-direction once.
   - intersection: A in B, B in A, shared same-direction once.
   - difference: A out of B, B in A reversed, shared opposite once.
   - xor: all non-shared, the "in" ones reversed; shared dropped.
7. Link at exact vertices; with several choices prefer the same operand's
   next edge (touching configurations), else refuse rather than guess.
8. Nest by containment depth (even = outer, odd = hole), exact tests.
9. Output: vertices rounded to f64 once, at the end; sub-arc bulges from the
   original circle and the rounded endpoints.

### Honesty boundary
Topology (crossings, order, in/out, linking, nesting) is exact for the given
f64 input. Output coordinates of constructed points are rounded to f64 once.
Self-intersecting input rings stay outside the contract, as before.

## M5 scope (#120)
General exact curved B-rep booleans are out of reach here; #120 lands as
`narrow` with the refused subset named. In scope: coaxial prismatic solids
with arc sections (cylinders, rounded walls, round openings), including
results with holes (currently refused), on the M4 overlay.
