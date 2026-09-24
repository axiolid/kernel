# PLAN: exact curves, then users (#155, #120)

Standing goal (Friedrich, 2026-09-24): finish axiolid-exact's "not done"
(nested square roots, conics), then #155 and #120 so the crate has users.

## Done (local commits, see git log)
- M1 Tower: nested radicals, depth <= 6, numeric filter via sqrt_enclosure.
- M2 IntPoly/RealRoot: Sturm isolation, exact compare, sign_of(q) at a root.
- M3 Conic: line/conic (Root2), conic/conic (shear + resultant + RealRoot),
  exact tangency, side_of_line, sign_of_conic, exact x()/y().

## M4 #155: exact arc overlay (replaces cavalier_contours)
Input stays `ArcRing` (f64 point + f64 bulge). Key fact: an arc through
f64 endpoints with f64 bulge b has a circle with EXACT dyadic implicit
coefficients after multiplying by 4b (centre = m + perp(d)(1-b^2)/(4b)).
So every edge is exactly a line or a circle; no rounding before output.

Algorithm (O(n^2) pairs; building profiles have tens of edges):
1. Edges -> exact supporting curve (Line / Conic circle) + span.
2. Pairwise intersections exactly (line/line dyadic-rational,
   line/circle Root2, circle/circle: radical line -> line/circle, so
   Root2 too). Keep points on both spans (exact span tests).
3. Split edges; order split points along each edge exactly; merge
   coincident points by exact equality (vertex identity).
4. Overlapping collinear/co-circular pieces: detect exactly (same
   supporting curve), keep one copy with both operands' insideness.
5. Classify each sub-edge by the other operand's region: status toggles
   at transversal crossings from an exactly-classified start vertex;
   tangencies do not toggle.
6. Select sub-edges by operation, assemble loops (leftmost turn at each
   vertex, angle order decided exactly), outer vs hole by exact signed
   area sign / nesting.
7. Output f64 ArcRing: vertices rounded once, bulges recomputed from the
   rounded endpoints and the exact circle. Claim: TOPOLOGY exact
   (existence/order/coincidence of vertices, in/out), coordinates
   rounded once at output.
Then drop cavalier_contours; ADR 0069 addendum; F1 -> implemented if
holes/multi-region/degenerate overlaps all covered.

## M5 #120 (scoped honestly)
General curved B-rep booleans (OCCT BOPAlgo scale) do not fit one
landing. Deliverable: coaxial arc-prism booleans on the exact overlay,
lifting today's refusals (result with holes, several components, base
not at z = 0). C9 stays NARROW with the refused subset named
(non-coaxial curved solids). #120 stays open.

## Risks
- Output rounding can make a result ring self-touch at 1 ulp: re-validate
  output with validate_arc_ring; refuse rather than return invalid.
- Circle/circle via radical line needs the circles non-concentric;
  concentric = no intersection or identical (overlap case 4).
