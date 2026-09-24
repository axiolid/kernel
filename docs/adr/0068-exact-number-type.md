# 0068 — Exact numbers: num-bigint integers under an owned filtered layer

- **Status:** Accepted
- **Date:** 2026-09-23
- **Deciders:** Friedrich, axiolid
- **Relates to:** [0016](0016-predicate-ownership-and-adopted-implementations.md), [0050](0050-arc-aware-planar-overlay.md), [0069](0069-arc-boolean-path.md); issues #154, #120, #155

## Context

Axiolid's exact arithmetic today is predicate-only: Shewchuk-style
expansions give proven signs for `orient2d`, `orient3d`, `incircle` and
`insphere`. Nothing can hold an exact *constructed* value, such as the
point where two segments cross, or where an arc meets a line. Curved
B-rep booleans (#120) and an exact arc overlay (#155) both need that.

Two questions: which big-number crate, and what to build over it.

## Measured

`docs/research/exact-arithmetic-bench/` (standalone package, not a
workspace member). Workload: intersect two segments exactly, then take the
sign of `orient2d(a, b, intersection)`. 20,000 random cases, full 53-bit
coordinates in [-1000, 1000), identical inputs everywhere. Median of 7
passes, two runs, Xeon w7-3565X, release build.

| Form | Crate | Per case |
| --- | --- | ---: |
| Normalised rationals | `num-rational` 0.4.2 | 76 us |
| Normalised rationals | `dashu-ratio` 0.6 | 6.6 us |
| Integers, division-free | `num-bigint` 0.4 | 0.84 us |
| Integers, division-free | `num-bigint` 0.5.1 | 0.85 us |
| Integers, division-free | `dashu-int` 0.6 | 0.58 us |

All five agree on all 20,000 signs.

Two readings matter more than which crate won:

1. **Formulation beats library by 10-90x.** A rational type pays for a gcd
   at every step. Clearing denominators once (the point is `(Nx/D, Ny/D)`,
   so `sign(orient) = sign(D) * sign(orient * D)`) leaves integer products
   only. The same crate got 90x faster by changing the formulation, not
   the library.
2. **Between integer crates the gap is 1.5x.** `dashu-int` is faster;
   `num-bigint` is within 50%.

## Decision

1. **Big integers: `num-bigint` 0.5**, behind our own types so the choice
   can be swapped. MIT/Apache, pure Rust, three crates total (`num-traits`
   already in our lock), the ecosystem default (about 600M downloads).
2. **No general rational type.** Constructions return homogeneous
   integer coordinates (numerators over a shared denominator); a sign is
   taken by clearing denominators, never by dividing.
3. **Filter first.** Every exact operation is guarded by an f64
   interval/error-bound evaluation, and escalates to big integers only
   when the filter cannot decide. That is how the predicates already stay
   cheap, and it removes the big-integer cost from almost every call.
4. **Square-root extension for arcs.** A line meets a circle at
   `a + b*sqrt(c)`. We own a small type for that (sign by squaring with
   case analysis, as CGAL's `Root_of_2`), built on the integers above.
   Without it, the number type does not reach the curved cases that
   motivated it.

## Alternatives considered

| Option | Why not |
| --- | --- |
| `dashu-int` 0.6 | 1.5x faster here, but 6 crates, ~2M downloads, one main maintainer. The filter (point 3) means the big-integer path rarely runs, so its speed matters less than its longevity. Revisit if profiles of real workloads show the exact path hot. |
| `num-rational` 0.4.2 | 76 us per case, 90x slower than integers. Also still requires `num-bigint` 0.4, so adopting it pins an old major. |
| `dashu-ratio` 0.6 | 8x slower than its own integers; the gcd cost is inherent to normalised rationals, not the crate. |
| `rug` (GMP) | Fastest at large sizes, but a C dependency: breaks pure Rust, cross-compilation and wasm. LGPL-3.0+. |
| `malachite` | LGPL-3.0-only and MSRV 1.90 (ours is 1.88). |
| Build our own big integers | Solved problem; nothing Axiolid-specific to gain. |

## Consequences

**Positive**

- #120 and #155 share one exact substrate instead of inventing two.
- The adopted part is the commodity part (big-integer multiply); the
  geometric parts (filters, homogeneous constructions, square-root
  extension) are ours, testable and documented.

**Negative / costs**

- The square-root extension is real work: sign determination for
  `a + b*sqrt(c)` compared against another such number needs care.
- Homogeneous coordinates leak into construction APIs; results must be
  rounded to f64 at a boundary we choose explicitly.

**Follow-ups / risks to watch**

- Degree growth: nested constructions multiply bit lengths. Each exact
  construction documents its maximum degree.
- The benchmark is one workload. Re-measure on #120's real cases before
  claiming anything about boolean performance.
- `num-bigint` sits behind our types, so switching to `dashu-int` later
  is a local change.

## Relation to existing code

The predicates in `axiolid-reference` stay as they are. The new layer
reuses their filter pattern; it does not replace them.
