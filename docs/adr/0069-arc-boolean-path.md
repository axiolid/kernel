# 0069 — Arc booleans: keep cavalier_contours, scaled; exact path later

- **Status:** Accepted
- **Date:** 2026-09-23
- **Deciders:** Friedrich, axiolid
- **Relates to:** [0047](0047-absorb-mesh-boolean.md), [0050](0050-arc-aware-planar-overlay.md), [0068](0068-exact-number-type.md); issue #155

## Context

ADR 0050 adopted `cavalier_contours` 0.9 for planar booleans with circular
arcs, behind `arc_overlay`. It works to machine precision on the probed
cases, but it is tolerance-based, not exact: it has roughly 545 fuzzy
comparisons, with fixed thresholds in drawing units (1e-5 for positions,
1e-8 for its default fuzzy equality, 1e-5 and 1e-3 in its debug
self-checks).

Fixed thresholds mean the answer depended on units. Measured before the
fix, with each drawing using its own 1 um tolerance:

| Case (drawn in mm, then in m) | mm | m |
| --- | --- | --- |
| Union of squares 5 um apart | gap kept | gap lost |
| Intersection overlapping by 5 um | 0.005 mm^2 | 0 |
| Disc crossing a wall edge by 4 um | closed form | 0.004 mm^2 short |

The row F1 (exact polygon booleans) stays narrow while arcs go through
this path. The question was whether to keep, absorb, or replace it.

## Decision

1. **Now: keep the crate, and make it unit-independent.** `arc_overlay`
   scales the drawing by a power of two so the backend's native 1e-5 lands
   on the caller's linear tolerance, then scales the result back. All of
   the backend's thresholds move together, which is why this works where
   overriding its `pos_equal_eps` option alone does not: that trips its
   own debug consistency checks. The scale is capped so coordinates stay
   below 1e6, where f64 still resolves its finest (1e-8) comparisons.
2. **Later: an exact arc overlay of our own**, on ADR 0068's integers
   and square-root extension. #155 is blocked by #154 for that reason.
3. **Not absorbing.** We do not vendor `cavalier_contours` into the tree.

## Alternatives considered

| Option | Why not |
| --- | --- |
| Absorb, as boolmesh (ADR 0047) | boolmesh was absorbed to keep and own its algorithm. Here the part we would own is the ~13k lines of tolerance-based logic we intend to replace with exact arithmetic. Absorbing buys maintenance of code we do not want. |
| Set `pos_equal_eps` to our tolerance | Tried. The backend's other thresholds stay fixed, and its debug assertions fire (`EndPointOnFinalOffsetVertex`) on the 5 um cases. |
| Leave it as it was | Answers depended on drawing units; three tests prove it. |
| Another crate | None found: `i_overlay` (our polygon backend) has no arcs, and no other Rust crate does arc booleans. |

## Consequences

**Positive**

- The same scene gives the same answer in mm and m (tests
  `*_in_any_unit` in `crates/algorithms/planar/overlay/tests/arc_overlay.rs`).
- A mutation that forces the scale back to 1 fails all three; removing
  the extent cap makes the backend panic on the existing 1e-18 audit test.

**Negative / costs**

- Still not exact: results are within floating-point error of the true
  answer, not bit-exact like the polygon path. F1 stays narrow and says so.
- The effective tolerance is within a factor of sqrt(2) of the requested
  one (power-of-two scaling), and a tolerance finer than about 1e-11 of the
  drawing's extent is honoured only down to that floor.

**Follow-ups / risks to watch**

- Replace the backend once ADR 0068's number layer exists (#155).
- Upstream: offering the crate a single scale-relative epsilon would help
  every user; worth an issue there.
