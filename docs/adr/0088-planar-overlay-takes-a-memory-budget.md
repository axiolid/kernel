# 0088 — The planar overlay takes a memory budget from ExecutionOptions

- **Status:** Accepted
- **Date:** 2026-10-20
- **Deciders:** Friedrich Schrödter
- **Supersedes:** —

## Context

Issue #338 (consumer: axioval/engine). A soup of grown triangles whose
union is one polygon needed gigabytes in `union_soup` and `overlay`: 128
facets took 4 GB, and 256 facets and the soup's self-union failed to
allocate, which aborts the host process. Profiling found two causes, both
fixed exactly with unchanged results (see the crate's changelog): the
chain reduction gave up on the whole soup when one cycle of its single
walk crossed itself, and the subdivision recorded every pair of edges
meeting at a shared vertex and every pair of repeated segments. A caller
still cannot tell in advance how large an arrangement a soup makes, and
an allocation failure is not a refusal it can handle. ADR 0017 and #308
give boolmesh a memory budget under `ExecutionOptions`; the overlay had
none, and `axiolid-overlay` did not depend on `axiolid-contracts`, where
`ExecutionOptions` lives.

## Decision

`axiolid-overlay` depends on `axiolid-contracts` and adds
`overlay_within` and `union_soup_within`, which take `OverlayOptions` and
an `&ExecutionOptions`: its tolerance, and its memory budget. With a
budget, the subdivision charges what it is about to hold before
allocating it (the operands by vertex, the edges, every recorded meeting
point, the pieces, each piece's rings) and refuses with
`OverlayRefusal::BudgetExceeded { needed, budget }` as soon as the charges
pass it. Without a budget nothing new is refused, and every existing entry
point is unchanged.

- The refusal is a new `#[non_exhaustive]` error, `OverlayRefusal`, with
  `Invalid(OverlayError)` for every refusal the unbudgeted calls make.
  `OverlayError` is exhaustive, and a variant added to it would break
  every match on it.
- The charges are upper bounds measured with a counting allocator
  (`tests/grown_soup.rs`): under any budget, admitted or refused, a call
  holds no more heap than the budget, and a budget of four times the
  unbudgeted high-water mark is admitted.
- No other execution option is read; a cancellation token is not polled.

## Alternatives considered

| Option | Why not |
| --- | --- |
| A `BudgetExceeded` variant on `OverlayError` | A breaking change to an exhaustive enum every consumer matches on. |
| A budget in bytes on `OverlayOptions` | A second budget vocabulary beside `ExecutionOptions`; the triage of #338 puts it under `ExecutionOptions`, as boolmesh's is. |
| Estimate from the candidate pairs alone and refuse up front | Box overlaps far outnumber crossings in a soup of long, nearly parallel edges, so the estimate would refuse calls that fit. Charging each point before recording it stops as early, at the real count. |

## Consequences

**Positive**

- A caller with a memory limit gets a named refusal instead of an aborted
  process, before the subdivision allocates past the budget.
- The closure of every profile that holds `axiolid-overlay` already held
  `axiolid-contracts`, so no declared closure changes.

**Negative / costs**

- `axiolid-overlay`'s direct dependents (`axiolid-project`,
  `axiolid-brep-boolean`, and `axiolid-measure` with its `exact` feature)
  now build `axiolid-contracts` too.
- The budget covers the subdivision, not the operands' own storage or the
  output; snapping and the chain reduction run before the subdivision is
  charged.

**Follow-ups / risks to watch**

- The charges' constants are measured on grown soups, crossing squares and
  strip intersections. A scene whose exact arithmetic builds unusually
  large coefficients could hold more per point; the test is the place to
  add it.

## Relation to existing code

- `crates/algorithms/planar/overlay/src/budget.rs`: `OverlayRefusal`, the
  charges and their constants.
- `crates/algorithms/planar/overlay/src/exact_arc/arrangement.rs`:
  `build_charged`, where the subdivision charges.
- `crates/algorithms/planar/overlay/src/lib.rs`: `overlay_within`,
  `union_soup_within`.
- `crates/contracts/common/base/src/execution.rs`: `ExecutionOptions`.
