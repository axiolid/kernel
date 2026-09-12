# Certified boundary roots — the gate on curved-surface booleans

Standing scope note. Measured 2026-09-09 against `b69f3dc`.

## Why this document exists

An external reviewer observed that Axiolid has "no curved-surface boolean".
That is accurate. This note records WHERE the blocker actually is, because
the obvious guesses are wrong and cost real time to rule out.

## What is already built

The numerics are not the gap. Measured sizes in `axiolid-nurbs`
(10,236 LOC total, 6,799 of it certified):

| file | LOC | role |
|---|---|---|
| `certified_surface_bezier.rs` | 1244 | outward-rounded surface enclosures |
| `certified_curve_intersection.rs` | 1125 | certified curve/curve |
| `certified_curve_surface_intersection.rs` | 698 | **the blocker lives here** |
| `certified_bezier.rs` | 649 | interval arithmetic, `Interval {lo, hi}` |
| `certified_surface_surface_intersection.rs` | 640 | SSI driver |
| `certified_surface_inversion.rs` | — | point inversion |
| `certified_refinement.rs` | 509 | `RefinementBudget` |

Krawczyk root proofs, interval arithmetic, explicit work budgets, and
`BudgetExceeded` / `ProjectionStatus::BudgetExhausted` refusals all exist.

## Where the wall is

`split_surface_pair_certified` refuses the dual-boundary case — two patches
whose intersection segment ends on a boundary of BOTH — with
`IntersectionUnresolved`. It never reaches topology classification.

Measured chain:

```
split_surface_pair_certified
  -> intersect_surface_surface_certified
     -> trace_affine_pair
        -> collect_boundary_roots   (8 boundary curve/surface queries)
           -> intersect_curve_surface_certified
              -> Unresolved  =>  any_unresolved = true
        -> endpoints=0, any_unresolved=true
        -> AffineTraceOutcome::Unresolved(vec![])
```

Instrumented output for `xy_plane(-1,1)` against `xz_plane(-1,1)`:

```
SSI endpoints=0 any_unresolved=true      <- dual-boundary case
SSI endpoints=2 any_unresolved=false     <- the 9 supported cases
```

## It is a proof limit, not a budget limit

This is the part worth recording, because it rules out the cheap fix.
Escalating every knob leaves the verdict unchanged:

```
tol=1e-7   nodes=100000 depth=64 -> UNRESOLVED IntersectionUnresolved
tol=1e-10  nodes=100000 depth=64 -> UNRESOLVED IntersectionUnresolved
tol=1e-12  nodes=100000 depth=64 -> UNRESOLVED IntersectionUnresolved
tol=1e-14  nodes=100000 depth=64 -> UNRESOLVED IntersectionUnresolved
tol=1e-16  nodes=100000 depth=64 -> UNRESOLVED IntersectionUnresolved
tol=1e-7   nodes=100000 depth=32 -> UNRESOLVED IntersectionUnresolved
tol=1e-7   nodes=1000   depth=64 -> UNRESOLVED IntersectionUnresolved
```

(`max_refinement_work` is hard-capped at 100000, so larger budgets are
rejected by options validation rather than tried.)

Invariant at 1e-16 — below double epsilon. More work cannot help.

## Root cause

In `certified_curve_surface_intersection.rs`, a Krawczyk root is accepted
only when `certificate_meets_resolution` holds. A root lying exactly ON a
patch boundary sits at the edge of the parameter box, so the Krawczyk
operator cannot prove strict containment in the box interior — which is
what its contraction argument requires. The solver then contracts, finds
`contracted == current.parameters`, and pushes to `unresolved`:

```rust
if contracted == current.parameters {
    push_result(&mut unresolved, contracted)?;
    continue;
}
```

Subdivision cannot separate a root from a boundary it lies on, at any
depth. The refusal is correct given the method; the method is the limit.

## What closing it requires

Not a tuning change. A boundary-aware root certificate: restrict the
system to the boundary (one parameter pinned at a knot extreme, reducing
curve/surface to a 1-D problem in the remaining parameter) and prove the
root there with its own certificate, then confirm the pinned parameter is
exactly at the extreme rather than near it. The exact-arithmetic substrate
for that already exists — `exact_sum_is_zero` in the SSI file uses a
Shewchuk nonoverlapping expansion to decide the affine cross-term identity
without any tolerance.

Sequence, each step independently gateable:

1. **Boundary-restricted certificate** in `certified_curve_surface_intersection.rs`.
   Pin one parameter, certify the reduced system, prove the pin exactly.
   Verify with a mutation that a near-boundary root is NOT accepted as on-boundary.
2. **Plumb it** through `collect_boundary_roots` so an on-boundary root
   stops setting `any_unresolved`.
3. **Dual-chord topology** in `trimmed_intersection_classify.rs`: the
   `(true, true, _, _)` arm, which partitions BOTH patches into 4 faces.
   This is where `UnsupportedEndpointOwnership` stops being the answer.
4. **Lift the affine restriction**. `is_exact_single_span_affine` requires
   `u_degree == 1 && v_degree == 1`, a single span, and no weights. Genuinely
   curved booleans need the non-affine path: traces become curves, not the
   line segments `add_line_edge` currently assumes.

Steps 1–3 unlock the dual-boundary PLANAR case. Step 4 is the actual
curved-surface boolean and is the largest of the four by a wide margin.

## Honest status

Do not describe the curved-surface boolean as close. Steps 1 and 2 are
contained and well-understood. Step 3 is ordinary topology work. Step 4 is
a project: every assembly path that assumes a straight intersection edge
(`add_line_edge`, `add_pcurve` over `Interval::UNIT`) has to learn curved
trace geometry first.


## Step 1 outcome (measured)

Step 1 is DONE and the proof gap is closed. `krawczyk_root_on_edge` pins
the parameter that sits on a domain edge and certifies the reduced
two-unknown system, where the root is interior.

Measured before: the dual-boundary case reached
`endpoints=0 any_unresolved=true` -- no boundary root could be certified.

Measured after: `endpoints=4 any_unresolved=false`. Four roots certified,
zero unresolved. The certification blocker is gone.

### The next blocker is NOT what this doc predicted

`trace_affine_pair` refuses `endpoints > 2`:

```rust
if any_unresolved || endpoints.len() == 1 || endpoints.len() > 2 {
    return Ok(AffineTraceOutcome::Unresolved(Vec::new()));
}
```

Four endpoints is CORRECT for the dual-boundary case: two planar patches
crossing symmetrically meet four domain edges, because the chord runs
boundary-to-boundary on BOTH patches. The guard assumes a single trace has
exactly two endpoints, which holds only when one patch contains the
chord's interior.

So the remaining work is trace ASSEMBLY, not certification: pair the four
certified endpoints into the correct chord per patch, then split both
patches instead of one. That also means `CertifiedTrimmedSurfacePair3`'s
`split_surface: SurfacePairMember` / `unsplit_face` shape must generalise --
in a dual chord there is no unsplit face. That is a breaking type change
and belongs with steps 2-3, not smuggled into step 1.

### Honest coverage note

The discarded-row verification inside `krawczyk_root_on_edge` is NOT
covered: every off-patch input reachable today is rejected earlier by
`residual_excludes_zero`, so deleting the check fails no test. It is
retained as defence and marked as unverified in the source.

## Step 2 outcome (measured)

Shipped: the dual-boundary chord now splits BOTH patches, four closed
trimmed faces sharing one intersection edge.

My step 1 prediction was WRONG in two ways, both corrected by measuring:

1. I predicted the blocker was assembling a 4-endpoint dual chord.
   Measurement showed the four endpoints were only TWO distinct points,
   each reported twice -- once from scanning each surface's boundaries.
   A chord ending on a boundary of BOTH patches is found from both
   sides. Before step 1 nothing was certifiable there, so the
   duplication had never surfaced. The fix was de-duplication, not
   assembly.

2. I predicted this forces a breaking change to
   `CertifiedTrimmedSurfacePair3`. It does not. The enum is
   `#[non_exhaustive]`, so a new `DualSplit` variant carrying a new
   `CertifiedDualSplitSurfacePair3` was added alongside `Split`.
   The single-split type keeps its exact meaning -- `unsplit_face` and
   `embedded_curve` stay honest because they are absent from the dual
   type, where no unsplit face exists.

### What the dual type does NOT claim

`CertifiedDualSplitSurfacePair3` has no `embedded_curve`: with both
patches partitioned there is no containing face to embed a dangling
edge into. Every one of the four loops uses the shared edge, asserted
in test.

### Remaining

Still planar-only: `is_exact_single_span_affine` requires degree 1,
single span, no weights. Genuinely curved surface/surface intersection
is unchanged by steps 1 and 2 and remains the large piece of work.
