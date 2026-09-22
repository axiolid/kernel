# RTC step 1: does survey-scale quantisation actually break anything?

Probe result for the question left open by `fix(measure)`: the input mesh's
own vertices are quantised to ~2 nm at base 1e7, and no downstream fix can
recover that. Is it a real problem?

**Answer: not for the reason I predicted. Two of my three claims were wrong.**

## Claim 1: "coincident faces get snapped to different values, creating slivers"

**Refuted.** 0 of 5 magnitudes lost exact coincidence.

```
COINCIDENCE base=0e0   drift=0.000e0   ulp=2.220e-16  sign=Zero
COINCIDENCE base=1e3   drift=0.000e0   ulp=1.137e-13  sign=Zero
COINCIDENCE base=1e5   drift=1.455e-11 ulp=1.455e-11  sign=Zero
COINCIDENCE base=1e6   drift=0.000e0   ulp=1.164e-10  sign=Zero
COINCIDENCE base=1e7   drift=0.000e0   ulp=1.863e-9   sign=Zero
```

Even where two *different* arithmetic paths to the same intended coordinate
drifted by a full ULP (base 1e5), `orient3d` still certified `Zero`. The
reason is structural: rounding is deterministic and the predicate is exact
for any finite binary64 input, so coincidence authored consistently survives
regardless of magnitude.

The end-to-end check agrees. Two boxes sharing a face exactly, union at every
magnitude:

```
FLUSH base=0e0 .. 1e7   volume=2.00000000000000000  rel_err=0.000e0
FLUSH summary: 0 of 5 magnitudes wrong
```

Exact at 1e7. The sliver failure mode does not occur.

## Claim 2: "rebasing makes the filter hit its cheap path more often, so it's faster"

**Refuted.** Escalation rate is flat.

```
ESCALATION base=0e0 .. 1e7   1/2000 = 0.1% fell back to exact
```

Identical at every magnitude. The filter's error bound scales with the
operands, so its decision threshold scales too -- magnitude does not push it
toward the slow path. **The performance argument for RTC was wrong, and the
claim should not be repeated.**

## What IS real: thin features at large coordinates

```
PLATE base=1e6  t=1e-4  rel_err=1.667e-5   <-- WRONG
PLATE base=1e7  t=1e-2  rel_err=1.758e-3   <-- WRONG
PLATE base=1e7  t=1e-4  rel_err=1.681e-5   <-- WRONG
```

But the attribution probe shows this is **not the boolean**:

```
ATTR base=1e7  t=1e-2  input_vol_rel_err=6.913e0    corner_drift=0.000e0
ATTR base=1e7  t=1e-4  input_vol_rel_err=1.614e5    corner_drift=0.000e0
ATTR base=1e7  t=1e-6  input_vol_rel_err=2.965e8    corner_drift=0.000e0
```

`corner_drift = 0` everywhere: the authored geometry is **exact**. The plate's
vertices are where they should be. What fails is *measuring* it -- the
divergence-theorem sum over a 1e-6 m feature at 1e7 m cancels catastrophically,
giving errors up to 2.9e8 relative.

Note the direction: the measurement error (2.9e8) is many orders worse than
the boolean's output error (1.8e-3). The boolean is more accurate than the
tool used to check it.

`fix(measure)` already re-bases to the mesh's first vertex, which fixes this
for a mesh whose own extent is small. It does not help when the *feature* is
tiny relative to the mesh -- a 1e-6 m plate is below the 1.9e-9 m grid only
by a factor of 500, so its own corners span just a few hundred ULPs.

## Conclusion

The originally-stated motivation for an RTC/local-origin facility -- exact
predicates breaking on coincident geometry -- **does not reproduce**. Neither
does the performance argument.

The real limit is narrower and different: *mass properties of features whose
size approaches the representable grid at their location*. That is a
measurement conditioning problem, already partly addressed, and it does not
require a coordinate-frame facility, an API change, or an offset ownership
model.

**Recommendation: do not build the RTC facility.** Document the measured
limit, keep these probes as regression tests, and revisit only if a real
workload produces a failure that these probes do not cover.


## Reproducing

```
cargo test -p axiolid-predicates --test survey_scale -- --nocapture
cargo test -p axiolid-mesh-boolean-boolmesh --test survey_scale -- --nocapture
```

## Correction: the per-triangle follow-up was measured against the wrong baseline

An earlier revision of this document recommended re-basing the volume sum
per-triangle instead of per-mesh, citing an eleven-orders-of-magnitude
improvement. **That recommendation was wrong and is withdrawn.**

The comparison used `boolmesh/tests/support.rs::volume` as the "per_mesh"
baseline. That helper sums about the WORLD origin. Production
`volume_properties` has summed about a local origin since the RTC fix
(`fix(measure): sum mass properties about a local origin`), so the numbers
compared a per-triangle variant against code that no longer exists in the
measurement path.

Measured against production directly, on the same thin-plate fixtures:

```
PLATE base=1e7 t=1e-2 production=0.000e0    per_triangle=9.313e-10
PLATE base=1e7 t=1e-4 production=1.355e-16  per_triangle=5.069e-10
PLATE base=1e7 t=1e-6 production=0.000e0    per_triangle=1.563e-10
WORST                 production=1.355e-16  per_triangle=9.313e-10
```

Production is exact or within one ULP at every magnitude and thickness
tested. The per-triangle variant is **worse at every non-zero magnitude**,
by up to six orders of magnitude.

Why: per-mesh re-basing subtracts a nearby origin once, so every coordinate
entering the cross product is already edge-sized. The per-triangle variant
still forms `a . (ab x ac)` with `a` at full world magnitude -- it shrinks
two operands and leaves the third large. Per-mesh re-basing shrinks all
three.

**There is no remaining measurement defect to fix.** The thin-plate error
reported in the step-1 probe belongs to the test helper, not to any shipped
code path. The conclusion of this document stands unchanged and is now
stronger: do not build the RTC facility, and do not change
`volume_properties` either.
