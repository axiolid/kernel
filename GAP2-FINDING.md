# Gap 2 finding: the production boolean is NOT the weak link

## What I claimed in the audit

> "Ours computes intersections in f64 ... our exactness stops at the test
> boundary ... wiring exact predicates into boolmesh's production path is the
> priority, because it decides whether the speed advantage we publish is a
> fair comparison."

That recommendation was wrong, and the measurements below are why. Recording
it here rather than quietly dropping it.

## What the probes actually measured

### Probe 1 — near-coincident faces at origin scale

Swept a cutter face across a host face, 1e0 down to 1e-15, comparing against
the exact analytic volume:

```
eps=1e-1   err=1.110e-16
eps=1e-8   err=1.110e-16
eps=1e-15  err=1.110e-16
```

One ULP at every separation. No lost slab, no flipped predicate. The
`interpolate`/`intersect` pair carries an upstream comment claiming they are
"carefully designed to minimize rounding error and to remove it at edge
cases"; at this scale the claim holds up.

Of the ten sign decisions in the boolean kernels, six are direct coordinate
comparisons (already exact -- comparing two f64s is not a rounding problem),
and only four involve an interpolated value.

### Probe 2 — large coordinates

```
base=0     half-cut exact
base=1e3   half-cut rel err 1.192e-8
base=1e6   half-cut rel err 1.626e-2
base=1e7   half-cut rel err 9.068e-1   <-- 91% wrong
```

This looks like the smoking gun for gap 2. It is not.

### Probe 3 — attribution

Measuring the INPUT mesh before any boolean runs:

```
base=0     input_mesh_volume rel_err=4.337e-16
base=1e3   input_mesh_volume rel_err=5.862e-10
base=1e6   input_mesh_volume rel_err=5.668e-3
base=1e7   input_mesh_volume rel_err=2.528e-1   <-- 25% wrong ALREADY
```

A quarter of the error at base 1e7 exists before the boolean is called. The
divergence tracks ULP-at-magnitude almost exactly, which is the signature of
catastrophic cancellation in the divergence-theorem volume sum: it adds terms
of order 1e21 to produce an answer of order 1e-3.

## Conclusion

The large-coordinate failure is real and worth fixing, but it is NOT
"the boolean needs exact predicates". It is coordinate magnitude: geometry
far from the origin loses relative precision everywhere -- in the volume
measure, in the mesh, and only then in the boolean.

The standard fix is a local origin / RTC (relative-to-centre) offset, which
ifc-lite implements (`router/rtc_offset.rs`) and axiolid does not. That is a
different and cheaper change than a filtered-exact cascade through the
boolean kernels, and it fixes the measure and the mesh too.

## Revised recommendation

1. **Do not** rewrite the boolean kernels for exactness on this evidence.
   Measured at origin scale, they are already at one ULP.
2. **Do** add an RTC/local-origin facility, and a gate that fails when a
   fixture's own input mesh cannot be measured to tolerance -- that would
   have caught this class of bug without anyone reading the boolean at all.
3. Keep the probes as a regression test with HONEST thresholds: exact at
   origin scale, documented degradation at 1e6+, so the limit is stated
   rather than discovered by a user with a site survey in national grid
   coordinates.
