# Curve evaluation as a capability contract

Status: accepted
Date: 2026-09-16
Issue: axiolid/kernel#106

## Context

Curve evaluation was the only kernel capability without a contract
crate. A consumer needing "a point at a distance along this curve" had to
depend on `axiolid-evaluate`, which names itself a backend:

```rust
backend: BackendId::new("axiolid-evaluate"),
```

So there was no way to NAME the capability without binding an engine to
it. Downstream (`openbimrs/ifc`), an architecture gate classifies kernel
crates as neutral representations or execution providers; bridges may
depend only on the former. Depending on the evaluator failed that gate,
correctly: it chooses an engine rather than requesting a capability.

The concrete blocked feature is IFC4x3 `IfcLinearPlacement` -- placing a
sign or drainage structure at a distance along an alignment.

## Decision

Add `axiolid-curve-evaluate-contract`, shaped exactly like the mesh
contracts: a `CurveEvaluator: Backend` trait, a capability id
`org.axiolid.geometry.curve-evaluate.v1`, and a conformance suite.
`axiolid-evaluate` implements it as `ReferenceCurveEvaluator`.

Two departures from the shape proposed in the issue follow.

## This is not only packaging

The issue states every function already exists and that this is
"packaging, not new mathematics". That holds for `Intrinsic` and
`Elevated`, which are authored against distance. It does not hold for
the rest of `Curve3`.

`evaluate3` takes a native PARAMETER: an angle for a circle, a knot
value for a B-spline, a segment index for a polyline. Distance and
parameter coincide only for an arc-length-parameterised family. No
arc-length reparameterisation exists in the kernel, and none is added
here.

The contract therefore carries a `DistanceConvention`, and providers
report which distance they measure:

| Family | Convention | Why |
| --- | --- | --- |
| `Line` | `ArcLength3d` | `d / |direction|`, exact |
| `Circle` | `ArcLength3d` | `d / radius`, exact |
| `Intrinsic` | `ArcLength3d` | already arc length |
| `Elevated` | `PlanDistance` | authored against plan distance |
| `Ellipse` | `Unsupported` | needs elliptic integrals |
| `BSpline` | `Unsupported` | needs numeric inversion |
| `Polyline` | `Unsupported` | needs cumulative traversal |

A non-unit `direction` is the trap worth naming: import adapters
preserve it, so treating distance as the parameter scales every
placement by `|direction|`. A doubled direction halves the distance.

`PlanDistance` is reported rather than converted. Plan distance and 3D
arc length differ by `sqrt(1 + grade^2)`: 0.125 m per 100 m at 5%,
0.499 m per 100 m at 10%. A station on a drawing IS a plan distance, so
converting would be wrong; the caller must know which they asked for.

## The frame is reference-up, NOT Frenet

The issue calls `frame_at` "the one genuine API request" and is right
that the convention must be decided once, upstream. The obvious choice
is the Frenet frame -- `frenet_frame` already exists and returns one.
It is the wrong choice for placement, and measurably so.

The Frenet normal points toward the centre of curvature. On an
alignment's vertical profile that direction FLIPS at an inflection:

```text
crest    -> normal z = -1.0     (points down)
sag      -> normal z = +1.0     (points up)
straight -> UNDEFINED           (zero curvature, 0/0)
```

A sign placed with a Frenet frame is upright on one side of a crest and
upside-down on the other, and undefined on the tangent between them --
which is most of a real alignment. This is exactly the "wrong convention
tilts a road sign rather than failing loudly" failure the issue warns
about.

The contract instead defines:

```text
x = unit tangent
z = normalise(tangent x up_reference)     (right)
y = z x x                                 (up, re-orthogonalised)
```

with `up_reference` defaulting to global +Z. This is stable across
crest, sag and straight alike, and is the convention surveying and
highway software already use. Verified: up stays within 0.1 degrees of
vertical across a crest-to-sag inflection where the Frenet normal
inverts.

When the tangent is parallel to `up_reference` the cross product
vanishes and roll is genuinely undefined. That REFUSES rather than
picking an arbitrary roll, because a silently-rolled placement is
undetectable downstream.

`frenet_frame` remains available for curve analysis, where the
centre-of-curvature direction is the point. The two frames answer
different questions and both are kept.

## Consequences

- A consumer can name curve evaluation without linking an engine, which
  is what unblocks `IfcLinearPlacement` in `openbimrs/ifc`.
- `Ellipse`, `BSpline` and `Polyline` refuse by name at a distance API.
  They remain fully evaluable through `evaluate3` at their native
  parameter; only the DISTANCE question is refused. Closing that gap
  means arc-length reparameterisation, which is its own decision.
- The closure profiles `c-abi-profile`, `cad-exact`,
  `parametric-curves` and `rust-facade-application` grow by exactly one
  contract crate. Acknowledged deliberately in
  `architecture/closure-profiles.toml` rather than by relaxing the gate.
- `capability_ids::ALL` did not previously list
  `POINTCLOUD_RECONSTRUCTION`. That is left as found rather than
  silently fixed here; it deserves its own change.

## Alternative rejected

The issue offers to hand callers the exact `Elevated3` and let the
application evaluate it. That keeps the bridge pure but fragments the
frame convention across consumers -- which, given the Frenet flip above,
means several of them would get it wrong in a way no test catches.
Deciding the convention once is the substance of this change; the
packaging is the cheap part.

## Addendum: distance versus native parameter

The first version of this contract took a bare `Scalar` distance. That
was incomplete, raised on #106 after the fact.

There are two independent axes, and the original design named only one:

1. WHAT a distance measures -- 3D arc length or plan distance. Answered
   by `DistanceConvention`.
2. WHETHER the caller's number is a distance at all. Not answered.

IFC4x3 makes the second axis explicit:
`IfcPointByDistanceExpression.DistanceAlong` is an
`IfcCurveMeasureSelect`, so an authored value is EITHER an
`IfcNonNegativeLengthMeasure` or an `IfcParameterValue`, and a file says
which. STEP carries the same distinction.

With a bare `Scalar` a consumer holding a parameter had nothing to
prevent passing it as a distance. On a circle of radius 4, `1.5` as a
parameter is 1.5 rad round; as a distance it is 0.375 rad -- about 86
degrees apart, both finite and plausible, undetectable downstream.

### Decision

The three methods take a `CurveMeasure` instead of a `Scalar`:

```rust
pub enum CurveMeasure {
    Distance(Scalar),
    Parameter(Scalar),
}
```

Rejected: a `DistanceConvention::NativeParameter` variant, which was the
other option offered. That enum answers what a DISTANCE measures, and a
parameter is not a distance -- the variant would make the enum answer two
questions and let `distance_convention` describe something that is not
one. Rejected too: separate `*_at_parameter` methods, which leave the
caller branching, and a caller that branches can branch wrongly.
Carrying the method of measurement in the value makes the mistake
unrepresentable rather than merely documented.

### The refusal table was also too strict

A second defect surfaced while implementing this. The original reply to
#106 told consumers that `Ellipse`, `BSpline` and `Polyline` remain
evaluable "through `evaluate3` at their native parameter". That advice
was unusable: `evaluate3` lives in `axiolid-evaluate`, the engine crate,
and the whole point of the contract is that a bridge may not depend on
it. The parameter route existed only outside the contract.

`CurveMeasure::Parameter` fixes that too. It needs no conversion and no
convention, so it is answered for EVERY family, including those whose
arc length has no closed form:

| Family | `Distance` | `Parameter` |
| --- | --- | --- |
| `Line`, `Circle` | exact | yes |
| `Intrinsic` | exact (identity) | yes, same value |
| `Elevated` | exact, plan distance | yes, same value |
| `Ellipse`, `BSpline`, `Polyline` | refused | yes |

So a refused DISTANCE no longer means a curve is unreachable through the
contract; only that one method of measurement is unavailable for it.

### Evidence

- The same number gives points ~1.5 apart on a circle depending on the
  method of measurement, each pinned against its own closed form.
- On an arc-length-parameterised curve the two routes AGREE exactly,
  pinned so the distinction is not over-enforced.
- An ellipse refuses a distance and answers a parameter, with a frame.
- Mutation 3/3: collapsing `Parameter` into `Distance`, routing a
  parameter through the distance conversion, and dropping the finiteness
  guard all die.

Landed before the crate was published, so no consumer saw the bare
`Scalar` signature.
