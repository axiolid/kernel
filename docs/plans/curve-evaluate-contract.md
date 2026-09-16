# Curve evaluation contract (issue #106)

Status: done (ADR 0063)

## Goal

Let a consumer NAME curve evaluation without binding an engine, so an
IFC bridge can request "point and frame at a distance" behind a feature
gate, exactly as it already does for mesh compilation and booleans.

## What the issue got right

- Curve evaluation is the only capability with no contract crate.
- `Backend`, `CapabilityId` and the conformance pattern already exist.
- A contract-only dependency is genuinely cheap (20 vs 51 crates).
- `frame_at` is the real request: a tangent alone does not fix roll.

## What it got wrong: this is NOT only packaging

The proposed trait takes a DISTANCE. Three findings say that cannot be
a thin re-export of existing functions.

### 1. The Frenet frame is unsafe for placement

Verified numerically (`scratch/framecheck.py`). On an alignment whose
plan is straight and whose profile is a crest then a sag:

```
crest     -> Frenet normal z = -1.0
sag       -> Frenet normal z = +1.0
straight  -> UNDEFINED (zero curvature)
```

The Frenet normal FLIPS at an inflection and is undefined on any
straight run. A sign placed with it would be upright on the crest,
upside down in the sag, and unplaceable on the straight between them.
This is precisely the "wrong convention tilts a road sign" failure the
issue warns about -- so the convention must NOT be Frenet.

### 2. The convention that works is reference-up

`right = normalise(tangent x up)`, `up' = right x tangent`, with `up`
the global +Z. Same check, same alignment:

```
crest/sag/straight -> up_z stays ~0.9988..1.0, right stays [0,-1,0]
```

Stable through inflections and across straights. Degenerate only when
the tangent is parallel to the reference (a truly vertical curve),
where `|tangent x up| = 0` and the frame must be REFUSED, not guessed.

### 3. Distance is not the parameter, and for Elevated3 it is not even
    3D arc length

No arc-length reparameterisation exists anywhere in the kernel. For
`Elevated3` the native parameter is PLAN distance, and the 3D arc
length differs by `sqrt(1 + grade^2)`:

```
grade  2%  ->  0.020 m drift per 100 m
grade  5%  ->  0.125 m drift per 100 m
grade 10%  ->  0.499 m drift per 100 m
```

Silently treating one as the other misplaces a drainage structure by
half a metre per 100 m. The contract must SAY which it means.

## Design

New crate `crates/contracts/operations/curve-evaluate`,
`axiolid-curve-evaluate-contract`, mirroring the mesh-boolean layout
(`lib.rs` + `contract.rs` + `conformance.rs`).

```rust
pub trait CurveEvaluator: Backend {
    fn distance_convention(&self, curve: &Curve3) -> DistanceConvention;
    fn point_at(&self, curve: &Curve3, distance: Scalar) -> GeomResult<Point3>;
    fn tangent_at(&self, curve: &Curve3, distance: Scalar) -> GeomResult<Vec3>;
    fn frame_at(&self, curve: &Curve3, distance: Scalar) -> GeomResult<Frame3>;
}
```

`DistanceConvention` is the honesty valve: `ArcLength3d`,
`PlanDistance`, or `Unsupported`. A caller that needs true 3D arc
length on a graded alignment can see it will not get it, instead of
discovering the drift on site.

Exactness tiers, following the house two-tier rule:

| Family | Distance recoverable | How |
| --- | --- | --- |
| `Line` | exact | `t = d / |direction|` (direction may be non-unit) |
| `Circle` | exact | `angle = d / radius` |
| `Polyline` | exact | walk segments, interpolate the remainder |
| `Intrinsic` | exact, identity | the parameter already IS arc length |
| `Elevated` | exact in PLAN distance | both halves authored against it |
| `Ellipse` | refused | needs elliptic-integral inversion |
| `BSpline` | refused | needs numeric arc-length inversion |

Refusing by name beats returning a parameter-as-distance lie.

## Validation

- Frame stays upright across a crest/sag inflection and on a straight.
- Vertical tangent refuses rather than returning a degenerate frame.
- Circle: distance d lands exactly at angle d/r; full turn = TAU*r.
- Line with a NON-unit direction lands at true distance.
- Ellipse and BSpline refuse, and say why.
- Conformance suite runnable by any future provider.
- Mutation testing on every exactness and refusal claim.
