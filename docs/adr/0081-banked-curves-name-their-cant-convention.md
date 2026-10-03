# 0081 — Banked curves name their cant convention

- **Status:** Accepted
- **Date:** 2026-10-03
- **Deciders:** Friedrich Schrödter
- **Supersedes:** —

## Context

Issue #240 (consumer: openbimrs/ifc#93). A canted track centreline is
authored as three laws over plan distance: the plan, the vertical profile
and a cant law `D(d)`, the height of the left rail head above the right.
With the rail-head distance `b`, the bank angle is `psi = asin(D / b)`.
Before this, no Axiolid curve carried a roll about its tangent:
`Elevated3` pairs a plan with a profile only, `Intrinsic3` carries
curvature and torsion but no roll. So the section frame of a segmented
reference curve (IFC 4.3 `IfcSegmentedReferenceCurve`, cant from
`IfcAlignmentCant`) could not be represented.

On a grade the cant is ambiguous. IFC 4.3 gives `psi = asin(D / b)` and
the rail heights left and right of the profile, but does not say whether
`D` is a height difference measured vertically or the cross-fall of the
section square to the 3D tangent. The two readings differ by
`D (1 - cos theta)` at grade angle `theta`: `1.9994e-4 D` at 2%, 0.03 mm
on a 150 mm cant. Small, but systematic, and it compounds into every
object placed on the track.

## Decision

We will add `Curve3::Banked(Banked3)` and make the reading an explicit,
required field, `convention: BankConvention`, with no default. Both
readings are supported, exactly, and both keep the section frame
orthonormal:

- `base: Elevated3`, evaluated by plan distance `d`;
- `cant: CantLaw`, pieces laid end to end, each in its own
  `xi = s / length`: polynomial in `xi` (constant, linear, Bloss, Helmert
  as two quadratic pieces), half-cosine `D1 + dD (1 - cos(pi xi)) / 2`,
  sine `D1 + dD (xi - sin(2 pi xi) / (2 pi))`, and the Viennese bend
  `psi1 + dpsi xi^4 (35 - 84 xi + 70 xi^2 - 20 xi^3)`, which gives the
  bank angle directly (`D = b sin psi`);
- `pivot: CantLaw` of height pieces, the elevation `e(d)` of the rotation
  point above the profile (`0` about the centreline, `D / 2` about the low
  rail of a left-high cant);
- `rail_head_distance` `b`; `|D| > b` is refused by name;
- the point is `base(d) + e(d) z`, and `t` is that point's own unit
  tangent (its grade includes `e'(d)`).

With `n` the horizontal left normal and `u = t x n`, the section is
rolled about `t` by `rho`: lateral `l = cos(rho) n + sin(rho) u`, section
up `v = -sin(rho) n + cos(rho) u`, rail heads at `point +- (b / 2) l`.

| Convention | Roll `rho` | Vertical rise across the rail heads | Exact quantity |
| --- | --- | --- | --- |
| `TangentRotation` | `psi = asin(D / b)` | `D cos(theta)` | the bank angle about the tangent is `psi`, as on level track |
| `VerticalRise` | `asin(D / (b cos theta))` | `D` | the rail heads, `b` apart square to the tangent, differ by `D` vertically; refused when `|D| > b cos(theta)` |

`TangentRotation` is the reading of the IFC 4.3 cant as a rotation of the
gradient curve's frame about its tangent (the frame the consumer issue
writes, `(t, cos psi n + sin psi u, -sin psi n + cos psi u)`).
`VerticalRise` is the reading of `StartCantLeft`/`StartCantRight` as rail
heights measured vertically, as track geometry states cant. On level
track they coincide. The evaluator's `frame_at` returns this section frame
for a banked curve (`x` tangent, `y` section up, `z = -l`), which at zero
cant is the reference-up frame against `+Z`.

## Alternatives considered

| Option | Why not |
| --- | --- |
| Pick one reading and document the other | The issue's own maintainer decision: both have consumers, and a kernel choice would be silently wrong for the other (ADR 0066). |
| A default convention | Hides the choice behind a name that implies there is none; a caller must say which `D` it has. |
| A cant field on `Elevated3` | Changes every `Elevated3` literal (not additive), and an elevated curve with no cant is a meaningful value of its own. |
| `VerticalRise` keeping the rail-head line in the vertical plane through `n` at `psi` | Rise is exactly `D` and the angle is IFC's `psi`, but the lateral axis is then not square to the tangent: the frame is not orthonormal. |
| Frame on the centreline's tangent, ignoring the pivot's rate | `frame_at` would no longer agree with `tangent_at` where the pivot moves. |

## Consequences

**Positive**

- A segmented reference curve lowers exactly, with every cant form the
  consumer named, and the convention is visible in the value.
- Both conventions are pinned against hand-computed frames, including
  the `D (1 - cos theta)` difference at 2%.

**Negative / costs**

- A low-rail rotation through a Viennese bend has `e = (b / 2) sin(psi)`,
  which a pivot law of height pieces cannot state exactly; such a pivot
  is refused rather than approximated.
- Flattening a banked curve is measured, not certified (#232):
  `certifies_flattening3` reports it uncertified.
- A rigid motion would tilt `+Z`; `ExactBRep::transformed` refuses a
  banked curve, as it does an elevated one.

**Follow-ups / risks to watch**

- If a consumer needs the pivot as a fraction of the cant (to cover angle
  pieces), add it as a new pivot form rather than reinterpreting the law.

## Relation to existing code

- `crates/representations/analytic/curve/src/banked.rs`: the value, laws,
  conventions and `BankError`.
- `crates/algorithms/parametric/evaluate/src/banked.rs`: point, tangent,
  `BankedSection`; wired into `curve.rs`, `provider.rs` and `bound.rs`.
- `crates/representations/brep/src/transform.rs`: refusal by name.
