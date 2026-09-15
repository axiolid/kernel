# Torsion: space curves from their natural equations

Status: accepted
Date: 2026-09-15

Follows ADR 0060, which deferred torsion as "not what an alignment states".

## Context

ADR 0060 added `Curve3::Elevated`: a planar layout plus an elevation law. That
covers alignments, where the vertical profile is authored independently of the
plan. It does not cover a curve whose *frame* twists — a helical stair stringer,
a spiral ramp handrail, a coil — where the shape is stated as curvature and
torsion along arc length.

`Curve2::Intrinsic` already carries curvature laws exactly. The obvious move is
to add a torsion field to it. That is wrong, and the reason is the whole content
of this ADR.

## The 2D trick does not generalise

In the plane the frame is one angle and the shape equation integrates:

```text
theta(s) = theta_0 + int_0^s k
```

Angles commute, so heading has an elementary closed form and only *position*
needs quadrature. That is exactly what ADR 0060 shipped.

In space the frame is a rotation obeying Frenet–Serret:

```text
R'(s) = R(s) Omega(s),   Omega = [[0, -k, 0], [k, 0, -tau], [0, tau, 0]]
```

This is a matrix ODE on SO(3). Its solution is a *product integral*, and

```text
R(s) != exp(int_0^s Omega)
```

except when generators at different arc lengths commute — which happens exactly
when `tau/k` is constant, i.e. a helix or a plane curve. So in 3D the **frame
itself** needs integration, not just the position. A torsion field bolted onto
`Intrinsic2` would invite precisely the wrong formula.

Measured, on `k` linear and `tau` constant over 40 m at the panel budget the
code picks: using `exp(int Omega)` lands 2.7e-2 from the reference where the
correct scheme lands 5.2e-6 — a factor of **5000**.

## Decision

Add `Intrinsic3` as its own value type, and `Curve3::Intrinsic` to carry it.
Both laws reuse `CurvatureLaw`, which is already the general "scalar function of
arc length" in this crate — a second near-identical enum named `TorsionLaw`
would be duplication with a different label.

Evaluate with a **fourth-order Magnus expansion** on SO(3), in
`axiolid-evaluate::frenet`:

```text
M = (h/2)(Omega_1 + Omega_2) - (sqrt(3) h^2/12)[Omega_2, Omega_1]
R(s0 + h) = R(s0) exp(M)
```

with `Omega_i` sampled at the two-point Gauss nodes and `exp` by Rodrigues.

Two properties follow structurally, not by tuning:

- **Orthonormality at any step size.** `exp` of a skew matrix is a rotation and
  a product of rotations is a rotation, so the frame cannot drift off SO(3).
  At 100 panels, `|R^T R - I|` is 2.6e-15 here versus 8.4e-10 for RK4.
- **Zero torsion reproduces the planar answer exactly**, because the generator
  then has no `tau` component and the motion stays in the start plane.

Position integrates the tangent `T(u) = R(u) e_x` by 8-point Gauss–Legendre on
the same panels.

## A latent bug this surfaced in ADR 0060

The 2D panel budget used `total_turning()` — the **signed** integral of `k`.
Over a whole number of periods of a zero-mean oscillation that is zero, so a
violently wiggling curve was budgeted **one panel**.

Measured on `k(s) = 2 sin(10 s)` over `[0, pi]`: signed turning is 0, one panel,
endpoint 2.1e-1 wrong. The fix budgets from the total variation `int |k|`,
added as `Intrinsic2::turning_variation_bound`, which buys 26 panels and lands
the endpoint to 4.9e-15. Exact for `Constant`; a triangle-inequality upper bound
for the other families, which is what a budget wants — overestimating costs
panels, underestimating costs correctness.

This was shipped and wrong in ADR 0060, and is fixed here with a regression test.

## Consequences

- `Curve3` now carries helices and general space spirals as exact values.
- A start frame that is not a right-handed orthonormal triad is **refused**, not
  silently re-orthonormalised: the integrator propagates it by rotations, so a
  start frame that is not a rotation makes every downstream frame meaningless.
- Arc length outside `[0, length]` is refused rather than extrapolated.
- A malformed law that would need more than `MAX_PANELS` refuses rather than
  hanging.
- `Intrinsic3::is_helical` names the commuting case: constant `k` and `tau`,
  where the product integral collapses to one exponential and the curve has an
  elementary closed form.
- **Not done here:** no `CurveRelation` composition of an `Intrinsic3` with
  anything, no trimming, no offsetting, and no lowering from any file format.
  `Curve3::Intrinsic` is a representable value with an evaluator, nothing more.

## Evidence

- Helix position and tangent against the **Darboux closed form** — a screw about
  `(tau, 0, k)/|omega|`, not about `z`. Agreement better than 1e-9 at six
  stations over two full turns.
- Zero torsion against the 2D path, itself pinned to Fresnel: agreement to 1e-9,
  and `|z| < 1e-12`.
- A varying law against 20,000 locally-exact screw steps — a reference that
  never uses the Magnus expansion, so it cannot move when the scheme changes.
- Mutation 6/6 killed: commutator dropped, torsion ignored, tangent frozen at
  the panel start, non-orthonormal start accepted, handedness check removed, and
  the budget reverted to signed turning.
