# axiolid-evaluate instructions

Purpose: the scalar evaluation oracle for parametric geometry (ADR 0012, ADR 0036).

Allowed internal dependencies: `axiolid-core`, `axiolid-contracts`,
`axiolid-curve`, `axiolid-surface`. Do not add mesh, spatial, measure, or
provider dependencies — the point of this package is that a parametric consumer
(NURBS, CAD) acquires evaluation without the `axiolid-reference` umbrella graph.

## Module ownership

`curve.rs` native-domain evaluation, derivatives, jets, adaptive flattening;
`surface.rs` evaluation, partials, normals, jets, elementary inversion;
`arc_length.rs` arc-length evaluation of intrinsic curves and of the
plan-plus-elevation composition;
`frenet.rs` frame, point and tangent of a space curve from curvature AND
torsion, by Magnus expansion on SO(3);
`nurbs.rs` shared private spline-axis machinery.

## Invariants

- A 3D frame law does NOT integrate like a 2D one. In the plane, heading is
  `theta_0 + int k` because angles commute. In space the frame obeys a matrix
  ODE on SO(3) and `R(s) != exp(int Omega)` unless `tau/k` is constant. Use the
  Magnus step in `frenet.rs`; dropping its commutator term costs ~5000x accuracy
  on a clothoid-with-torsion (ADR 0061).
- A quadrature budget must come from TOTAL VARIATION of the law, never from its
  signed integral: a zero-mean oscillation integrates to zero and would buy one
  panel for a curve that swings through many radians (ADR 0061).
- An intrinsic curve's HEADING is exact in closed form; its POSITION is not
  elementary and is Gauss-Legendre quadrature, subdivided by total turning and
  bounded so a malformed law refuses rather than hangs. Pin any change to it
  against an independent closed form (Fresnel for the clothoid, the elementary
  arc for constant curvature) — never against another run of the quadrature
  (ADR 0060).
- No feature gates. An oracle that varies by feature is not an oracle.
- No intrinsics or threading; stay obviously correct in preference to fast.
- `axiolid-reference` re-exports `curve` and `surface` unchanged. Renaming or
  reshaping a public item here silently breaks `axiolid_reference::curve::*`
  callers, so treat those paths as part of this package's public surface.
