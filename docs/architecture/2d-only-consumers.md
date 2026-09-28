# 2D-only consumers

A plan viewer or drawing tool does **not** need Axiolid's solid-modeling stack. The
machine-checked `2d-curves` profile resolves exactly three internal packages:

- `axiolid-core` — points, vectors, frames, scalar values, and affine transforms;
- `axiolid-linear` — line/polyline value types required transitively by curves;
- `axiolid-curve` — 2D line, circle, ellipse, polyline, and B-spline vocabulary.

Depend directly on `axiolid-core` and `axiolid-curve` with default features
disabled. `axiolid-linear` is their measured transitive dependency; applications
only need to name it directly when importing its public types.

The executable fixture at
`tests/consumers/2d-curves` converts source coordinates into the application's
canonical unit, applies a `Transform2`, and constructs a `Curve2`. Unit
conversion is scalar preprocessing: Axiolid deliberately does not impose a
global model unit or add another package for it.

## What stays out

The profile fails if it acquires any other Axiolid package. In particular, it
excludes:

- `axiolid-mesh`, `axiolid-spatial`, and mesh operation contracts/providers;
- `axiolid-model`, `axiolid-topology`, and `axiolid-brep`;
- `axiolid-construct`, `axiolid-reference`, and `axiolid-nurbs`;
- the `axiolid` application facade and CPU/GPU provider bundles.

This is a dependency-graph guarantee, not merely a feature convention. CI
resolves the isolated fixture as its own workspace, compares the complete
resolved package set with `docs/architecture/closure-profiles.toml`, compiles and
runs it, and mutation-tests the closure gate.

## Choose a larger profile only for a larger job

| Need | Profile | Internal packages |
| --- | --- | ---: |
| 2D values, affine transforms, and curve vocabulary | `2d-curves` | 3 |
| Certified 2D line/segment intersections | `linear-intersection-minimal` | 5 |
| Analytic curve/surface evaluation | `parametric-curves` | 7 |
| Exact CAD topology, B-rep, and NURBS | `cad-exact` | 11 |
| Portable reference workflows through the application facade | `rust-facade-application` | 28 |

The generated [closure profile reference](./closure-profiles.md) lists every
resolved and forbidden package. Run the same proof locally with:

```bash
cargo xtask architecture closure explain 2d-curves
cargo xtask architecture closure check
bash scripts/probe_closure_gate.sh
```

A downstream adapter such as `ifc-geometry` should therefore feature-gate its
solid path and bind its plan-view path to the leaf packages above. Pulling the
application facade for 2D representation selection would intentionally select
the broad provider bundle and defeat this closure.
