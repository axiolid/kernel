# 0037 — Intersection and inversion results are verified in mapped 3D

- **Status:** Accepted
- **Date:** 2026-09-03
- **Deciders:** Friedrich Schrödter
- **Supersedes:** —

## Context

Every certified query in `axiolid-nurbs` returns *native parameter* enclosures:
curve intervals, surface `(u, v)` boxes, four-parameter surface-pair boxes. The
existing tests assert on those parameters and on residual bounds computed by the
same subdivision and interval machinery that produced the result.

That is a closed loop. If the Bézier conversion, the hull bound, or the
Krawczyk step shares a wrong assumption with the test, the test agrees with the
bug. Native parameters also carry no model units — a shifted knot domain and
non-uniform rational weights make a "small" parameter box arbitrary in 3D. A
consumer cares about one thing the parameter box never states: *do these two
pieces of geometry actually meet in model space, and how far off is the answer?*

Milestone v0.3 issue #17 requires that intersection and inversion results be
checked in mapped 3D, by something independent of the implementation under
test, with adversarial fixtures, and that failures report the 3D deviation.

## Decision

We will verify certified intersection and inversion results with a separate
`axiolid-oracle` package that maps claimed parameter boxes back into model space
through the portable scalar evaluator and measures the deviation there.

- The oracle depends on `axiolid-core`, `axiolid-contracts`, `axiolid-curve`,
  `axiolid-surface`, and `axiolid-evaluate`. It deliberately does **not** depend
  on `axiolid-nurbs`, so it shares no subdivision, hull, interval, or root
  isolation code with the implementations it checks.
- It is a **falsifier, not a prover**. `closer_point_refutation` returning a hit
  is a sound refutation of a claimed global minimum. Returning nothing is not a
  proof of global minimality — the certified path remains the only source of
  proof.
- Search is a deterministic uniform closed grid with an explicit
  `SampleDensity`. There is no hidden default and no adaptive heuristic to
  reason about; a coarse scan that reports false agreement is a policy error the
  caller can see.
- Every helper returns the measured 3D deviation and the two mapped points, so a
  failing assertion reports how far off the result was.
- The oracle lives under `tools/`, not `crates/`. It is verification
  infrastructure, not a shipped kernel capability, and it is `public = false`.

## Alternatives considered

| Option | Why not |
| --- | --- |
| Assert on native parameters and residual bounds only | The closed loop above: the implementation grades its own homework, and the assertion carries no model-space meaning. |
| Put the oracle inside `axiolid-nurbs` as a test module | It would compile against the same internals it is meant to be independent of, and reviewers could not tell shared machinery from independent evidence. |
| Reuse `axiolid-reference` as the oracle | The umbrella pulls mesh, spatial, and measure into a curve/surface check, and it is itself a reference *implementation* rather than an independent evaluator of parametric claims. |
| Adaptive/Newton-refined sampling for tighter witnesses | Reintroduces exactly the machinery under test, and makes coverage impossible to reason about by inspection. |

## Consequences

**Positive**

- A certified result is now checked against an independently derived model-space
  witness, satisfying the "independent oracle" requirement of #17.
- Deviations are reported in metres, not parameter units, so a regression says
  how wrong it is.
- The oracle can refute an inflated distance claim, which is what makes the
  passing cases meaningful rather than vacuous.

**Negative / costs**

- One more workspace package, and a dev-dependency edge from `axiolid-nurbs`.
- Sampling cost grows with the product of the scanned spans; surface-pair scans
  are quartic in density and must stay at modest densities in the gate.
- A no-hit scan proves nothing on its own, so tests must be read with that
  asymmetry in mind.

**Follow-ups / risks to watch**

- Curve/curve ownership and tangency classification (#5), curve/surface and
  surface/surface curve construction (#6), and globally certified surface
  inversion (#7) each land with mapped-3D fixtures using this oracle.
- If a future certified path needs a witness the oracle cannot express, extend
  the oracle rather than weakening the independence rule.

## Relation to existing code

- `tools/oracle/src/grid.rs` — explicit span/density policy and endpoint-exact sampling.
- `tools/oracle/src/contact.rs` — curve/curve, curve/surface, surface/surface mapped deviations.
- `tools/oracle/src/distance.rs` — sound refutation of claimed global minima.
- `tools/oracle/tests/oracle.rs` — validates the oracle against closed-form values.
- `crates/algorithms/parametric/nurbs/tests/mapped_oracle.rs` — applies it to certified results.
