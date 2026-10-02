# 0079 — Least-squares fitting lives in axiolid-nurbs and depends on axiolid-numeric

- **Status:** Accepted
- **Date:** 2026-10-02
- **Deciders:** Friedrich Schrödter
- **Supersedes:** —

## Context

Issue #150 (ledger B11) asked for least-squares B-spline curve and
surface approximation, distinct from the exact-interpolation fit
`axiolid-nurbs::fit::interpolate_curve3` already has. The solves
involved are the standard ones for this problem: an `m x n` collocation
matrix least squares (Piegl and Tiller's global curve/surface
approximation), optionally subject to linear equality constraints
(endpoint interpolation) and an optional Tikhonov-style smoothing term.
`interpolate_curve3` itself solves its square system with a hand-rolled
Gaussian elimination local to `fit.rs`.

Issue #136 (ledger A5, closed by the 0.1.0 release of `axiolid-numeric`)
built exactly this general substrate: column-pivoted QR with numerical
rank, `least_squares`, `constrained_least_squares` with condition
estimates, each refusing non-finite input and rank-deficient systems by
name rather than solving into infinities. Building a second, narrower
copy of that solver inside the fitting code would duplicate tested,
numerically-reasoned work and silently drift from it, which is the
dependency failure mode ADR 0036 already names for `axiolid-predicates`
in the `linear-intersection-minimal` closure.

`axiolid-nurbs` had never depended on `axiolid-numeric` before, so
taking the edge changes four declared minimal closures
(`docs/architecture/closure-profiles.toml`): `cad-exact`,
`rust-facade-application`, `c-abi-profile` and `full`, each of which
already carries `axiolid-nurbs`. `cargo xtask architecture closure
check` fails until `expected_internal` is updated to include
`axiolid-numeric` for all four, and the root `AGENTS.md` requires that
change to be deliberate and recorded in an ADR, not just edited past
the gate.

## Decision

Least-squares fitting (`fit_curve3`, `fit_curve3_to_tolerance`,
`fit_surface_grid`) is added to `axiolid-nurbs` as a new module
(`src/approximate.rs`), beside the existing interpolation in `fit.rs`,
and it solves through `axiolid-numeric`'s `least_squares` and
`constrained_least_squares` rather than hand-rolling normal equations
or another Gaussian elimination.

Why `axiolid-nurbs` is the right place:

- It already owns the B-spline representation's algorithms
  (differential geometry, transforms, interpolation, lofting,
  projection, intersection) and the Cox-de Boor basis evaluation
  (`span_of`, `basis_at` in `fit.rs`) that both the existing
  interpolation and the new approximation need. Approximation and
  interpolation are two outputs of the same basis-evaluation and
  knot-placement machinery; splitting them into separate crates would
  either duplicate that machinery or force one crate to depend on the
  other for no reason beyond file organisation.
- The dependency rule already allows `axiolid-nurbs` to sit in the
  `algorithms` layer depending on representation crates
  (`axiolid-curve`, `axiolid-surface`) and other algorithm crates
  (`axiolid-evaluate`, `axiolid-exact`, `axiolid-predicates`);
  `axiolid-numeric` is a peer algorithm crate with a `numeric.substrate`
  role and no dependencies of its own, so the new edge is exactly the
  "algorithms depending on algorithms" shape the layer already permits,
  not a layering violation.

Why the `axiolid-numeric` edge is acceptable to take broadly rather
than behind a feature flag:

- `axiolid-numeric` has zero dependencies (`allowed-internal-dependencies
  = []` in its own manifest) and operates on plain `f64` slices and a
  small row-major `Matrix` — no geometry types, no additional transitive
  crates, no native or optional backend. Every closure profile that
  gains it gains exactly one small, dependency-free package.
- It is the shared substrate #136 built for precisely this purpose: its
  own README says new code should use it in place of "the interpolation
  solve in axiolid-nurbs" by name. Declining the dependency here to
  avoid a closure diff would mean #136 shipped a crate that its first
  real consumer still doesn't use, and would leave the new fitting code
  as a second, divergent implementation of the same refusal behaviour
  (rank deficiency, non-finite input) that `axiolid-numeric` already
  gets right with condition estimates this crate would not reproduce.

## Alternatives considered

| Option | Why not |
| --- | --- |
| A new `axiolid-fitting` crate for the approximation code, depending on `axiolid-nurbs` and `axiolid-numeric` | Would need to re-expose `axiolid-nurbs`'s private basis/knot helpers (`span_of`, `basis_at`, `chord_parameters`, `collapse`) as a public internal API, or duplicate them, purely to satisfy a crate boundary with no other consumer in view. No downstream closure wants approximation without interpolation or vice versa, so the split buys no narrower application closure, only an extra crate to version and an extra edge in the graph (`axiolid-fitting -> axiolid-nurbs` as well as `-> axiolid-numeric`). Revisit if a consumer needs approximation without the rest of `axiolid-nurbs`. |
| Keep approximation in `axiolid-nurbs` but gate it behind a Cargo feature so closures that don't want it can exclude it | `axiolid-numeric` has no dependencies and no runtime cost when unused (dead code is simply not called), so a feature flag would trade a one-package closure growth for a second axis of conditional compilation that `cargo xtask architecture check`, the closure gate and every downstream consumer now has to account for, for a saving the dependency-free crate doesn't need. The workspace does gate real optional cost this way (`axiolid-dispatch?/parallel` in the `full` profile); this is not that case. |
| Hand-roll the least-squares and constrained least-squares solves locally in `approximate.rs`, as `fit.rs` does for its square system | Duplicates `axiolid-numeric`'s QR, rank detection and condition estimation, which #136 built and tested specifically so this kind of code would not re-derive it; a second copy drifts (different rank threshold, different refusal messages) and was explicitly what the task required against: "do not hand-roll normal equations." |

## Consequences

**Positive**

- One, already-reviewed least-squares implementation is used by both
  curve and surface fitting, with condition estimates and rank
  detection the local Gaussian elimination in `fit.rs` does not have.
- `axiolid-numeric`'s existence is now justified by a real consumer
  instead of being closure-dead weight.

**Negative / costs**

- Four closures (`cad-exact`, `rust-facade-application`,
  `c-abi-profile`, `full`) grow by one package: `axiolid-numeric`. Each
  has `expected_internal` updated with a comment in
  `docs/architecture/closure-profiles.toml` pointing at this ADR.
- A consumer that wants `axiolid-nurbs`'s transform/projection/
  intersection algorithms but needs to avoid pulling in a least-squares
  solver (for code size on a constrained target, say) can no longer get
  that narrower closure without the alternative above (a feature flag
  or a crate split) being revisited.

**Follow-ups / risks to watch**

- If a future consumer genuinely needs `axiolid-nurbs` without
  `axiolid-numeric` (a narrow embedded closure profile), reconsider the
  feature-flag alternative above rather than accepting the edge as
  unconditional.
- `axiolid-nurbs`'s own interpolation solve (`fit.rs::solve`) still
  hand-rolls Gaussian elimination; migrating it onto `axiolid-numeric`
  too (its `Lu` or `Qr`) would remove the one remaining duplicate and is
  tracked as follow-up work, not done here to keep this change to the
  new capability only.

## Relation to existing code

- `crates/algorithms/parametric/nurbs/src/approximate.rs`: the new
  module; `fit_curve3`, `fit_curve3_to_tolerance`, `fit_surface_grid`.
- `crates/algorithms/parametric/nurbs/Cargo.toml`:
  `axiolid-numeric.workspace = true` and the matching
  `allowed-internal-dependencies` entry.
- `docs/architecture/closure-profiles.toml`: `expected_internal` for
  `cad-exact`, `rust-facade-application`, `c-abi-profile` and `full`
  each gain `axiolid-numeric`, with a comment pointing here.
- `docs/architecture/capability-ledger.toml` row B11: evidence updated
  to the new functions.
