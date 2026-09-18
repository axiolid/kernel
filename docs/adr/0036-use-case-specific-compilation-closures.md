# ADR 0036: Use-case-specific compilation closures

- Status: Accepted
- Date: 2026-09-02
- Supersedes: none
- Amends: [ADR 0012](./0012-scalar-reference-ownership.md)

## Context

Axiolid is a large geometry kernel. A downstream application that only needs
line queries should not compile curves, surfaces, NURBS, meshes, B-rep,
topology, providers, or execution machinery.

The pre-existing structure made this *almost* possible but not verifiable:

- `Line`/`Polyline` lived inside `axiolid-curve`, so linear values could not be
  obtained without the general curve aggregate.
- The certified `orient2d`/`orient3d` predicates lived inside the broad
  `axiolid-reference` umbrella, whose declared graph includes `axiolid-mesh`,
  `axiolid-surface`, `axiolid-spatial`, and `axiolid-measure`. Any narrow
  package needing a certified sign acquired all of it.
- The architecture checker validated declared ownership from unresolved
  metadata. It could not prove what a real consumer resolves and compiles.
- The feature matrix exercised facade features but never asserted that
  forbidden packages were **absent**.

## Decision

### 1. Hard exclusion is a package boundary, not a feature

Cargo features are additive and unified across a dependency graph. A feature
cannot promise that a package is absent, because any other dependency may
re-enable it. Only omitting the package proves omission.

### 2. Extract focused packages

- `axiolid-linear` (L1 `representations`): `Line`, `Segment`, `Ray2`,
  `Polyline`. Data only, empty default features.
- `axiolid-predicates` (L2 `algorithms`): error-free transformations, expansion
  arithmetic, `orient2d`/`orient3d`, `incircle`/`insphere`, `StaticFilter`, and
  the degeneracy scene generator. Depends only on `axiolid-core` and
  `axiolid-guarantees`.
- `axiolid-linear-intersection` (L2 `algorithms`): certified line/line and
  segment/segment classification.

### 3. Preserve every existing public path

`axiolid-curve` depends on `axiolid-linear` and re-exports it, so
`axiolid_curve::Line2` and `axiolid_linear::Line2` remain the same type.
`axiolid-reference` depends on `axiolid-predicates` and re-exports its full
surface, so `axiolid_reference::orient2d` still resolves. No consumer breaks.

### 4. Shared predicates over a literal three-package chain

The blueprint offered a literal `core → linear → linear-intersection` closure by
privately duplicating the exact arithmetic. We chose the shared-predicate
option: duplicating certified exact arithmetic would create two sources of
truth for a *sign*, which is the one value that must never disagree.

The resulting closure is five internal packages — `axiolid-core`,
`axiolid-guarantees`, `axiolid-linear`, `axiolid-predicates`,
`axiolid-linear-intersection`. Success is defined as **only relevant low-level
packages**, not a decorative package count.

### 5. Classification, never `Option<Point2>`

Linear intersection returns distinct variants for crossing, endpoint contact,
parallel-disjoint, coincident, collinear-disjoint, and overlap. Invalid input
is a typed refusal naming the operand at fault, never a plausible answer.

The parallel/coincident decision uses certified predicates that escalate to
exact arithmetic. The caller's `Tolerance` governs residual acceptance of the
computed coordinate only — it can never flip a topological branch.

### 6. Closures are machine-checked compatibility promises

`architecture/closure-profiles.toml` declares each minimal closure with its
expected and forbidden internal packages. `cargo xtask architecture closure
check` resolves an isolated consumer fixture under `tests/consumers/`
(its own workspace root, so workspace feature unification cannot mask the
result), compares the real package set, and additionally runs `cargo check`
because a resolved graph is not a compiled program.

## Amendment to ADR 0012

ADR 0012 requires that every runnable operation has an available, non-feature-
gated portable reference. That doctrine is preserved and sharpened:

> A portable reference implementation is mandatory **per operation** and lives
> in a focused, no-feature package. `axiolid-reference` is a convenience
> umbrella, not a substrate dependency. A narrow production package must depend
> on the specific reference package for the operation it implements, never on
> the umbrella.

## Consequences

### Positive

- A line-query application compiles five internal packages instead of the
  kernel.
- Certified exact arithmetic has exactly one source of truth.
- Closure regressions fail a gate instead of being discovered downstream.
- The umbrella remains available for replay, diagnostics, and broad tests.

### Negative

- The workspace has three more packages.
- Adding a package to a declared closure now requires a deliberate profile edit
  and review. This friction is the point.

### Measured

Cold builds on bbv-dev, `cargo 1.88.0`, dev profile, best of three runs with a
fresh `CARGO_TARGET_DIR` each run:

| Build | Wall time | `target/` |
|---|---:|---:|
| `axiolid` with `--no-default-features --features linear-intersection` | 3.57 s | 41 MB |
| `axiolid` with default features (`mesh`, `cpu`) | 3.64 s | 39 MB |
| `axiolid` with `--features advanced` | 7.54 s | 322 MB |

A predicate consumer's internal closure drops from 15 packages (via
`axiolid-reference`) to 3 (via `axiolid-predicates`), no longer compiling
`axiolid-curve`, `axiolid-mesh`, `axiolid-surface`, `axiolid-spatial`,
`axiolid-measure`, `axiolid-primitive`, `axiolid-contracts`, and the three mesh
operation contracts.

### Not claimed

- The `linear-intersection` and default facade builds are within noise of each
  other; the meaningful separation is against `advanced` (2.1x wall time, 7.8x
  `target/`). Package omission is not claimed to help a build that was already
  small.
- These are dev-profile numbers on one machine with a warm registry. They are
  not a release-profile or CI-cold-cache result.
- No binary-size result is claimed: the fixture is a `println!` harness, so its
  linked size measures the harness, not the geometry.
- `axiolid-linear-intersection` covers 2D line/line and segment/segment only.
  3D linear intersection and curve/surface intersection remain future work.

---

## Phase 2 amendment: archetype profiles and the evaluation split

Phase 1 proved one closure. A single verified profile does not show that the
architecture generalises, and the design document names four application
archetypes.

### The umbrella leak

`axiolid-nurbs` depended on `axiolid-reference` only for curve/surface
evaluation, but the umbrella also carries mesh, spatial, measure, primitive,
and the mesh contracts. Measured consequence: a CAD application resolving
curves, surfaces, topology, B-rep, and NURBS compiled **18** internal packages,
of which seven were discrete-geometry packages it never called.

### Decision

Extract the self-contained `curve` / `surface` / `nurbs` evaluation cluster from
`axiolid-reference` into `axiolid-evaluate` (L2, `algorithm.parametric`), and
depend on it from `axiolid-nurbs`. `axiolid-reference` re-exports `curve` and
`surface` unchanged, so `axiolid_reference::curve::evaluate3` and friends keep
resolving to the same items.

`axiolid-reference` remains the convenience oracle umbrella; it did not lose
capability, only exclusivity.

### Measured result

| Profile | Internal packages | Note |
| --- | --- | --- |
| `linear-intersection-minimal` | 5 | unchanged |
| `mesh-rule-checker` | 4 | discrete stack only |
| `parametric-curves` | 7 | evaluation without the umbrella |
| `cad-exact` | 11 | was 18 before this amendment |

The CAD closure no longer contains `axiolid-mesh`, `axiolid-spatial`,
`axiolid-measure`, `axiolid-primitive`, `axiolid-reference`, or the mesh
contract packages.

### Governance

All four profiles are verified by `cargo xtask architecture closure check`,
which resolves each fixture, compiles it, and additionally fails when
`docs/architecture/closure-profiles.md` drifts from the declarations.
`scripts/probe_closure_gate.sh` injects a forbidden dependency into every
fixture and requires the checker to reject each one, so the gate is proven
capable of failing rather than assumed to be.

### Not claimed

- No claim that four archetypes cover every construction-industry application.
  They are the archetypes the design document names; more can be added.
- No claim that `axiolid-evaluate` is a complete parametric kernel. It is the
  scalar evaluation oracle; certified NURBS inversion and general intersection
  remain future work.
- Package count is not runtime performance. It bounds what a consumer compiles,
  not how fast the resulting code runs.

---

## Phase 3 amendment: 2D plan-geometry boundary

Issue #1 asked whether a plan viewer can use Axiolid curves and transforms without
implicitly compiling the solid/CSG stack. The answer is yes, with leaf-package
dependencies rather than the application facade.

### Decision

Declare `2d-curves` as a compatibility profile. Its consumer depends directly on
`axiolid-core` and `axiolid-curve`; `axiolid-linear` is the curve vocabulary's
only transitive Axiolid package. Unit conversion remains scalar preprocessing at
the adapter boundary and therefore requires no geometry dependency.

The profile explicitly forbids every other Axiolid workspace package. This
includes `axiolid-model`, `axiolid-topology`, `axiolid-brep`, `axiolid-mesh`,
`axiolid-construct`, the application facade, execution packages, and providers.
A downstream plan viewer can therefore select 2D representations, normalize
units, apply affine transforms, and retain analytic 2D curves without acquiring
a solid kernel.

### Measured result

`cargo xtask architecture closure explain 2d-curves` resolves exactly three
internal packages and one external package (`glam`). The isolated fixture
constructs and transforms real public values. `scripts/probe_closure_gate.sh`
then injects `axiolid-model` and requires the closure check to reject the graph.

### Not claimed

- `axiolid-curve` contains both 2D and 3D curve vocabulary; this profile proves a
  small dependency graph, not dimensional erasure within that crate.
- Curve evaluation is a separate choice. Adding `axiolid-evaluate` intentionally
  selects the broader `parametric-curves` closure, including surface values.
- Full solid construction should use the application facade or the exact CAD
  packages documented by the other profiles; it is not part of `2d-curves`.

### Measured again, repeatably (2026-09-18, kernel#12)

The numbers above were a spot check on two profiles. `scripts/closure-bench.py`
now measures every declared profile, and reports a noise floor so a
difference can be called a result or not. Dev profile, 3 reps, cold
target dir per profile, cargo 1.88.0, axiolid 0.3.0:

| Profile | Packages | Cold build | `target/` |
|---|---:|---:|---:|
| `2d-curves` | 5 | 2.41 s | 47 MB |
| `mesh-rule-checker` | 6 | 2.65 s | 74 MB |
| `linear-intersection-minimal` | 7 | 2.37 s | 47 MB |
| `parametric-curves` | 16 | 3.21 s | 120 MB |
| `cad-exact` | 20 | 4.30 s | 185 MB |
| `rust-facade-application` | 64 | 6.49 s | 660 MB |
| `c-abi-profile` | 65 | 7.40 s | 856 MB |

Measured noise floor: 2.6% of a 2.39 s median.

The spread is real and large at the ends: 13x in package count and 18x
in `target/` size between the smallest profile and `c-abi-profile`.
But `2d-curves` and `linear-intersection-minimal` differ by 1.7%,
BELOW the floor -- they are indistinguishable here, not ranked. This
restates the original finding: closure profiles separate a big
consumer from a small one; they do not usefully separate two small
ones. Disk size, unlike wall time, does separate cleanly at every step.

