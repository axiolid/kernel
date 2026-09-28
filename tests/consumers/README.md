# Consumer fixtures

Each directory is a small downstream application that depends on Axiolid the
way a real consumer would. Each one backs one profile in
[`docs/architecture/closure-profiles.toml`](../docs/architecture/closure-profiles.toml)
([ADR 0036](../../docs/adr/0036-use-case-specific-compilation-closures.md)). The
profile names the internal packages that must resolve and the ones that must
not. A fixture proves that a focused application can build without pulling in
the rest of the kernel.

| Fixture | What it proves |
| --- | --- |
| `core-only` | `axiolid-core` is usable alone: points, vectors, frames, intervals, tolerances |
| `linear-data` | Alignments, centrelines and polylines can be stored and measured with no query algorithm |
| `linear-intersection-minimal` | Line queries need only core, linear, linear intersection and predicates |
| `2d-curves` | 2D plan geometry and transforms, with unit conversion left to the application, pull in no solids or CSG |
| `parametric-curves` | Curve and surface evaluation works without the reference umbrella |
| `spatial-rule-checker` | Proximity rules over points need a spatial index and core values, with no discrete geometry |
| `mesh-rule-checker` | Mesh rule checks need only mesh values, spatial acceleration and measurement |
| `cad-exact` | Analytic curves and surfaces, topology, exact B-rep and NURBS resolve together |
| `rust-facade-application` | The `axiolid` facade with portable providers runs the reference workflows |
| `c-abi-profile` | The internal Rust closure behind the C ABI. It calls only the version symbol. The C compile, link and run probe is `crates/facade/axiolid-capi/tests/c/smoke.c` |
| `full` | Every facade feature at once: the upper bound the narrow profiles are measured against |

## Rules for a fixture

- Its `Cargo.toml` has an empty `[workspace]` table, so it is its own
  workspace root. Inside the kernel workspace, feature unification would hide
  what a real consumer resolves.
- It depends on leaf packages by relative path with `default-features = false`.
  The facade and C ABI fixtures depend on `axiolid` or `axiolid-capi` with the
  features their profile names.
- It exercises real behaviour, not only type names. A fixture that only names
  a type can pass while the package is unusable.
- Adding a dependency changes a declared compatibility promise. Change the
  profile on purpose, then regenerate the closure docs.
- A fixture's `.gitignore` ignores its `Cargo.lock`. The checker resolves with
  `--offline` and writes the lock on first resolution. Build output is ignored
  by the root `**/target/` rule and never committed.

## Checks

```bash
cargo xtask architecture closure check
cargo xtask architecture closure explain <profile>
cargo xtask architecture closure docs            # after a profile change
bash scripts/probe_closure_gate.sh               # each closure rule can fail
python3 scripts/closure-bench.py --reps 3        # package count, cold build, target size
```

`closure-bench.py` measures the noise floor first by rebuilding one profile
several times, and it does not rank profiles whose gap is smaller than that
floor. Each profile builds cold in its own `CARGO_TARGET_DIR` under `--bench-root`,
outside the repository.

## Downstream gate

`scripts/test-downstream-consumers.py` copies six of these fixtures
(`linear-intersection-minimal`, `mesh-rule-checker`, `2d-curves`,
`parametric-curves`, `cad-exact`, `rust-facade-application`) into separate
temporary workspaces. It replaces each path dependency with an exact version
pinned to one immutable Git commit, then checks the source identities Cargo
resolved. `tests/downstream/test_downstream_consumers.py` unit-tests that
manifest rewriting and policy:

```bash
python3 -m unittest tests/downstream/test_downstream_consumers.py
python3 scripts/test-downstream-consumers.py
```
