# Benchmarking foundation — plan

Status: in progress. Written before execution so it survives context compaction.

## Goal

Build a trustworthy measurement system BEFORE any CPU/SIMD/Rayon/GPU work.
We have architectural seams (`CpuFeatures`, `ExecutionTarget`, dispatch) but no
numbers, so we cannot identify high-value workloads, quantify regressions,
determine scaling, or validate acceleration decisions.

Explicitly NOT in this goal: implementing CPU/GPU acceleration. No git submodule.

## Verified starting state (measured, not assumed)

- Kernel has 3 pre-existing bench files, inconsistently wired:
  - `crates/algorithms/predicates/benches/predicates.rs` — 118 lines, NO
    `[[bench]]` declaration (builds only via cargo auto-discovery), no criterion,
    hand-rolled `Instant` timing. Its doc comment says "run with
    `-p axiolid-reference`" which is the WRONG crate.
  - `crates/algorithms/reference/benches/clash.rs` — 90 lines, declared, criterion.
  - `crates/providers/mesh/boolmesh/benches/subtract_many.rs` — 162 lines, declared, criterion.
- `criterion 0.5` already in root `[workspace.dependencies]`.
- **No CI workflow runs `cargo test` or `gate.sh`.** Workflows are docs, native,
  publish, roadmap, issue-triage only. Correctness CI is a prerequisite gap for
  a *performance* gate.
- `gate.sh` runs clippy `--all-targets -D warnings` → benches must be lint-clean.
- valgrind NOT installed (needed for iai-callgrind). Installable; sudo works.
- Architecture gate (`tools/xtask`) requires EVERY workspace member to declare
  `[package.metadata.axiolid]`. Layer `tools` must live under `tools/` and must
  be `public = false`. This DECIDES the location of the bench crate.
- `axiolid-fixtures` is dev-dependency-only (used by boolmesh dev-deps).
- Sibling repo `../benchmarks` is mature: cross-kernel comparison vs boolmesh,
  ifc-lite, Manifold, CGAL, with volume validation, mutation-tested verifier,
  determinism probe. NOT a workspace member by design (absolute/relative paths
  to an arbitrary kernel checkout). Keep it that way.

## Design

### Location: `tools/benchmark/` (kernel repo)

Layer `tools`, `publish = false`, so:
- architecture gate accepts it (path prefix `tools/` matches layer `tools`),
- it never enters the crates.io publish plan (50-crate release unaffected),
- it may depend on any layer (`"tools" => true` in the layering matrix).

Rejected `crates/benchmark/`: layer rules force a `crates/<layer>/` path to a
non-tools layer, which would make the bench crate publishable and put it in the
release cascade.

### Contents

```
tools/benchmark/
  src/lib.rs            shared utilities: workload generation, validation
  benches/              criterion microbenchmarks + iai-callgrind regression
  data/                 small deterministic synthetic corpora (in-repo, tiny)
```

Key principle carried from `../benchmarks/AGENTS.md` (learned the hard way):
**a kernel that declines to answer is not faster than one that answers.** Every
timed result must be validated against a derived ground truth, and the harness
must fail on mismatch rather than report a fast wrong answer.

### Split of responsibility

- **kernel `tools/benchmark/`** — internal, single-kernel: microbenchmarks,
  deterministic regression counts, scaling curves, end-to-end scenarios.
- **`../benchmarks` (unchanged role)** — cross-kernel comparison, larger corpora,
  base-vs-head of arbitrary kernel checkouts.

## Workstreams

1. Scaffold `tools/benchmark` crate + architecture metadata; gate stays green.
2. Shared utilities: deterministic workload generators (seeded, no HashMap
   iteration order), validation helpers.
3. Criterion microbenchmarks over real hot paths.
4. iai-callgrind deterministic instruction-count benchmarks (needs valgrind).
5. Scaling tests (complexity/growth, not just single-point timings).
6. Correctness validation wired into every benchmark.
7. CI: correctness workflow first, then a performance gate.
8. Consolidate the 3 orphaned bench files into the new area.

## Validation strategy

- `gate.sh` green after every step (clippy `--all-targets` covers benches).
- `cargo bench --no-run` builds all benches.
- iai-callgrind: instruction counts stable across two consecutive runs.
- Mutation-check the validators: perturb expected values, confirm benches fail.

## Risks

- iai-callgrind needs valgrind on every machine that runs the gate → make it
  opt-in / skipped-when-absent, never a hard gate failure on a dev box.
- Criterion wall-clock in CI is noisy → wall-clock benches inform, iai-callgrind
  gates. Do not gate on wall-clock.
- Adding a workspace member changes `cargo metadata` → must re-run publish plan
  and verify the 50-crate release list is unchanged (expect exactly 50 still).
