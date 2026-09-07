# tools/benchmark

The kernel's internal measurement system. Single-kernel benchmarks only —
cross-kernel comparison against OCCT, CGAL, Manifold, ifc-lite, and Truck lives
in the sibling `benchmarks` repository, which can point at an arbitrary kernel
checkout to compare base against head.

## Why it is here and not in `crates/`

Layer `tools`, `publish = false`. The architecture gate requires a `tools` layer
package to live under `tools/`, and keeping it unpublished keeps it out of the
50-crate release cascade entirely — measurement scaffolding is not a public
capability.

## Layout

```
src/workload.rs     deterministic inputs with derived ground truth
src/validate.rs     the correctness checks every benchmark must pass
benches/micro.rs        criterion wall-clock, for finding where time goes
benches/regression.rs   iai-callgrind instruction counts, for CI gating
examples/filter_probe.rs  which orient3d inputs actually escalate
```

## The two benchmark kinds answer different questions

**`micro.rs` (criterion, wall-clock).** Reports elapsed time and throughput.
Informative for finding hot paths, useless as a gate: frequency scaling and
noisy neighbours move elapsed time by more than most real regressions, so any
threshold loose enough to avoid false alarms catches nothing.

**`regression.rs` (iai-callgrind, instruction counts).** Counts instructions
under valgrind emulation. Independent of machine load and reproducible
run-to-run — measured bit-identical across consecutive runs here — which is what
makes automated comparison trustworthy. This is the CI gate.

Neither replaces the other. Instruction count ignores cache behaviour and memory
latency, so an optimisation that improves locality without removing instructions
will not show up in `regression.rs` at all.

## Running

```bash
cargo bench -p axiolid-benchmark --bench micro       # wall-clock
scripts/bench-regression.sh                          # instruction counts
```

`bench-regression.sh` skips cleanly when valgrind or `iai-callgrind-runner` is
absent, so a developer without them can still run the gate.

```bash
sudo apt-get install valgrind
cargo install iai-callgrind-runner --version 0.16.1   # must match Cargo.toml
```

## The rule every benchmark obeys

**A result that is not validated is not a measurement.** Three ways to look fast
while being wrong: declining the work, failing silently, and being optimised
away. Every workload therefore carries a ground truth *derived from its
construction* — never read back from the code under test — and the harness fails
rather than reporting a fast wrong answer.

This is carried from the sibling repo, where an unvalidated column once reported
a kernel returning zero volume in 0.2 ms as the fastest result in the table.

## Pitfalls found while building this

- **`[profile.bench]` inherits `strip = true` from `[profile.release]`.** A
  stripped binary silently defeats callgrind's function toggle: the run
  succeeds, reports `Collected: 0`, and prints zero instructions for every
  benchmark. That looks like a working harness. `strip = false` is now explicit
  in the root manifest with a comment saying why.
- **`#[library_benchmark]` rejects `///` doc comments** on the benchmarked
  function. Use `//`.
- **The `iai-callgrind-runner` binary version must equal the `iai-callgrind`
  library version.** The runner has no plain `--version` output — it reports a
  diagnostic instead — so do not try to parse one.
- **Guessing at a degenerate input does not produce a degenerate input.** Two
  attempts to construct a near-degenerate `orient3d` case (`z = 1e-30`, then
  `1e-17`) both measured identically to the well-separated case, because
  shrinking the offset shrinks the predicate's error bound with it. Only exact
  zero and sub-normal magnitudes escalate. `examples/filter_probe.rs` asks the
  filter directly instead of assuming; use it before adding a predicate case.

## Measured baseline (Xeon w7-3565X, 20 cores)

Recorded so a later reader can tell drift from noise. Instruction counts, not
time — reproducible on any machine.

| benchmark | instructions |
|---|---|
| `volume_box` | 8,714 |
| `volume_sphere` (1,024 samples) | 2,624,625 |
| `orient3d_filtered` | 2,002 |
| `orient3d_exact` | 7,379 |
| `point_index_build` (1,024 points) | 335,610 |

The `orient3d` pair is the interesting one: certification costs ~3.7x the
filtered path, and that gap is what any future optimisation of the predicates
has to move.

## Scenarios and scaling

`benches/scenario.rs` drives `axiolid::application::Application` -- the public
facade, not internal crates -- so a scenario measures what a consumer actually
pays: dispatch, provider selection, validation, and the geometry itself. The
wall-subtraction rows sweep opening count; the section rows sweep mesh density.
Both validate their output (signed volume against a derived ground truth,
contour count against the expected cut) so a fast wrong answer cannot look
like a win.

`tests/scaling.rs` asserts complexity rather than recording it. A benchmark
says an operation took 40us; these say the cost grows linearly with triangle
count and that a bounded radius query does not grow with cloud size. They run
in the normal test suite, so a change of complexity class fails the gate
instead of appearing as a slow drift on a chart nobody reads.

Both are mutation-checked: perturbing the query radius so it scans the whole
cloud, or changing a workload seed so runs stop being reproducible, must fail.

## Rejected: radix sort in the mesh audit

Recorded so it is not attempted twice. `callgrind_annotate` showed 68.8% of
`volume_sphere` inside the audit's `sort_unstable_by_key`. Replacing it with
an LSD radix sort over the `(low, high)` key:

| metric | comparison sort | radix sort |
|---|---|---|
| instructions | 2,512,504 | 2,037,744 (**-18.9%**) |
| D1 misses | 11,607 | 50,648 (**+336%**) |
| wall clock, 399k tris | 41.2 ms | 70.3 ms (**2.2x slower**) |

Fewer instructions, more than twice the time. The 256-bucket scatter defeats
the cache; the comparison sort's access pattern does not. **The audit is
memory-bound, not compute-bound** -- which is also why SIMD would not help it.

The lesson generalises: instruction count is a deterministic *regression*
signal, not a speed metric. Always confirm a win in wall clock and cache
counters before claiming one.


## Where the boolean path's cost actually is

Profiled `boolean_subtract` (8.25M instructions, the most expensive
benchmark here) with `callgrind_annotate`:

| region | share |
| --- | --- |
| `boolmesh` + libc + libm | 66.6% |
| allocator (`malloc`/`free`/`realloc`) | 11.3% |
| `shadows01` (single hottest fn, inside `boolmesh`) | 16.0% |
| transcendental (`acos` etc.) | 1.0% |
| **Axiolid's own crates** | **0.0% (884 instructions)** |

Two consequences for acceleration work:

1. **SIMD or rayon inside Axiolid cannot speed this path up.** There is
   almost no Axiolid code executing. Optimisation here means changing the
   provider, contributing upstream, or adding a second provider -- not
   vectorising kernel code.
2. **The grouping path exposes no parallelism on this workload.**
   `examples/grouping_probe.rs` measures it: 4/16/64/256 openings all
   collapse to ONE group, because the tools are mutually disjoint and get
   fused into a single cut. Groups run sequentially by necessity (each cuts
   the previous result), so parallel-across-groups would gain nothing here.
   Re-run the probe before assuming otherwise on a different workload.

The 11.3% allocator share is the one Axiolid-side lever visible: it comes
from per-operation `TriMesh` clones and temporaries. Reducing it needs a
measured before/after like any other change, and both wall-clock and
instruction counts, given the audit result above.
