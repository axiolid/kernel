# axiolid-benchmark

The kernel's own measurement system: deterministic workloads with a ground truth
derived from how they were built, validation that every benchmark must pass
before its number counts, wall-clock and instruction-count benchmarks, and
tests that assert how cost grows. It measures this kernel only. Comparison
against OCCT, CGAL, Manifold, ifc-lite and Truck lives in the sibling
`benchmarks` repository, which can point at any kernel checkout and compare a
base against a head. The crate is in the `tools` layer and is never published,
so it stays out of the release.

## Layout

| Path | What it answers |
| --- | --- |
| `src/workload.rs` | Seeded inputs with derived ground truth, byte-identical on every machine |
| `src/validate.rs` | The correctness checks a result must pass before it is reported |
| `src/dataset.rs` | Small inputs that once broke something, built in code with the reason beside them |
| `benches/micro.rs` | Where wall-clock time goes (criterion) |
| `benches/regression.rs` | Instruction counts under callgrind (iai-callgrind), the CI regression signal |
| `benches/scenario.rs` | End-to-end cost through the public `axiolid::application` facade |
| `benches/audit.rs` | How the mesh audit scales with triangle count |
| `benches/exact_openings.rs` | Cost per placed opening cut exactly from a placed wall, by openings already cut (#228) |
| `benches/boundary_distance.rs` | Exact boundary distance between building elements that touch or cross, and a control pair apart (#273) |
| `benches/minkowski_plan.rs` | Time and peak memory of the exact Minkowski sum and erosion of a floor's walls by a 0.9 m square, against corner count (#292) |
| `tests/scaling.rs` | Growth rate asserted from operation counts, run by `cargo test` |
| `examples/filter_probe.rs` | Which `orient3d` inputs actually escalate past the floating-point filter |
| `examples/grouping_probe.rs` | How much independent work the grouped wall subtraction exposes |

## Running

```bash
cargo bench -p axiolid-benchmark --bench micro      # wall clock
cargo bench -p axiolid-benchmark --bench scenario   # end to end
cargo bench -p axiolid-benchmark --bench exact_openings  # exact placed openings
cargo bench -p axiolid-benchmark --bench boundary_distance  # touching elements
cargo bench -p axiolid-benchmark --bench minkowski_plan  # plan morphology, time and peak memory
scripts/bench-regression.sh                         # instruction counts
cargo test -p axiolid-benchmark                     # scaling assertions
```

`scripts/bench-regression.sh` needs valgrind and an `iai-callgrind-runner`
whose version equals the `iai-callgrind` version in `Cargo.toml`. It skips
cleanly when either is missing:

```bash
sudo apt-get install valgrind
cargo install iai-callgrind-runner --version 0.16.1
```

On every pull request, `.github/workflows/performance.yml` runs the regression
benchmarks on the base commit and on the head, so the base is measured on the
same runner rather than read from a recorded baseline.

## Method

- **An unvalidated result is not a measurement.** A kernel that declines the
  work, fails silently, or is optimised away looks fast. Every case checks its
  output against the ground truth from `src/workload.rs` and fails instead of
  reporting a fast wrong answer.
- **Instruction count shows regressions, not speed.** It ignores cache
  behaviour and memory latency, so a change can remove instructions and still
  run slower. Before claiming a speedup, confirm it in wall clock and in cache
  counters (`valgrind --tool=cachegrind` or `callgrind_annotate`).
- **Wall clock is not a gate.** Frequency scaling and noisy neighbours move it
  by more than most real regressions.
- **Ask the predicate, don't guess the input.** A near-degenerate `orient3d`
  case built by shrinking an offset does not escalate, because the filter's
  error bound shrinks with it. Run `examples/filter_probe.rs` before adding a
  predicate case that claims to exercise the exact path.
- **Re-run a probe before assuming parallelism.** On the wall workload,
  `examples/grouping_probe.rs` shows 4 to 256 disjoint openings fuse into one
  group, so parallelism across groups gains nothing there.
