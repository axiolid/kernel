# axiolid-backend-cpu instructions

Purpose: Portable/runtime-specialized CPU execution context.

Allowed internal dependencies: `axiolid-contracts` and operation contracts plus L2 algorithms. Follow parent `../AGENTS.md`. Do not read
`PLAN.md` unless assigned implementation or roadmap work.

## Module ownership

features.rs; topology.rs; config.rs; execution.rs. Split a module before unrelated data, validation, and algorithms grow
together. Add no empty placeholder files.

## Invariants

This crate is an execution **context** (ISA detection, worker pool, policy). It
is explicitly **not** the correctness oracle: per `docs/adr/0012` the scalar
reference implementation is owned by `axiolid-reference`, and the scalar
implementation of an operation lands before any optimized implementation of it.

Portable path is the differential oracle's target, not its owner. SIMD requires runtime detection. Optional Rayon uses a
local bounded pool. Operation providers compose this context and implement a
capability trait only when the algorithm works. Feature-gated tests must prove
default scalar selection, SIMD runtime selection, disabled-parallel rejection,
and configured local-pool worker counts. Never compile the whole workspace for
the build host only.

Public values derive `Debug` and `Clone`; add other standard traits only when
semantically valid. Validate unsupported and unavailable paths in tests.

## topology.rs — tuning inputs, not capability gates

`CpuFeatures` answers *what can this machine execute* (a correctness
question: running AVX-512 where it is absent faults). `CpuTopology` answers
*what shape is this machine* -- cache sizes, line size, sharing, logical CPU
count, heterogeneous cores. Getting topology wrong costs speed, never
correctness, so the two must not be merged.

Every field is `Option`. An undetectable value stays `None` rather than
defaulting: a fabricated 32 KiB L1 would silently mistune a strategy choice,
and `None` lets a caller apply its own policy explicitly. `heterogeneous()`
returns `Option<bool>` for the same reason -- on this Xeon it is `None`
(neither signal present), which is *undetermined*, not "homogeneous".

Detection uses sysfs only: no dependency, no `unsafe`, no CPUID, and it
works identically on x86_64 and aarch64 Linux. `logical_cpus` comes from
`available_parallelism`, so cgroup limits and CPU affinity are honoured --
a container pinned to 2 cores must not be told it has 20.

Intended use is **strategy selection**: pick the algorithm whose working set
fits the measured cache, rather than assuming one path is universally best.
`fits_in_l1/l2/llc(bytes)` exist for exactly that query. A named machine
preset may only ever *narrow* what detection reports, never widen it, so a
preset written for one machine cannot claim capabilities on another.

Run `cargo run --release -p axiolid-backend-cpu --example host_profile` to
print what the current host actually reports.
