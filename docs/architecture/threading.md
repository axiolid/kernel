# Threading

This page explains Axiolid's CPU thread-pool model, why `boolmesh`'s rayon
feature is deliberately OFF (with the measurements that decided it), and
how a caller sizes worker count. It is documentation of a shipped design,
not a proposal.

## Three separate axes

It is easy to conflate three things that are actually independent:

1. **Build-machine concurrency** -- `cargo build -j 20`. Controls how many
   `rustc` processes run in parallel while COMPILING Axiolid's crates. Has
   no effect on the resulting binary's runtime behaviour.
2. **Compile-time capability** -- the `parallel` Cargo feature on
   `axiolid-backend-cpu` (and, transitively, `axiolid-dispatch`). Decides
   whether the shipped binary CAN run multi-threaded at all: on adds
   `rayon` and the pool-construction path; off strips both.
3. **Runtime worker count** -- `CpuExecution`'s configured thread count,
   chosen per-process at startup via `CpuExecutionBuilder::threads(n)`.
   Orthogonal to (1) and (2): a `parallel`-enabled binary can still be told
   to run with 1 worker, exactly like `IfcConvert --threads 1` on a build
   that supports `-j 20`.

## Should (1) auto-toggle (2)?

No -- asked and answered explicitly during this work. "Parallel compile
turns rayon on, single-file compile turns it off" is not a coherent
policy: `cargo build -j N` is a property of the machine doing the build,
chosen fresh on every invocation, while the `parallel` feature is baked
into the artifact and shipped. Syncing them would mean a developer's
laptop core count silently decided whether a *released* binary could ever
use more than one thread -- a build-environment accident becoming a
permanent capability of the shipped crate. `parallel` is a distribution
decision (own it explicitly per release channel); worker COUNT is a
runtime decision (own it via `CpuExecutionBuilder`, mirroring how
`IfcConvert -j N` already works on the consuming side).

## `boolmesh`'s rayon feature: measured, and deliberately OFF

`crates/providers/mesh/boolmesh` does **not** enable `boolmesh`'s optional
`rayon` feature. It was enabled, benchmarked end-to-end, found to be a
net regression on realistic IFC, and reverted. This section records the
evidence so the decision is not re-litigated from the upstream headline
number alone.

There are two independent reasons the feature stays off, and only
one of them is about speed:

1. **Determinism.** The parallel path has not been audited for
   run-to-run byte stability, so `determinism()` drops to
   `BestEffort` and `Plan::admit` refuses a stronger request rather
   than returning a schedule-dependent result. This is the reason
   recorded on the feature itself in `Cargo.toml`.
2. **Measured performance.** Even setting determinism aside, it is a
   net regression on realistic IFC. That is what the rest of this
   section documents.

Upstream reports roughly 2x wall-clock with it on (8s -> 4s, Apple M4,
upstream README) via `rayon::join`/`par_iter_mut` on the internal
Morton-collider and intersection-table construction. That benchmark is a
**Menger sponge**: one enormous mesh, where per-solve work dwarfs
thread-coordination cost.

IFC is the opposite shape -- many small elements, each with modest
boolean work. Measured with real `IfcConvert` runs, same source tree and
flags, only the axiolid commit differing (best-of-N wall clock):

| fixture | subtraction tools / element | `-j1` | `-j20` |
| --- | --- | --- | --- |
| Sverchok facade, 608 elements | ~1 | **-19.5%** | **-5.1%** |
| synthetic, 8 elements x 150 openings | 150 | **-9.5%** | **-4.5%** |
| synthetic, 2 elements x 600 openings | 600 | +1.6% | +3.7% |

Negative = slower with rayon on. The crossover sits at roughly **600
subtraction tools on a single element**, and even past it the win is
~2-4% wall clock for **2.2x the CPU** (user 0.510s vs 0.231s for the same
work). Real models -- facades, MEP, sprinkler networks -- live far to the
left of that crossover, in the losing regime.

So this is not feature-flagged, it is off. A feature that is never
correct to enable on representative workloads is untested surface area,
not an option.

The lever that *did* help element-level throughput is the consumer's own
element parallelism (`IfcConvert -j N`, ~15-20%), which parallelises
across elements rather than inside one solve. If intra-solve threading is
revisited, re-run the table above before turning it back on.

## Where the scoping lives, and why

`MeshBooleanRegistry::with_execution` (in `axiolid-dispatch`, gated behind
this crate's own `parallel` feature) wraps every dispatched provider call
in a `CpuExecution`'s local rayon pool via `ThreadPool::install`, so any
provider's internal rayon work runs inside that scoped pool instead of the
ambient global one.

This is orthogonal to the decision above and remains supported: it bounds
whatever parallelism a provider does use, which is a policy question for
the embedding application. It is useful regardless of `boolmesh`'s own
feature state -- and with that feature off, it costs nothing.

This does **not** live on `BoolmeshBoolean` itself. The first
implementation attempt put `with_execution` directly on the provider; it
failed `cargo xtask architecture check` immediately:

```
axiolid-mesh-boolean-boolmesh (providers) must not depend on
axiolid-backend-cpu (execution)
```

That is not an arbitrary rule. `providers` implement a narrow contract
(`MeshBoolean`, ...) and are dispatched BY the `execution` layer
(`axiolid-dispatch`); `execution` selects among, orders, and now scopes
providers. A provider depending on execution would point the dependency
arrow backwards -- the layer that is supposed to be swappable underneath
dispatch would instead reach up into the thing dispatching it. Fixing it
at the dispatch layer instead is also strictly more useful: the SAME
`with_execution` call scopes every future `MeshBoolean` provider that gets
registered, not just `boolmesh`.

## Usage

```rust
use axiolid_backend_cpu::CpuExecutionBuilder;
use axiolid_dispatch::MeshBooleanRegistry;
use axiolid_mesh_boolean_boolmesh::BoolmeshBoolean;
use std::num::NonZeroUsize;

let execution = CpuExecutionBuilder::new()
    .threads(NonZeroUsize::new(4).unwrap())
    .build()?;

let mut registry = MeshBooleanRegistry::new().with_execution(execution);
registry.register(0, BoolmeshBoolean::new());
// Every `registry.boolean(...)` / `registry.subtract_many(...)` call now
// scopes its dispatched provider's rayon work to 4 threads.
```

Requires `axiolid-dispatch`'s `parallel` feature. Without it,
`with_execution` does not exist and dispatch is unchanged from before this
work: providers run against rayon's process-global pool, same as always.

## Measured evidence (what this does and does not fix)

Same-fixture IfcConvert benchmark (`SverchokFacadesIfcSverchok.ifc`, 608
elements, this VM: 20 cores), axiolid vs manifold vs passthrough kernels,
IfcConvert `-j 1` vs `-j 20`:

| kernel | -j 1 | -j 20 | speedup |
| --- | --- | --- | --- |
| axiolid | 0.155s | 0.138s | 1.12x |
| manifold | 0.142s | 0.115s | 1.23x |
| passthrough | 0.128s | 0.108s | 1.19x |

`IfcConvert -j N` parallelizes ELEMENT mapping (independent products
converted concurrently); it does not thread the geometry inside a single
element's boolean solve. The ~15-20% gain above is real but modest at
this fixture's scale (608 elements, sub-200ms total -- thread-pool
spin-up/sync overhead eats most of the theoretical win).

Intra-solve threading -- the complementary axis, reached via `boolmesh`'s
rayon feature -- WAS subsequently benchmarked end-to-end (see the table
under "`boolmesh`'s rayon feature" above) and turned out to be a net
regression on representative IFC, so it is off. Element-level parallelism
in the consumer is currently the only threading lever that pays.

The remaining gap to faster kernels is therefore not a threading problem:
it is per-element overhead and algorithmic cost in the single-threaded
path. That is where further optimisation work belongs.

## Addendum: inter-boolean batch threading (`parallel-batch`)

Everything above measures threading INSIDE one solve. `union_many`'s tree
reduction opened a second, independent axis: the pairs at one level of the
tree are mutually independent, so they can run concurrently without any
coordination inside a boolean. That is a different question from the one
the table above answered, and it needed its own measurement.

`parallel-batch` is therefore a SEPARATE feature from `parallel`. Enabling
one does not enable the other.

Measured on a 20-core box, disjoint grids, best-of-3, speedup against the
SAME BINARY at one thread (so codegen is held constant):

| solids | 1 | 2 | 4 | 8 | 16 |
|---|---|---|---|---|---|
| 64 | 1.00x | 1.51x | 1.95x | 2.07x | 2.07x |
| 125 | 1.00x | 1.63x | 2.18x | 2.56x | 2.59x |
| 216 | 1.00x | 1.58x | 1.94x | 2.37x | 2.28x |

**It saturates at roughly 2.1-2.6x and 16 threads buys nothing** -- at 216
solids it is slightly slower than 8. That ceiling is structural, not a
defect, and the per-level profile shows why:

```
level  unions       ms    share      (125 solids, serial)
    1      62     1.25    12.6%
    2      31     1.27    12.8%
    3      16     1.49    15.0%
    4       8     1.48    14.9%
    5       4     1.40    14.1%
    6       2     1.41    14.2%
    7       1     1.63    16.4%
```

Cost per level is nearly FLAT while the available parallel width collapses
62 -> 1: operand size doubles as operand count halves, so each level does
about the same total work. The last three levels hold 44.7% of the runtime
and can use at most 4, 2 and 1 threads.

Feeding those measured level costs into Amdahl predicts 1.72x / 2.45x /
2.95x / 3.18x at 2 / 4 / 8 / 16 threads, against 1.63x / 2.18x / 2.56x /
2.59x measured -- consistently 10-20% below the model, which is thread-pool
overhead. **The asymptotic ceiling with infinite threads is 3.28x.**

So: worth enabling for large disjoint batches, not worth expecting linear
scaling from, and not a substitute for the algorithmic win. The tree
reduction itself already bought 8.7x-18x over the sequential fold at these
sizes -- an order of magnitude more than threading adds on top.

Correctness is unaffected: `union` takes `&self` on a unit struct, the
per-level map is order-preserving (`into_par_iter().enumerate()`, never
`par_bridge`, which does not preserve order), and `tests/union_batch.rs`
gates that the batch path agrees with the sequential fold under both
feature states. `determinism()` stays `Topological`, which the general path
already was for reasons unrelated to threading.

