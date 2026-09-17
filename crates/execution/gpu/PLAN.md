# axiolid-backend-gpu plan

Design constraints for GPU execution. Status lives on GitHub, not here
(kernel#25).

## Standing invariants

The generic adapter validates device and precision policy, graph-owned
roots, and one-result-per-root cardinality before accepting executor
output. An executor that cannot satisfy those is refused rather than
trusted.

## Shape of the work

A wgpu graph compiler stays separately feature-gated, and only earns its
place with real batched compute kernels plus CPU differential tests --
a GPU path that cannot be differentially checked against the CPU one is
not evidence of anything.

Further GPU operation executors arrive as separate traits and adapters,
never as methods on one god backend.

## Exit evidence

Targeted tests, feature-isolated compile where applicable, mutation-verified
architecture/validation gates, and benchmarks before performance claims.
