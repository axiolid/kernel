# Exact arithmetic benchmark (ADR 0068)

Standalone package, not a workspace member. Reproduces the numbers in
[ADR 0068](../../adr/0068-exact-number-type.md).

```bash
cd docs/research/exact-arithmetic-bench
cargo run --release
```

Workload: intersect two segments exactly, then take the sign of
`orient2d(a, b, intersection)`, over 20,000 deterministic random cases.
Five implementations run on identical inputs; the program asserts that all
signs agree, so a wrong implementation fails instead of looking fast.
