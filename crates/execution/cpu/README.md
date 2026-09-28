# axiolid-backend-cpu

A CPU execution context for Axiolid providers: runtime instruction-set
detection (`CpuFeatures`), measured cache and core topology for tuning
(`CpuTopology`), and an optional context-owned Rayon pool. The default build
is portable and single-threaded; the `simd` feature lets providers select an
instruction set at run time, and `parallel` adds a bounded local pool instead
of touching Rayon's global one. It bundles no geometry algorithm and is not
the correctness oracle: that is `axiolid-reference` (ADR 0012). Operation
providers compose this context.

```bash
cargo add axiolid-backend-cpu
```

- API documentation: [docs.rs/axiolid-backend-cpu](https://docs.rs/axiolid-backend-cpu)
- Source and issues: [github.com/axiolid/kernel](https://github.com/axiolid/kernel)

To see what the current host reports:

```bash
cargo run --release -p axiolid-backend-cpu --example host_profile
```
