# axiolid-backend-gpu

An API-neutral seam for GPU graph compilation. `GpuCompiler` adapts any
`GpuGraphExecutor` to the `MeshCompiler` contract: it validates device,
precision, residency and roots before submitting one batch, and checks the
executor's results afterwards. The crate chooses no GPU API (CUDA, Metal,
Vulkan, WebGPU) and ships no executor, so default builds carry no driver
stack. Concrete executors live in their own crates, including out of tree
(ADR 0011).

```bash
cargo add axiolid-backend-gpu
```

- API documentation: [docs.rs/axiolid-backend-gpu](https://docs.rs/axiolid-backend-gpu)
- Source and issues: [github.com/axiolid/kernel](https://github.com/axiolid/kernel)

## Design notes

- Each further GPU operation gets its own narrow executor trait and
  adapter, never a method on one catch-all backend.
- A GPU path is evidence only if a CPU differential test checks it. Maturing
  the executor is tracked in
  [#22](https://github.com/axiolid/kernel/issues/22).
