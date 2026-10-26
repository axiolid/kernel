# C ABI v0.4

The generated public header is [`axiolid.h`](https://github.com/axiolid/kernel/blob/main/crates/facade/axiolid-capi/include/axiolid.h). Generate it with `cargo +1.99.1 xtask ffi header`; verify freshness with `cargo +1.99.1 xtask ffi check`.

## Compatibility contract

- **Symbols:** all v0.4 exports use `axiolid_v0_4_`. Existing symbols and layouts remain stable throughout the 0.4 line.
- **Integers:** handles are opaque 64-bit values; zero is invalid. Sizes/counts use C `size_t`. Indices are unsigned 32-bit values.
- **Coordinates:** positions are triples of finite `double` values in a right-handed Cartesian frame. Transform arrays are column-major affine 4×4 matrices.
- **Tolerance:** every geometric operation receives positive finite linear and angular tolerances. No default is hidden in the ABI.
- **Exactness:** `AxiolidGeometryKind_ExactBrep` and `AxiolidGeometryKind_TriangleMesh` are distinct. `AxiolidStatus_UnsupportedExact` is a refusal, never an implicit mesh fallback.
- **Provider selection:** `AxiolidContextConfig.provider_profile` accepts the stable `AXIOLID_PROVIDER_PORTABLE` value and rejects unknown integers; no concrete Rust provider type crosses the boundary. Capability records expose stable capability/provider identifiers. The table includes only operations callable through this ABI (`Healing`, `Measurement`, `MeshBoolean`, and exact `Sweep`); Rust-only section and ray-query entry points are not advertised.

## Ownership and lifetime

| Value | Owner | Lifetime / release |
| --- | --- | --- |
| `AxiolidContextHandle` | caller | Created by `context_create`; released exactly once by `context_destroy`. Destroying it drops all child objects. |
| `AxiolidMeshHandle` | context | Created by import or `result_take_mesh`; released by `mesh_destroy` or context destruction. |
| `AxiolidResultHandle` | context | Created by an operation; released by `result_destroy`, consumed by `result_take_mesh`, or dropped with the context. |
| Import arrays | caller | Borrowed only for the duration of the call and copied before return. |
| Export arrays | caller | Caller allocates using `mesh_counts`; Axiolid writes only within supplied capacities. |
| Error text buffer | caller | `context_last_error_message` reports required bytes including NUL, then copies into caller memory. |
| Error record | context | Replaced by the next recorded failure in that context; copied into `AxiolidErrorInfo` on query. |

Handles are globally unique within one process and context-scoped. Passing a child handle to another context returns `InvalidHandle`; stale or twice-destroyed handles never alias a later object. Never reinterpret handles as pointers.

## Errors and refusal

Every function returns `AxiolidStatus`. Operation failures also publish structured context: operation, provider identifier, tolerance, status, and a separately copied message. Call the error-info and error-message functions before another failing call on that context if the exact diagnostic matters.

Unknown integer operation values, non-finite coordinates, invalid affine transforms, invalid tolerances, null required pointers, oversized inputs, stale handles, wrong result kinds, and insufficient output buffers are rejected. For a size query, a null message buffer with zero capacity intentionally returns `BufferTooSmall` after writing the required size.

## Concurrency and foreign exceptions

Calls are safe from multiple native threads. Context lookup briefly serializes on the global handle registry; geometry work then serializes only with calls sharing that context, so independent contexts may execute concurrently. A concurrent context destruction cannot free state held by an in-flight call, but any child handle returned by that call is unreachable after destruction; callers must coordinate shutdown rather than relying on that race. Since the ABI accepts no callbacks, a C++ exception cannot enter Rust. Consumers must still catch C++ exceptions before crossing any wrapper they build around these C functions.

All exported Rust functions catch panics and return `AxiolidStatus_Panic`. No Rust unwind crosses a C frame.

## Resource budgets

`AxiolidContextConfig` sets independent hard maxima for vertices per imported mesh, triangles per imported mesh, live meshes, and live results. Values must all be non-zero. Checked multiplication precedes slice construction, so count overflow and budget excess fail before dereferencing caller memory.

## Verification

```bash
cargo +1.99.1 test -p axiolid-capi
cargo +1.99.1 xtask ffi check
scripts/check-capi.sh
```

The C smoke program compiles the generated header with C11 warnings-as-errors and covers success paths, exact refusal, invalid input, repeated destruction, explicit leak counts, and concurrent independent contexts. Rust tests additionally cover panic containment, wrong-context handles, buffer sizing, ownership transfer, Boolean and batch operations, mesh audit, bounds, measurements, transforms, and exact-result classification.
