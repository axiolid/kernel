# C ABI v0.4

The generated public header is [`axiolid.h`](https://github.com/axiolid/kernel/blob/main/crates/facade/axiolid-capi/include/axiolid.h). Generate it with `cargo +1.99.1 xtask ffi header`; verify freshness with `cargo +1.99.1 xtask ffi check`.

## Compatibility contract

- **Symbols:** all v0.4 exports use `axiolid_v0_4_`. Existing symbols and layouts remain stable throughout the 0.4 line.
- **Integers:** handles are opaque 64-bit values; zero is invalid. Sizes/counts use C `size_t`. Indices are unsigned 32-bit values.
- **Coordinates:** positions are triples of finite `double` values in a right-handed Cartesian frame. Transform arrays are column-major affine 4×4 matrices.
- **Tolerance:** every geometric operation receives positive finite linear and angular tolerances. No default is hidden in the ABI.
- **Exactness:** `AxiolidGeometryKind_ExactBrep` and `AxiolidGeometryKind_TriangleMesh` are distinct. `AxiolidStatus_UnsupportedExact` is a refusal, never an implicit mesh fallback.
- **Provider selection:** `AxiolidContextConfig.provider_profile` accepts the stable `AXIOLID_PROVIDER_PORTABLE` value and rejects unknown integers; no concrete Rust provider type crosses the boundary. Capability records expose stable capability/provider identifiers. The table includes only operations callable through this ABI (`Healing`, `Measurement`, `MeshBoolean`, `MeshHiddenLines`, `MeshVisibleFaces`, and exact `Sweep`); Rust-only section and ray-query entry points are not advertised.

## Ownership and lifetime

| Value | Owner | Lifetime / release |
| --- | --- | --- |
| `AxiolidContextHandle` | caller | Created by `context_create`; released exactly once by `context_destroy`. Destroying it drops all child objects. |
| `AxiolidMeshHandle` | context | Created by import or `result_take_mesh`; released by `mesh_destroy` or context destruction. |
| `AxiolidResultHandle` | context | Created by an operation; released by `result_destroy`, consumed by `result_take_mesh`, or dropped with the context. |
| Drawing result (`AxiolidGeometryKind_HiddenLines`) | context | Created by `hidden_lines`; read with `hidden_lines_count`, `hidden_lines_copy` and `hidden_lines_evidence`; released by `result_destroy` (`result_take_mesh` returns `WrongResultKind`). |
| Visible-face result (`AxiolidGeometryKind_VisibleFaces`) | context | Created by `visible_faces`; read with `visible_faces_counts`, `visible_faces_copy` and `visible_faces_evidence`; released by `result_destroy` (`result_take_mesh` returns `WrongResultKind`). |
| Drawing arguments (`AxiolidViewMesh` array, per-triangle group ids, view, edge selection, region, limits) | caller | Borrowed only for the duration of the call. |
| Import arrays | caller | Borrowed only for the duration of the call and copied before return. |
| Export arrays | caller | Caller allocates using `mesh_counts`; Axiolid writes only within supplied capacities. |
| Error text buffer | caller | `context_last_error_message` reports required bytes including NUL, then copies into caller memory. |
| Error record | context | Replaced by the next recorded failure in that context; copied into `AxiolidErrorInfo` on query. |

Handles are globally unique within one process and context-scoped. Passing a child handle to another context returns `InvalidHandle`; stale or twice-destroyed handles never alias a later object. Never reinterpret handles as pointers.

## Errors and refusal

Every function returns `AxiolidStatus`. Operation failures also publish structured context: operation, provider identifier, tolerance, status, and a separately copied message. Call the error-info and error-message functions before another failing call on that context if the exact diagnostic matters.

Unknown integer operation values, non-finite coordinates, invalid affine transforms, invalid tolerances, null required pointers, oversized inputs, stale handles, wrong result kinds, and insufficient output buffers are rejected. For a size query, a null message buffer with zero capacity intentionally returns `BufferTooSmall` after writing the required size.

## Drawings

`axiolid_v0_4_hidden_lines` draws the edges of borrowed meshes in an orthographic view with hidden-line removal. Each `AxiolidViewMesh` names a context mesh, a caller-chosen unique `id` that every segment echoes, a role (`AXIOLID_MESH_ROLE_SOLID` or `SHEET`) and a group mode: `SINGLE`, `PER_TRIANGLE` (`groups` holds one id per triangle) or `PLANAR` (planar patches derived exactly). `AxiolidEdgeSelection` chooses boundary, silhouette, crease (only when `has_crease` is 1) and group-boundary edges and whether hidden pieces are reported. `AxiolidView` carries the frame, `AXIOLID_PROJECTION_ORTHOGRAPHIC` and an optional depth range, each bound read only when its `has_*` flag is 1. A null `region` is no crop. The limits are the caller's hard bounds for the call.

Every float must be finite, every integer enum known, every flag 0 or 1 and every reserved byte zero; otherwise the call returns `InvalidArgument` and records the error. A mesh count above `AxiolidContextConfig.max_meshes` is `LimitExceeded` before any entry is read, and an unknown mesh handle is `InvalidHandle`. A per-triangle `group_count` must equal the mesh's triangle count, checked (`InvalidArgument`, naming the mesh) before the ids are read. Frame orthonormality, the depth range and the limits are checked by the drawing contract, and its refusal is `OperationFailed`, except a capability the provider lacks or an input family it does not support (`Unsupported`) and an exhausted budget (`LimitExceeded`). No substitute drawing is ever returned.

A refusal that concerns one input mesh keeps that mesh's `id`: `axiolid_v0_4_context_last_error_mesh` writes `has_mesh` 1 and the id, so the caller can drop the mesh and retry. The id is meaningful only right after a status this call recorded; a later success does not clear it. It returns `NoError` when no error was recorded. `hidden_lines_copy` writes in the outcome's order (by mesh position, source edge, then position along the edge) and returns `BufferTooSmall`, writing nothing, when the capacity is below the count. The drawing providers verify their conformance suites once per process on the first call; a provider that fails refuses every call.

## Visible faces

`axiolid_v0_4_visible_faces` finds, for the same borrowed meshes and orthographic view as `hidden_lines`, the regions of the image where each face group is frontmost (ADR 0094). It takes no edge selection and has its own `AxiolidVisibleFaceLimits`: the six bounds are meshes, source vertices, source triangles, candidate pairs, output regions and output ring vertices. Arguments are validated, and refusals recorded with their mesh id, exactly as for `hidden_lines`; the recorded operation is `MeshVisibleFaces`. A group mode of `SINGLE` puts every triangle in group 0, `PER_TRIANGLE` takes the caller's ids and `PLANAR` derives planar patches exactly (ADR 0106).

The result is a set of flat tables read with `visible_faces_counts` (`regions`, `rings`, `points`, `triangles`) and `visible_faces_copy`, which fills four caller arrays:

- `AxiolidVisibleRegion` (48 bytes): the caller's `mesh_id`, a run of rings (`first_ring`, `ring_count`), a run of triangle indices (`first_triangle`, `triangle_count`), the face `group` and the `source`. The first ring of a region is its outer ring, counter-clockwise; the remaining rings are holes, clockwise. `AxiolidFaceSource` is `Group` (the region of one face group; `group` names it) or `Cap` (the near-plane section of a solid; `group` is 0 and there are no triangles). Triangle indices are the mesh's own, ascending.
- `AxiolidRegionRing` (16 bytes): `first_point` and `point_count` into the point table. Rings start at their lexicographically smallest point and do not repeat it.
- `points_xy`: two `f64` per point, x then y; `point_cap` counts points. The regions are snap-rounded to the binary64 grid (ADR 0100), and the evidence bounds how far any boundary moved.
- `triangles`: `u32` indices.

Regions are ordered as the contract orders them: by mesh position in the call, then source, then outer ring. `visible_faces_copy` returns `BufferTooSmall` and writes nothing unless all four capacities cover their counts; a table with zero capacity may be null. `AxiolidVisibleFaceEvidence` (136 bytes) carries the sixteen counts of the contract and `snap_displacement`. `exact_escalations` can differ between identical calls and is not part of the result.

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
