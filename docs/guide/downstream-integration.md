# Downstream integration quickstarts

Four supported v0.4 paths consume Axiolid without checking out this repository or importing implementation-private crates. Every example below is exercised by CI against real code in this repository.

## Choose a profile

| Profile | Use when | Guide |
| --- | --- | --- |
| Narrow Rust crate | You need one small vocabulary (linear geometry, 2D curves, mesh values) and want the smallest dependency closure | [Narrow Rust crates](#narrow-rust-crates) |
| Rust facade | You want the supported portable-provider application boundary with one coherent API | [Rust facade](#rust-facade) |
| Plain C | Your consumer is C, or a language with a C FFI | [Plain C](#plain-c) |
| C++ / CMake | Your consumer is C++ built with CMake | [C++ and CMake](#c-and-cmake) |

Every profile is pinned to one immutable Git commit (or a released native archive) — never a branch or mutable tag.

## Capability and refusal matrix

| Request | Narrow Rust | Rust facade | C ABI |
| --- | --- | --- | --- |
| Mesh validation / measurement | crate-specific | `app.validate_mesh` / `app.measure_mesh` | `axiolid_v0_4_mesh_*` |
| Mesh Boolean (approximate) | `axiolid-mesh-boolean-boolmesh` | `app.boolean` / `app.subtract_many` | `axiolid_v0_4_*_boolean` |
| Ray/mesh nearest hit | `axiolid-ray-mesh` | `app.nearest_mesh_hit` | not exposed over the ABI |
| Exact profile extrusion | `axiolid-construct` | `app.extrude_profile_exact` | `axiolid_v0_4_exact_extrude_rectangle` |
| Exact Boolean | not implemented | typed refusal | `AxiolidStatus_UnsupportedExact` |
| Mesh section | `axiolid-mesh` + contract | not exposed on the facade `application` boundary | not exposed over the ABI |

A cell reading "typed refusal" or `AxiolidStatus_UnsupportedExact` means the operation is requested, recognized, and explicitly rejected — never silently downgraded to an approximation. See [downstream integration profiles](/architecture/downstream-integration) for the full versioned capability vocabulary and [C ABI v0.4](/architecture/c-abi-v0.4) for the native error/ownership contract.

## Narrow Rust crates

Pick the smallest closure that satisfies your use case. Every closure below is machine-checked; see [Closure profiles](/architecture/closure-profiles) for the complete, generated list including exact forbidden-package sets.

```toml
[dependencies]
# 2D plan geometry only: no solids, no CSG, no mesh.
axiolid-core = { git = "https://github.com/axiolid/kernel.git", rev = "<verified-commit>", default-features = false }
axiolid-curve = { git = "https://github.com/axiolid/kernel.git", rev = "<verified-commit>", default-features = false }
```

```rust
use axiolid_core::{Point2, Transform2};
use axiolid_curve::{Curve2, Line2};

let line = Curve2::Line(Line2 { origin: Point2::ZERO, direction: axiolid_core::Vec2::X });
let moved = Transform2::from_translation(axiolid_core::Vec2::new(1.0, 0.0));
```

This is the profile documented in full at [2D-only consumers](/architecture/2d-only-consumers); it resolves exactly three internal packages and no solid/CSG kernel. Other narrow profiles (`linear-intersection-minimal`, `mesh-rule-checker`, `parametric-curves`, `cad-exact`) follow the same pattern — see the [closure profiles table](/architecture/closure-profiles) to pick the smallest one that includes what you need.

The executed source for these fixtures lives at [`tests/consumers/`](https://github.com/axiolid/kernel/tree/main/tests/consumers); `cargo xtask architecture closure check` compiles and verifies every one on every push.

## Rust facade

Use the `application` feature when you want one coherent API over explicit portable providers, without assembling registries and implementation crates yourself.

```toml
[dependencies]
axiolid = { git = "https://github.com/axiolid/kernel.git", rev = "<verified-commit>", features = ["application"] }
```

```rust
use axiolid::application::ApplicationBuilder;
use axiolid::contracts::ExecutionOptions;
use axiolid::core::{BooleanOperator, Tolerance};

let app = ApplicationBuilder::new()
    .with_portable_boolean()?
    .with_portable_section()?
    .build();
let tolerance = Tolerance::new(1e-9, 1e-9).expect("finite tolerance");
let options = ExecutionOptions::new(tolerance);
// app.boolean(&subject, &tool, BooleanOperator::Difference, &options)?;
# Ok::<(), Box<dyn std::error::Error>>(())
```

Provider selection stays explicit; the facade never falls back from exact to mesh silently. `ApplicationError` preserves the requested operation, selected provider, tolerance, and typed underlying refusal.

The complete, CI-executed program is [`tests/consumers/rust-facade-application`](https://github.com/axiolid/kernel/tree/main/tests/consumers/rust-facade-application/src/main.rs); it exercises validation, measurement, Boolean, batched subtraction, ray queries, and strict exact-profile extrusion in one run.

## Plain C

Vendor `native/cmake/AxiolidFetch.cmake` or extract a verified [release archive](/architecture/native-distribution#release-archive), then link `axiolid.h` and the `axiolid_v0_4_*` ABI.

```c
#include <axiolid.h>
#include <stdio.h>

int main(void) {
  AxiolidVersion version = {0};
  axiolid_v0_4_version(&version);

  AxiolidContextConfig config = {AXIOLID_PROVIDER_PORTABLE, 8, 8, 64, 64};
  AxiolidContextHandle context = AxiolidContextHandle_INVALID;
  axiolid_v0_4_context_create(&config, &context);

  AxiolidTolerance tolerance = {1e-9, 1e-12};
  AxiolidResultHandle exact = AxiolidResultHandle_INVALID;
  axiolid_v0_4_exact_extrude_rectangle(context, 2, 3, 4, tolerance, &exact);
  /* exact.kind == AxiolidGeometryKind_ExactBrep */

  axiolid_v0_4_result_destroy(context, exact);
  axiolid_v0_4_context_destroy(context);
  return 0;
}
```

Every owned handle (`AxiolidContextHandle`, `AxiolidMeshHandle`, `AxiolidResultHandle`) must be released exactly once through its matching destructor; destroying a context drops every child object. See [C ABI v0.4](/architecture/c-abi-v0.4) for the complete ownership table, threading model, and error/refusal contract.

The full compiled, CI-executed source — including the typed `AxiolidStatus_UnsupportedExact` refusal path — is [`tests/native/cmake-consumer/main.c`](https://github.com/axiolid/kernel/tree/main/tests/native/cmake-consumer/main.c).

## C++ and CMake

The same ABI, linked through one CMake target regardless of consumption path:

```cmake
find_package(Axiolid 0.1 CONFIG REQUIRED PATHS "/opt/axiolid-native-v0.1.1-x86_64-unknown-linux-gnu")
target_link_libraries(your_target PRIVATE Axiolid::axiolid)
```

Or pin an immutable source commit directly:

```cmake
include(AxiolidFetch.cmake)
axiolid_fetch(
  GIT_REPOSITORY https://github.com/axiolid/kernel.git
  GIT_COMMIT <verified-commit>
  LINKAGE SHARED
)
```

`AXIOLID_LINKAGE` selects `SHARED` (default, recommended — it contains Rust's implementation dependencies) or `STATIC`. See [Native distribution and CMake](/architecture/native-distribution) for header discovery, runtime library placement on each platform, and the full supported-target table.

The complete compiled, CI-executed C++ source is [`tests/native/cmake-consumer/main.cpp`](https://github.com/axiolid/kernel/tree/main/tests/native/cmake-consumer/main.cpp); it calls the same functions as the C example above through `extern "C"` declarations, with no C++ runtime dependency across the boundary.

## What is verified, and how

Every code block above corresponds to real, compiled, executed source in this repository, not illustrative pseudocode:

- Narrow Rust closures: `cargo xtask architecture closure check` resolves and compiles each fixture under `tests/consumers/`, then `scripts/probe_closure_gate.sh` proves the closure gate can detect a forbidden dependency.
- Rust facade: `tests/consumers/rust-facade-application` runs on every push as part of `scripts/gate.sh`.
- Plain C / C++: `scripts/test-native-cmake.py` builds and runs both consumers against the in-tree source, an installed build, and an extracted release archive, in Debug and Release, with shared and static linkage — then `--mutations` proves that removing a required symbol, header, or package config breaks the consumer.
- Black-box compatibility: `.github/workflows/native.yml` and the `rust-consumers` CI job copy every consumer above into a clean temporary tree outside this workspace before building, on Linux, macOS, and Windows. See [ADR 0042](/adr/0042-black-box-downstream-compatibility-gate).

No example on this page claims support for an operation that lacks this kind of executable evidence.
