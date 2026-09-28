# Getting started

## Choose the narrowest dependency

Axiolid is a workspace, not a mandatory all-in-one dependency. Prefer a leaf crate when its public contract is sufficient; use the `axiolid` facade when the feature-gated composition is more convenient.

```bash
# Core scalar values, transforms, bounds, and tolerance policy.
cargo add axiolid-core

# Or: the facade with core values, meshes, and the portable CPU shell.
cargo add axiolid --features standard
```

The [crate reference](/reference/) lists every crate by layer, with its latest release and the facade feature that exposes it; [Selecting a package](/reference/selecting-packages) says which to start from.

## Start with core values

```rust
use axiolid::{Point3, Transform3};

let source = Point3::new(1.0, 2.0, 3.0);
let world = Transform3::IDENTITY.transform_point3(source);
assert_eq!(source, world);
```

`Point3`, `Vec3`, and transforms are double-precision `glam` aliases. Tolerance is explicit operation input; the kernel does not hide a global epsilon in geometry values.

## Opt into capabilities deliberately

```bash
# Mesh-oriented construction, triangulation, spatial operations, and the
# optional mesh-Boolean provider.
cargo add axiolid --features discrete

# Representation vocabulary plus general NURBS reference algorithms.
cargo add axiolid --features parametric

# Or select only curve/surface values and the general NURBS algorithms.
cargo add axiolid --features nurbs
```

The facade's default feature set is empty, so a build names what it needs. The [`axiolid` reference page](/reference/crates/axiolid#features) lists every feature and what it enables; the bundles are:

| Bundle | Includes | Does not imply |
| --- | --- | --- |
| `standard` | core values, mesh facade, portable CPU shell | every mesh algorithm |
| `discrete` | mesh-centric representations and operation contracts | a selected executable provider |
| `application` | supported portable provider selection plus v0.4 reference workflows | exact Boolean parity |
| `parametric` | curve, surface, topology, primitive, and graph vocabulary plus general NURBS reference algorithms | a complete CAD modeling/intersection kernel |
| `advanced` | `discrete` + `parametric` + healing vocabulary | production GPU computation |
| `full` | advanced facade plus optional parallel/SIMD/GPU seams | that each acceleration path is implemented or faster |

Read [Capabilities](/capabilities) before basing product behavior on a feature.

## Build and test locally

```bash
cargo fmt --all -- --check
cargo test --workspace
cargo clippy --workspace --all-targets --all-features -- -D warnings
```

For the architecture-specific feature matrix and mutation probes, see [Contributing](/guide/contributing).

## Use the supported application boundary

Enable `application` when a program wants one coherent portable provider path instead of assembling registries and implementation crates itself:

```bash
cargo add axiolid --features application
```

Provider choice is still explicit:

```rust
use axiolid::application::ApplicationBuilder;

let app = ApplicationBuilder::new()
    .with_portable_boolean()?
    .with_portable_section()
    .build();
let advertised = app.descriptor();
# Ok::<(), Box<dyn std::error::Error>>(())
```

The boundary provides mesh validation, surface/volume measurement, Boolean and batch-subtraction operations, ray/mesh queries, and strict exact-profile extrusion. `ApplicationError` preserves the requested operation, selected provider, tolerance, and typed underlying refusal. It never substitutes a mesh when exact output was requested.

The complete program is the [external Rust facade probe](https://github.com/axiolid/kernel/tree/main/tests/consumers/rust-facade-application). The closure budget is part of `cargo xtask architecture closure check`.
