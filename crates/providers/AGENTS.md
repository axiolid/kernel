# Providers

Concrete capability implementations and optional heavyweight dependencies. Provider policy must not leak into contracts.

## Naming

A provider is named `axiolid-<domain>-<engine>`: the capability it
implements, then the engine it wraps (`axiolid-mesh-boolean-boolmesh`,
`axiolid-pointcloud-reconstruction-sdf`). The engine suffix is required,
not decorative -- `MeshBoolean` has two implementations in this
workspace, so a crate named `axiolid-mesh-boolean` would claim to be the
only one.

`cargo xtask architecture check` derives the expected name from the
package's own `role` and `domain` metadata and fails on a mismatch. See
ADR 0064.
