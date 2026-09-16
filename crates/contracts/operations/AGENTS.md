# Operation contracts

Each child owns one portable request/result/evidence schema. Provider selection and fallback live in `crates/execution/dispatch`.

## Naming

A crate here is named `axiolid-<domain>-contract`, where `<domain>` is the
package's own `metadata.axiolid.domain` with dots replaced by hyphens. An
implementation NEVER takes the bare `axiolid-<domain>` name: most
capabilities here have more than one implementation, and the bare name
would let one claim to be the only one. Providers add an engine suffix
(`axiolid-mesh-boolean-boolmesh`).

`cargo xtask architecture check` enforces this; see ADR 0064. One
grandfathered exception (`axiolid-mesh-compile`, already published) is
listed explicitly in `tools/xtask/src/architecture/naming.rs`.
