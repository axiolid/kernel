# 0089 — Axiolid is dual-licensed with a commercial option

- **Status:** Accepted
- **Date:** 2026-10-28
- **Deciders:** Friedrich Schrödter
- **Supersedes:** —

## Context

Axiolid has been MPL-2.0 only. MPL-2.0 already lets a proprietary product
link the kernel, as long as the product points to the source of the
Axiolid files it ships. A product that ships *modified* Axiolid files has
to publish them. The copyright holder wants to offer such users a second
option: a commercial licence by written agreement.

The holder can only offer that for code it may license that way. Apart
from Dependabot's lockfile commits, every commit on `main` comes from the
holder. The exception is the mesh boolean absorbed into
`crates/providers/mesh/boolmesh/src/csg` (ADR 0047): that code is MPL-2.0
with a third-party copyright holder. External contributions accepted
later would otherwise come in under MPL-2.0 only and close the option
again.

## Decision

Every crate's licence is the SPDX expression
`MPL-2.0 OR LicenseRef-Commercial`, set once in
`[workspace.package]`. `axiolid-mesh-boolean-boolmesh` overrides it with
`MPL-2.0`.

- `LICENSE` stays the MPL-2.0 text. `LICENSE-COMMERCIAL` defines the
  `LicenseRef-Commercial` reference. It says a commercial licence is
  available from the copyright holder by signed written agreement, grants
  nothing itself, and excludes third-party code. It names no licensee:
  each grant is its own agreement.
- Contributions come in under contributor terms in
  `docs/guide/licensing.md`. The contributor keeps the copyright, licenses
  the contribution under MPL-2.0, and grants the holder a licence to use
  and sublicense it under other terms. Accepted contributions stay
  available under MPL-2.0.
- The native release bundle ships both licence files.
- No dedicated release: each crate carries the new metadata to crates.io
  at its next normal release (#358).

## Alternatives considered

| Option | Why not |
| --- | --- |
| Stay MPL-2.0 only | Linking is fine, but a user who ships modified files has no way to avoid publishing them. |
| Relicense to a permissive licence | Gives up the requirement that distributed kernel improvements stay open. |
| Name the licensee in the repository | A grant is a private agreement; the repository only has to say that the option exists. |
| A separate signed CLA (bot, signatures) | More process than a single-maintainer project needs today. The terms apply on submission and can become a signed CLA if outside contributions grow. |
| Release every crate now | Sixty releases for metadata alone; the next normal releases carry it (#358). |

## Consequences

**Positive**

- A commercial user can get terms that also cover modified kernel files.
- The expression is machine-readable for `cargo deny` and SBOM tools.

**Negative / costs**

- Contributors grant more than MPL-2.0 asks of them, which can put some
  off.
- Until each crate is released again, crates.io still shows `MPL-2.0`.

**Follow-ups / risks to watch**

- Before the first release under the new expression, check that
  crates.io accepts `LicenseRef-Commercial` (#358).
- Code absorbed from a third party later stays outside the commercial
  option unless its licence allows it; such a crate overrides `license`
  the way `axiolid-mesh-boolean-boolmesh` does.
- The contributor terms are not legal advice; have them reviewed before
  relying on them for a large outside contribution.

## Relation to existing code

`Cargo.toml` (`[workspace.package] license`),
`crates/providers/mesh/boolmesh/Cargo.toml`, `LICENSE-COMMERCIAL`,
`docs/guide/licensing.md`, `docs/guide/contributing.md`,
`scripts/package-native.py`, `scripts/verify-native-package.py`,
`native/CMakeLists.txt`, `tools/xtask/src/context.rs` (allowed root files).
