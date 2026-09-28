# xtask

Developer automation for this workspace, run as `cargo xtask <command>` (the
alias lives in `.cargo/config.toml`). It checks the architecture metadata every
package declares, generates the architecture docs and the per-crate reference
pages of the docs site, keeps the C header in step
with the C ABI, reads the capability ledger, and guards where repository context
lives. It is never published.

| Command | What it does |
| --- | --- |
| `architecture check` | Validates `[package.metadata.axiolid]`: layers, roles, placement, exact internal dependency allowlists, naming, unsafe policy, generated-doc freshness |
| `architecture list \| graph \| docs` | Prints the model, or regenerates `docs/architecture/` from `cargo metadata` |
| `architecture closure check \| docs \| explain <profile>` | Resolves each fixture in `tests/consumers/` against `docs/architecture/closure-profiles.toml` |
| `gaps [next \| list \| show \| check]` | Reads `docs/architecture/capability-ledger.toml`; `check` keeps its evidence paths, reference paths and issue keys resolvable |
| `ffi header \| check` | Regenerates or checks `crates/facade/axiolid-capi/include/axiolid.h` |
| `docs [--check]` | Regenerates, or checks, every page derived from the crates: `docs/reference/` (one page per published crate, the index, the per-crate changelog), the sidebar facts in `docs/.vitepress/data/facts.json` and the architecture maps; checks each published README links its docs.rs and reference pages and makes no publication or version claim |
| `context check` | Enforces ADR 0078: one root `AGENTS.md`, no plan files, a `README.md` per crate, `TODO(#N)` markers only |

Every check runs in `scripts/gate.sh`. Output must be deterministic: generated
docs are compared byte for byte.

Each gate has a mutation probe in `scripts/` (`probe_layering_gate.sh`,
`probe_closure_gate.sh`, `probe_naming_gate.sh`, `probe_gaps_gate.sh`,
`probe_context_gate.sh`, `probe_docs_gate.sh`) that breaks the input and requires the check to fail.
A new rule gets a new mutation in its probe; a rule that cannot fail is not
a check.
