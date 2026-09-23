# Architecture xtask

Run `cargo xtask architecture check|list|graph|docs`. Metadata and generated docs must remain deterministic and mutation-verified.

`cargo xtask gaps [next|list|show|check]` reads `architecture/capability-ledger.toml` (`src/gaps/`). `check` is in `scripts/gate.sh`; `scripts/probe_gaps_gate.sh` proves each rule can fail. A new rule in `verify.rs` gets a new mutation in the probe.
