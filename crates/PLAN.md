# Geometry package plan

Status: active after ADR 0035
Last updated: 2026-09-01

Read `AGENTS.md` for standing ownership rules. Canonical implemented structure is in [ADR 0035](../docs/adr/0035-nested-ownership-and-capability-contracts.md) and the generated [crate map](../docs/architecture/crate-map.md).

## Standing invariants

- Nested ownership folders and 32 explicit Cargo packages.
- Foundation/representation packages stay independently consumable.
- Stable, typed operation contracts for tessellation, mesh Boolean, mesh section, and graph-to-mesh.
- Exact B-rep and mesh result domains remain separate; `MeshCompiler` is explicitly mesh-valued.
- `axiolid-dispatch` owns provider ordering/fallback; contracts never select providers.
- `axiolid-reference` is the portable oracle; adopted implementations are providers.
- `axiolid-field` values and `axiolid-field-ops` algorithms are separate.
- `openbim.geometry` claim mappings remain external; Rust contracts contain no Pkl runtime or schema types.

## Shape of the work

1. Add a portable provider only after its typed contract, evidence, refusal, and conformance behavior are executable.
2. Expand exact operations without silently tessellating exact intent.
3. Add measured SIMD/parallel/GPU providers with differential reference tests.
4. Add bounded execution caches and workload diagnostics without leaking plans into contracts.
5. Extend downstream format adapters while keeping all source-format interpretation outside Axiolid.

## Required gates

```bash
cargo +1.88.0 xtask architecture check
scripts/probe_layering_gate.sh
scripts/probe_boolean_contract.py
scripts/geometry-feature-matrix.sh
scripts/gate.sh
cd docs && npm ci --ignore-scripts && npm run docs:build
```

Performance claims require repeatable benchmarks. Capability claims require implementation, diagnostics, and conformance evidence; scaffolded schemas are not capabilities.
