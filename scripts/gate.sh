#!/usr/bin/env bash
# Full gate for Axiolid. Trusts exit codes, not parsed output.
set -uo pipefail
export PATH="$HOME/.cargo/bin:$PATH"
cd "$(dirname "$0")/.."
gate_out="$(mktemp "${TMPDIR:-/tmp}/axiolid-gate.XXXXXX")" || exit 1
trap 'rm -f "$gate_out"' EXIT
fail=0
step() { local name="$1"; shift; printf '%-46s' "$name"; if "$@" >"$gate_out" 2>&1; then echo ok; else echo "FAIL (exit $?)"; tail -25 "$gate_out" | sed 's/^/    /'; fail=1; fi; }
step "fmt --check" cargo fmt --all -- --check
step "architecture" cargo xtask architecture check
step "architecture mutation probe" scripts/probe_layering_gate.sh
step "naming mutation probe" scripts/probe_naming_gate.sh
step "isolated build mutation probe" scripts/probe_isolated_build_gate.sh
step "closure check" cargo xtask architecture closure check
step "closure mutation probe" scripts/probe_closure_gate.sh
step "roadmap freshness" python3 scripts/check-roadmap-freshness.py
step "integration contract" scripts/check-integration-contract.sh
step "integration contract mutation" python3 scripts/probe_integration_contract_gate.py
step "downstream probe unit tests" python3 -m unittest tests/downstream/test_downstream_consumers.py
step "black-box downstream consumers" python3 scripts/test-downstream-consumers.py
step "build --workspace" cargo build --workspace
step "test --workspace" cargo test --workspace
step "test --all-features" cargo test --workspace --all-features
step "clippy" cargo clippy --workspace --all-targets -- -D warnings
step "clippy --all-features" cargo clippy --workspace --all-targets --all-features -- -D warnings
step "doc" env RUSTDOCFLAGS="-D warnings" cargo doc --workspace --no-deps
step "feature matrix" scripts/geometry-feature-matrix.sh
step "Rust facade consumer" cargo run --quiet --manifest-path tests/consumers/rust-facade-application/Cargo.toml
step "C ABI header and smoke" scripts/check-capi.sh
step "native CMake/package integration" scripts/check-native-packaging.sh
step "release script tests" python3 -m unittest scripts.test_release_scripts
step "release publish plan" python3 scripts/publish-workspace.py
step "release package preflight" python3 scripts/verify-packages.py
# Derived from `cargo metadata`, never hand-maintained: a new publishable
# member joins this loop the moment it exists (kernel#38). A crate that
# escapes the check is a crate whose standalone manifest nobody proved.
isolated_targets="$(python3 scripts/isolated-build-targets.py)" || { echo "cannot derive isolated build targets"; exit 1; }
[ -n "$isolated_targets" ] || { echo "isolated build target list is empty"; exit 1; }
for c in $isolated_targets; do
  step "isolated build -p $c" cargo build -p "$c"
done
echo
[ "$fail" -eq 0 ] && echo "GATE PASSED" || echo "GATE FAILED"
exit "$fail"
