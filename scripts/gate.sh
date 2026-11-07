#!/usr/bin/env bash
# Gate for Axiolid. Trusts exit codes, not parsed output.
#
#   scripts/gate.sh            the full gate
#   scripts/gate.sh --branch   the checks, then build/test/clippy/doc for the
#                              crates changed since the merge base with
#                              origin/main and everything depending on them;
#                              runs the full gate instead when gate
#                              infrastructure changed (#443)
#
# One of GATE_SLOTS (default 3) machine-wide slots is held while it runs, each
# with its own persistent CARGO_TARGET_DIR. The gate's own mutation probes run
# when gate infrastructure changed or GATE_META=1 (the release flow sets it).
# Output is also written to $GATE_LOG (default <gate root>/gate-logs/).
set -uo pipefail
export PATH="$HOME/.cargo/bin:$PATH"
gate_self="$(cd "$(dirname "$0")" && pwd)/$(basename "$0")"
cd "$(dirname "$0")/.."
# shellcheck source=scripts/gate-lib.sh
. scripts/gate-lib.sh

branch_requested=0
for arg in "$@"; do
  case "$arg" in
    --branch) branch_requested=1 ;;
    *) echo "usage: scripts/gate.sh [--branch]" >&2; exit 2 ;;
  esac
done

GATE_ROOT_DIR="$(gate_root)" || { echo "cannot locate the gate root"; exit 1; }
GATE_LOCK_DIR="${GATE_LOCK_DIR:-${TMPDIR:-/tmp}/axiolid-gate-locks}"

# Re-run under tee so the log holds everything, the wait for a slot included,
# and the exit code stays the gate's own.
if [ -z "${GATE_LOGGING:-}" ]; then
  GATE_LOG="${GATE_LOG:-$GATE_ROOT_DIR/gate-logs/$(basename "$PWD")-$(date +%Y%m%dT%H%M%S).log}"
  mkdir -p "$(dirname "$GATE_LOG")" || exit 1
  echo "gate log: $GATE_LOG"
  GATE_LOGGING=1 GATE_LOG="$GATE_LOG" "$gate_self" "$@" 2>&1 | tee "$GATE_LOG"
  exit "${PIPESTATUS[0]}"
fi

# Tier and meta decision, before waiting for a slot.
gate_decide_tier "$branch_requested"
if [ "$branch_requested" -eq 1 ] && [ "$tier" = full ]; then
  echo "gate tier: full (escalated from --branch: gate infrastructure changed)"
else
  echo "gate tier: $tier"
fi
if [ "$infra" -eq 1 ]; then
  printf '%s\n' "$infra_files" | head -n 5 | sed 's/^/  infrastructure: /'
fi
if [ "$meta" -eq 1 ]; then
  echo "gate meta probes: run"
else
  echo "gate meta probes: skipped (no gate infrastructure change, GATE_META unset)"
fi

# The branch tier's crates: changed ones and everything that depends on them.
pflags=()
if [ "$tier" = branch ]; then
  affected="$(python3 scripts/gate-affected.py packages ${GATE_CHANGED_FILES_FROM:+--files-from "$GATE_CHANGED_FILES_FROM"})" \
    || { echo "cannot determine the affected crates"; exit 1; }
  for p in $affected; do pflags+=(-p "$p"); done
  # shellcheck disable=SC2086
  echo "gate crates: $(echo ${affected:-none})"
fi

gate_acquire_slot || exit 1

gate_out="$(mktemp "${TMPDIR:-/tmp}/axiolid-gate.XXXXXX")" || exit 1
trap 'rm -f "$gate_out"' EXIT
trap 'exit 130' INT
trap 'exit 143' TERM
fail=0
# First, before any build rewrites it: a stale Cargo.lock passes every
# local step and fails the release job's `cargo package --locked`.
step "Cargo.lock current" cargo metadata --locked --offline --format-version 1
step "fmt --check" cargo fmt --all -- --check
# Seconds, no tests: a refactor that moves code under a mutation probe must
# re-anchor the probe, not leave it quietly weaker (#282).
step "mutation probe anchors" scripts/check-probe-anchors.sh
step "architecture" cargo xtask architecture check
meta_step "architecture mutation probe" scripts/probe_layering_gate.sh
meta_step "naming mutation probe" scripts/probe_naming_gate.sh
meta_step "isolated build mutation probe" scripts/probe_isolated_build_gate.sh
step "closure check" cargo xtask architecture closure check
meta_step "closure mutation probe" scripts/probe_closure_gate.sh
step "context (ADR 0078)" cargo xtask context check
meta_step "context mutation probe" scripts/probe_context_gate.sh
step "generated docs" cargo xtask docs --check
meta_step "docs mutation probe" scripts/probe_docs_gate.sh
# The docs site build refuses application names (docs/scripts/check-docs-ui.mjs);
# run its source half here so a commit cannot land what the deploy rejects.
step "docs source check" node docs/scripts/check-docs-ui.mjs .
step "roadmap freshness" python3 scripts/check-roadmap-freshness.py
step "capability evidence" python3 scripts/check-capabilities.py
meta_step "capability evidence mutation probe" scripts/probe_capabilities_gate.sh
step "capability ledger (gaps)" cargo xtask gaps check
meta_step "capability ledger mutation probe" scripts/probe_gaps_gate.sh

if [ "$tier" = branch ]; then
  if [ "${#pflags[@]}" -eq 0 ]; then
    printf '%-46s%s\n' "build/test/clippy/doc" "skip (no workspace crate changed)"
  else
    step "build (affected crates)" cargo build "${pflags[@]}"
    step "test (affected crates)" cargo test "${pflags[@]}"
    step "test --all-features (affected)" cargo test "${pflags[@]}" --all-features
    step "clippy (affected)" cargo clippy "${pflags[@]}" --all-targets -- -D warnings
    step "clippy --all-features (affected)" cargo clippy "${pflags[@]}" --all-targets --all-features -- -D warnings
    step "doc (affected)" env RUSTDOCFLAGS="-D warnings" cargo doc "${pflags[@]}" --no-deps
  fi
else
  step "integration contract" scripts/check-integration-contract.sh
  meta_step "integration contract mutation" python3 scripts/probe_integration_contract_gate.py
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
  meta_step "capi soname mutation probe" scripts/probe_capi_soname.sh
  step "native CMake/package integration" scripts/check-native-packaging.sh
  step "release script tests" python3 -m unittest scripts.test_release_scripts
  step "crate release script tests" python3 -m unittest scripts.test_crate_release
  step "gate affected-crates tests" python3 -m unittest scripts.test_gate_affected
  step "gate tier and slot tests" python3 -m unittest scripts.test_gate_tiers
  meta_step "gate tier and slot mutation probe" python3 scripts/probe_gate_tiers.py
  step "release publish plan" python3 scripts/publish-workspace.py
  step "release package preflight" python3 scripts/verify-packages.py
  meta_step "minimal versions mutation probe" python3 scripts/probe_minimal_versions_gate.py
  step "semver policy" python3 scripts/check-semver.py
  meta_step "semver mutation probe" scripts/probe_semver_gate.sh
  # Derived from `cargo metadata`, never hand-maintained: a new publishable
  # member joins this loop the moment it exists (kernel#38). A crate that
  # escapes the check is a crate whose standalone manifest nobody proved.
  isolated_targets="$(python3 scripts/isolated-build-targets.py)" || { echo "cannot derive isolated build targets"; exit 1; }
  [ -n "$isolated_targets" ] || { echo "isolated build target list is empty"; exit 1; }
  for c in $isolated_targets; do
    # The facade compiles nothing without a capability feature by design
    # (kernel#9), so "builds standalone" means "builds standalone once a
    # consumer names what they want". `standard` is the migration target the
    # crate's own diagnostic recommends, so it is what gets proven here.
    if [ "$c" = "axiolid" ]; then
      step "isolated build -p $c (standard)" cargo build -p "$c" --features standard
    else
      step "isolated build -p $c" cargo build -p "$c"
    fi
  done
fi
echo
[ "$fail" -eq 0 ] && echo "GATE PASSED" || echo "GATE FAILED"
exit "$fail"
