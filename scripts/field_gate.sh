#!/usr/bin/env bash
# Full verification gate for the axiolid-field slice.
set -euo pipefail

unset CARGO_NET_OFFLINE CARGO_HOME GIT_DIR GIT_WORK_TREE 2>/dev/null || true
export RUSTUP_TOOLCHAIN=1.99.1

# Run against the tree this script lives in, not a hardcoded checkout.
# The absolute path made the gate test someone else's working copy: from a
# worktree it reported failures caused by unrelated uncommitted work, and
# would equally have reported success while the tree under test was broken.
cd "$(git -C "$(dirname "${BASH_SOURCE[0]}")" rev-parse --show-toplevel)"

echo "=== 1. fmt ==="
cargo fmt --all -- --check
echo "fmt OK"

echo "=== 2. field contract suites ==="
cargo test -p axiolid-field --all-features
cargo test -p axiolid-field-ops --all-features

echo "=== 3. field mutation matrix ==="
python3 scripts/probe_field_gate.py

echo "=== 4. architecture gate ==="
cargo xtask architecture check

echo "=== 5. workspace tests, all features ==="
cargo test --workspace --all-features

echo "=== 6. clippy ==="
cargo clippy --workspace --all-targets --all-features -- -D warnings
echo "clippy OK"

echo "=== 7. rustdoc ==="
RUSTDOCFLAGS='-D warnings' cargo doc --workspace --no-deps
echo "doc OK"

echo "=== 8. feature isolation ==="
cargo check -p axiolid-field --no-default-features
cargo check -p axiolid-field-ops --no-default-features
cargo check -p axiolid-field-ops --features navigation
cargo check -p axiolid --features field
cargo check -p axiolid --features field-ops
cargo check -p axiolid --features field-navigation
echo "feature isolation OK"

echo "=== 9. whitespace ==="
git diff --check
echo "whitespace OK"

echo "=== GATE COMPLETE ==="
