#!/usr/bin/env bash
# Mutation probe: can `cargo xtask context check` actually fail?
# Each mutation breaks one rule of ADR 0078; the check must reject every one,
# for the right reason, and accept the untouched tree.
set -uo pipefail
cd "$(dirname "$0")/.." || exit 1
cargo build --quiet -p xtask || exit 1
fail=0
created=()
BAK="${TMPDIR:-/tmp}/contextprobe.$$"
mkdir -p "$BAK"
cleanup() {
  for f in "${created[@]}"; do rm -f "$f"; rmdir -p "$(dirname "$f")" 2>/dev/null; done
  for f in "$BAK"/*.path; do [ -e "$f" ] && cp "${f%.path}" "$(cat "$f")"; done
  rm -rf "$BAK"
}
trap cleanup EXIT
check() {
  printf "  %-58s" "$1"
  local out got
  out=$(cargo xtask context check 2>&1)
  got=$?
  if [ "$2" = "$got" ] && { [ -z "${3:-}" ] || grep -qF -- "$3" <<<"$out"; }; then
    echo ok
  else
    echo "MISS (want=$2${3:+ \"$3\"} got=$got)"; fail=1
  fi
  cleanup
  created=()
  mkdir -p "$BAK"
}
save() { local key; key=$(echo "$1" | tr '/' '_'); cp "$1" "$BAK/$key"; echo "$1" >"$BAK/$key.path"; }
new() { mkdir -p "$(dirname "$1")"; printf '%s\n' "$2" >"$1"; created+=("$1"); }

CRATE=crates/foundation/core
echo "=== baseline ==="
check "untouched tree passes" 0
echo "=== mutations ==="
new "$CRATE/AGENTS.md" "# nested"
check "nested AGENTS.md" 1 "only the root AGENTS.md exists"
new "$CRATE/PLAN.md" "# plan"
check "crate PLAN.md" 1 "plans are not checked in"
new "PLAN-something.md" "# plan"
check "root session plan" 1 "plans are not checked in"
new "FINDING.md" "# finding"
check "stray root-level note" 1 "the root holds only"
new "docs/plans/x.md" "# plan"
check "docs/plans page" 1 "plans are not checked in"
save "$CRATE/README.md"; rm "$CRATE/README.md"
check "crate without README" 1 "every crate has a README.md"
save "$CRATE/README.md"; printf -- '- [ ] later\n' >>"$CRATE/README.md"
check "README checkbox" 1 "task checkbox"
save "$CRATE/README.md"; seq 200 >>"$CRATE/README.md"
check "oversized README" 1 "limit 150"
save "AGENTS.md"; seq 200 >>AGENTS.md
check "oversized root AGENTS.md" 1 "limit 120"
save "$CRATE/Cargo.toml"; sed -i 's/^readme = "README.md"/readme.workspace = true/' "$CRATE/Cargo.toml"
check "publishable crate without its readme" 1 'declare `readme = "README.md"`'
save "$CRATE/src/lib.rs"; printf '// TODO: later\n' >>"$CRATE/src/lib.rs"
check "TODO without an issue" 1 "without an issue"
save "$CRATE/src/lib.rs"; printf '// TODO(#1): later\n' >>"$CRATE/src/lib.rs"
check "TODO(#N) is accepted" 0
save "$CRATE/src/lib.rs"; printf '//! See `../gone/README.md`.\n' >>"$CRATE/src/lib.rs"
check "dangling README pointer" 1 "does not resolve"
save "$CRATE/README.md"; printf 'See `crates/x/AGENTS.md`.\n' >>"$CRATE/README.md"
check "pointer to a nested AGENTS.md" 1 "points at a nested AGENTS.md"
[ "$fail" = 0 ] && echo "CONTEXT_PROBE=PASS" || echo "CONTEXT_PROBE=FAIL"
exit "$fail"
