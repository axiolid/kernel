#!/usr/bin/env bash
# Mutation probe: prove each narrow closure gate can actually FAIL.
#
# A gate that never fails proves nothing. For every profile with forbidden
# packages we inject one forbidden dependency into its fixture and require the
# checker to reject it, then restore the fixture and require a clean pass.
set -uo pipefail
export PATH="$HOME/.cargo/bin:$PATH"
cd "$(dirname "$0")/.."
fail=0
backup="$(mktemp "${TMPDIR:-/tmp}/closure-probe.XXXXXX")"
trap 'rm -f "$backup"' EXIT

probe() {
  local profile="$1" injected="$2" path="$3"
  local manifest="tests/consumers/$profile/Cargo.toml"
  cp "$manifest" "$backup"
  printf '%s = { path = "../../../%s", default-features = false }\n' "$injected" "$path" >> "$manifest"
  printf '%-28s + %-26s ' "$profile" "$injected"
  if cargo xtask architecture closure check >/dev/null 2>&1; then
    echo "NOT DETECTED (gate is blind)"
    fail=1
  else
    echo "detected"
  fi
  cp "$backup" "$manifest"
}

probe linear-intersection-minimal axiolid-curve     crates/representations/analytic/curve
probe mesh-rule-checker           axiolid-nurbs     crates/algorithms/parametric/nurbs
probe 2d-curves                   axiolid-model     crates/representations/modeling/graph
probe parametric-curves           axiolid-reference crates/algorithms/reference
probe cad-exact                   axiolid-mesh      crates/representations/discrete/mesh
probe core-only                  axiolid-mesh      crates/representations/discrete/mesh
probe spatial-rule-checker       axiolid-mesh      crates/representations/discrete/mesh
probe linear-data                axiolid-linear-intersection crates/algorithms/query/intersection/linear
probe full                       axiolid-capi      crates/facade/axiolid-capi

printf '%-57s' "restored fixtures pass"
if cargo xtask architecture closure check >/dev/null 2>&1; then echo ok; else echo FAILED; fail=1; fi
[ "$fail" -eq 0 ] && echo "CLOSURE PROBE PASSED" || echo "CLOSURE PROBE FAILED"
exit "$fail"
