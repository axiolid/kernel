#!/usr/bin/env bash
# Mutation probe: can the semver gate actually fail? (kernel#34)
set -uo pipefail
cd "$(dirname "$0")/.." || exit 1
TARGET=crates/foundation/core/src/bounds.rs
BAK="${TMPDIR:-/tmp}/semverprobe.bak"
cp "$TARGET" "$BAK"
restore() { cp "$BAK" "$TARGET"; rm -f "$BAK"; }
trap restore EXIT
fail=0
check() {
  printf "  %-54s" "$1"
  if [ "$2" = "$3" ]; then echo ok; else echo "MISS (want=$2 got=$3)"; fail=1; fi
}

echo "=== baseline: clean tree reports no breakage ==="
python3 scripts/check-semver.py >/dev/null 2>&1
check "clean tree passes" 0 $?

echo "=== mutation: remove a public method from a published crate ==="
# Bounds3::is_finite is public API in axiolid-core, published at 0.2.0.
# Deleting it is unambiguously breaking, so the gate MUST fail.
python3 - <<'PY'
import re, pathlib
p = pathlib.Path("crates/foundation/core/src/bounds.rs")
t = p.read_text()
i = t.index('    pub fn is_finite(&self) -> bool {')
j = t.index('\n    }\n', i) + len('\n    }\n')
p.write_text(t[:i] + t[j:])
print("removed Bounds3::is_finite")
PY

python3 scripts/check-semver.py >/dev/null 2>&1
rc=$?
check "gate rejects a removed public method" 1 $rc

echo "=== restored ==="
restore
trap - EXIT
python3 scripts/check-semver.py >/dev/null 2>&1
check "clean tree passes again" 0 $?

[ "$fail" -eq 0 ] && echo "SEMVER GATE MATRIX PASSED" || echo "SEMVER GATE MATRIX FAILED"
exit "$fail"
