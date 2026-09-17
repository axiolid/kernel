#!/usr/bin/env bash
# Mutation probe: can the capabilities gate actually fail? (kernel#35)
set -uo pipefail
cd "$(dirname "$0")/.." || exit 1
PAGE=docs/capabilities.md
BAK="${TMPDIR:-/tmp}/capprobe.bak"
cp "$PAGE" "$BAK"
restore() { cp "$BAK" "$PAGE"; rm -f "$BAK"; }
trap restore EXIT
fail=0
check() {
  printf "  %-52s" "$1"
  if [ "$2" = "$3" ]; then echo ok; else echo "MISS (want=$2 got=$3)"; fail=1; fi
}

echo "=== baseline ==="
python3 scripts/check-capabilities.py >/dev/null 2>&1
check "every row resolves" 0 $?

echo "=== mutation 1: cite a crate that does not exist ==="
sed -i '0,/axiolid-inspect/s//axiolid-nonexistent-xyz/' "$PAGE"
python3 scripts/check-capabilities.py >/dev/null 2>&1
check "gate rejects an unresolvable crate" 1 $?
cp "$BAK" "$PAGE"

echo "=== mutation 2: strip a row of all evidence ==="
python3 scripts/_strip_capability_evidence.py
python3 scripts/check-capabilities.py >/dev/null 2>&1
check "gate rejects a row with no evidence" 1 $?
restore; trap - EXIT

python3 scripts/check-capabilities.py >/dev/null 2>&1
check "clean again after restore" 0 $?
[ "$fail" -eq 0 ] && echo "CAPABILITIES GATE MATRIX PASSED" || echo "CAPABILITIES GATE MATRIX FAILED"
exit "$fail"
