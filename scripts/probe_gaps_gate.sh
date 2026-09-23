#!/usr/bin/env bash
# Mutation probe: can `cargo xtask gaps check` actually fail?
# Each mutation breaks one invariant the ledger promises; the check must
# reject every one, and accept the untouched file.
set -uo pipefail
cd "$(dirname "$0")/.." || exit 1
LEDGER=architecture/capability-ledger.toml
BAK="${TMPDIR:-/tmp}/gapsprobe.bak"
cp "$LEDGER" "$BAK"
restore() { cp "$BAK" "$LEDGER"; }
trap 'restore; rm -f "$BAK"' EXIT
fail=0
check() {
  printf "  %-58s" "$1"
  cargo xtask gaps check >/dev/null 2>&1
  local got=$?
  if [ "$2" = "$got" ]; then echo ok; else echo "MISS (want=$2 got=$got)"; fail=1; fi
  restore
}
mutate() { python3 - "$LEDGER" "$1" "$2" <<'EOF'
import sys, pathlib
p = pathlib.Path(sys.argv[1]); t = p.read_text()
assert sys.argv[2] in t, f"probe anchor missing: {sys.argv[2]}"
p.write_text(t.replace(sys.argv[2], sys.argv[3], 1))
EOF
}

echo "=== baseline ==="
check "untouched ledger passes" 0
echo "=== mutations ==="
mutate 'offset.rs::offset_polygons' 'offset.rs::offset_polygons_gone'
check "evidence symbol that is not defined" 1
mutate 'crates/algorithms/planar/overlay/src/offset.rs' 'crates/algorithms/planar/overlay/src/gone.rs'
check "evidence file that does not exist" 1
mutate 'src/ModelingAlgorithms/TKGeomAlgo/IntPatch' 'src/ModelingAlgorithms/TKGeomAlgo/IntPatchX'
check "OCCT path that is not in the package list" 1
mutate 'issue = "straight-skeleton"' 'issue = "no-such-issue"'
check "row pointing at an undeclared issue" 1
mutate 'blocked_by = ["surface-intersection"]' 'blocked_by = ["nope"]'
check "blocked_by naming no issue" 1
mutate 'priority = "Low"' 'priority = "Soon"'
check "priority that is not a board option" 1
mutate 'level = "narrow"' 'level = "mostly"'
check "level outside implemented|narrow|absent" 1
mutate 'id = "G2"' 'id = "G1"'
check "a row id declared twice" 1
echo
if [ "$fail" = 0 ]; then echo "gaps probe: all mutations rejected"; else echo "gaps probe: FAILED"; fi
exit "$fail"
