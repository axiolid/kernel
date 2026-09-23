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
  local out
  out=$(cargo xtask gaps check 2>&1)
  local got=$?
  # An optional third argument pins WHY it failed, so a mutation that trips
  # an unrelated rule (e.g. a TOML parse error) cannot pass as a catch.
  if [ "$2" = "$got" ] && { [ -z "${3:-}" ] || grep -qF -- "$3" <<<"$out"; }; then
    echo ok
  else
    echo "MISS (want=$2${3:+ \"$3\"} got=$got)"; fail=1
  fi
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
check "evidence symbol that is not defined" 1 'is not defined in'
mutate 'crates/algorithms/planar/overlay/src/offset.rs' 'crates/algorithms/planar/overlay/src/gone.rs'
check "evidence file that does not exist" 1 'does not exist'
mutate 'src/ModelingAlgorithms/TKGeomAlgo/IntPatch' 'src/ModelingAlgorithms/TKGeomAlgo/IntPatchX'
check "OCCT path that is not in the package list" 1 'is not a package in the pinned tree'
mutate 'issue = "straight-skeleton"' 'issue = "no-such-issue"'
check "row pointing at an undeclared issue" 1 'is not declared in [[issue]]'
mutate 'blocked_by = ["surface-intersection"]' 'blocked_by = ["nope"]'
check "blocked_by naming no issue" 1 'is not another issue key'
mutate 'priority = "Low"' 'priority = "Soon"'
check "priority that is not a board option" 1 'is not a board option'
mutate 'level = "narrow"' 'level = "mostly"'
check "level outside the four grades" 1 'unknown variant'
mutate 'id = "G2"' 'id = "G1"'
check "a row id declared twice" 1 'declared twice'
mutate 'scope_rationale = "Periodic' '# scope_rationale = "Periodic'
check "scoped row with no scope_rationale" 1 'scoped needs a non-empty scope_rationale'
mutate 'scope_rationale = "Periodic' 'scope_rationale = "   " # "Periodic'
check "scoped row with an empty scope_rationale" 1 'scoped needs a non-empty scope_rationale'
mutate 'level = "scoped"
summary = "No such structures found."' 'level = "absent"
summary = "No such structures found."'
check "scope_rationale on a row that is not scoped" 1 'only scoped rows carry one'
mutate 'evidence = []
occt = []
cgal = ["Periodic_2' 'issue = "delaunay-3d"
evidence = []
occt = []
cgal = ["Periodic_2'
check "scoped row assigned to an issue" 1 'scoped but still assigned to issue'
echo
if [ "$fail" = 0 ]; then echo "gaps probe: all mutations rejected"; else echo "gaps probe: FAILED"; fi
exit "$fail"
