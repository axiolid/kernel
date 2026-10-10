#!/usr/bin/env bash
# Mutation probe: can `cargo xtask gaps check` actually fail?
# Each mutation breaks one invariant the ledger promises; the check must
# reject every one, and accept the untouched file.
#
# A mutation whose anchor is missing or not unique fails loudly instead of
# leaving the ledger untouched. `--check-anchors` only verifies the anchors,
# edits nothing and runs no check (scripts/check-probe-anchors.sh, #282).
set -uo pipefail
cd "$(dirname "$0")/.." || exit 1
LEDGER=docs/architecture/capability-ledger.toml
CHECK_ANCHORS=0
[ "${1:-}" = --check-anchors ] && CHECK_ANCHORS=1
if [ "$CHECK_ANCHORS" = 0 ]; then
  BAK="$(mktemp "${TMPDIR:-/tmp}/gapsprobe.XXXXXX")" || exit 1
  cp "$LEDGER" "$BAK"
  restore() { cp "$BAK" "$LEDGER"; }
  trap 'restore; rm -f "$BAK"' EXIT
fi
fail=0
anchor_ok=1
check() {
  if [ "$CHECK_ANCHORS" = 1 ]; then anchor_ok=1; return; fi
  printf "  %-58s" "$1"
  if [ "$anchor_ok" = 0 ]; then
    # The mutation did not land, so the check's verdict would say nothing.
    echo "ANCHOR (mutation not applied)"; fail=1; anchor_ok=1; restore; return
  fi
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
# mutate <anchor> <replacement> [first]
# The anchor must occur exactly once. `first` is for a mutation that breaks
# its rule wherever it lands (any row's priority, level, issue key): it then
# needs at least one occurrence and edits the first, so routine ledger edits
# do not stale it.
mutate() {
  python3 - "$LEDGER" "$1" "$2" "$CHECK_ANCHORS" "${3:-unique}" <<'EOF' || { anchor_ok=0; fail=1; }
import sys, pathlib
p = pathlib.Path(sys.argv[1]); t = p.read_text()
n = t.count(sys.argv[2])
if n == 0 or (n > 1 and sys.argv[5] != "first"):
    where = "not found" if n == 0 else f"found {n} times, not unique"
    sys.exit(f"ANCHOR {where} in {p}: {sys.argv[2]!r}")
if sys.argv[4] != "1":
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
mutate 'issue = "straight-skeleton"' 'issue = "no-such-issue"' first
check "row pointing at an undeclared issue" 1 'is not declared in [[issue]]'
mutate 'blocked_by = ["delaunay-3d"]' 'blocked_by = ["nope"]' first
check "blocked_by naming no issue" 1 'is not another issue key'
mutate 'priority = "Low"' 'priority = "Soon"' first
check "priority that is not a board option" 1 'is not a board option'
mutate 'level = "narrow"' 'level = "mostly"' first
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
if [ "$CHECK_ANCHORS" = 1 ]; then
  if [ "$fail" = 0 ]; then echo "anchors ok"; else echo "gaps probe: stale anchors"; fi
  exit "$fail"
fi
if [ "$fail" = 0 ]; then echo "gaps probe: all mutations rejected"; else echo "gaps probe: FAILED"; fi
exit "$fail"
