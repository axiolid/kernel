#!/usr/bin/env bash
# Gate step: every mutation probe's anchors still apply (#282).
#
# A mutant whose anchor text no longer occurs exactly once in the source
# cannot be applied, so a refactor that moves the code under a probe would
# otherwise leave the probe weaker than its count claims. This runs every
# Python probe in `--check-anchors` mode (anchors only, no tests, seconds in
# all) and the shell probes that edit by anchor, and fails when any anchor is
# missing or not unique. Re-anchor the mutant to the current source, or remove
# it with a stated reason in the probe.
set -uo pipefail
cd "$(dirname "$0")/.." || exit 1

fail=0
count=0
for probe in scripts/probe_*.py scripts/probe_gaps_gate.sh; do
  [ "$probe" = scripts/probe_lib.py ] && continue
  count=$((count + 1))
  # A probe with its own loop would ignore --check-anchors and run its tests.
  if [ "${probe%.py}" != "$probe" ] && ! grep -q '^probe_lib\.main(' "$probe"; then
    echo "FAIL $probe"
    echo "    does not hand its mutants to probe_lib.main (scripts/probe_lib.py)"
    fail=1
    continue
  fi
  case "$probe" in
    *.py) out=$(python3 "$probe" --check-anchors 2>&1) ;;
    *.sh) out=$(bash "$probe" --check-anchors 2>&1) ;;
  esac
  if [ $? -ne 0 ]; then
    echo "FAIL $probe"
    printf '%s\n' "$out" | sed 's/^/    /'
    fail=1
  fi
done

if [ "$fail" -eq 0 ]; then
  echo "probe anchors: all $count probes apply"
else
  echo "probe anchors: stale anchors found (see above)"
fi
exit "$fail"
