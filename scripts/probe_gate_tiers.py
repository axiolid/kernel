#!/usr/bin/env python3
"""Prove the gate's tier, meta-probe, slot and affected-crate logic can fail (#443).

Each mutant breaks one rule of `scripts/gate-lib.sh` or
`scripts/gate-affected.py`; the dry-run tests (`scripts/test_gate_tiers.py`,
`scripts/test_gate_affected.py`) must go red. `scripts/gate.sh` itself is
never edited: bash reads a running script incrementally, and the gate runs
this probe.
"""

from pathlib import Path

import probe_lib

LIB = "scripts/gate-lib.sh"
AFFECTED = "scripts/gate-affected.py"


def tiers(*tests):
    return [["python3", "-m", "unittest", *[f"scripts.test_gate_tiers.{t}" for t in tests]]]


AFFECTED_TESTS = [["python3", "-m", "unittest", "scripts.test_gate_affected"]]

MUTANTS = [
    (
        "escalation skipped: --branch stays branch after an infrastructure change",
        LIB,
        '    [ "$infra" -eq 1 ] && tier=full\n',
        "    :\n",
        tiers("TierTests.test_branch_tier_escalates_when_gate_infrastructure_changed"),
    ),
    (
        "slot not acquired: every run takes slot 1 at once",
        LIB,
        "      if flock -n 9; then\n",
        "      if true; then\n",
        tiers("SlotTests.test_a_run_waits_while_every_slot_is_taken"),
    ),
    (
        "slot not released: the lock descriptor leaks into the steps",
        LIB,
        'if "${gate_wrap[@]}" "$@" 9>&- >"$gate_out" 2>&1; then',
        'if "${gate_wrap[@]}" "$@" >"$gate_out" 2>&1; then',
        tiers("SlotTests.test_a_straggling_step_process_does_not_keep_the_slot"),
    ),
    (
        "affected set misses every reverse dependent",
        AFFECTED,
        "if name not in affected and wants & affected:",
        "if name not in affected and wants & set():",
        AFFECTED_TESTS,
    ),
    (
        "affected set misses transitive reverse dependents",
        AFFECTED,
        "                affected.add(name)\n                grew = True\n",
        "                affected.add(name)\n                grew = False\n",
        AFFECTED_TESTS,
    ),
    (
        "thread share not divided between slots",
        LIB,
        "  share=$(($(nproc) / slots))\n",
        "  share=$(nproc)\n",
        tiers("SlotTests.test_threads_are_split_between_slots_unless_preset"),
    ),
    (
        "steps not niced",
        LIB,
        "    gate_wrap=(nice -n 10)\n    command -v ionice >/dev/null 2>&1 && gate_wrap=(ionice -c3 nice -n 10)\n",
        "    gate_wrap=()\n",
        tiers("SlotTests.test_steps_run_niced_unless_disabled"),
    ),
    (
        "meta probes skipped although infrastructure changed",
        LIB,
        '  meta="$infra"\n',
        "  meta=0\n",
        tiers("TierTests.test_full_gate_runs_meta_probes_when_infrastructure_changed"),
    ),
    (
        "meta probes ignore GATE_META=1",
        LIB,
        '  [ "${GATE_META:-0}" = 1 ] && meta=1\n',
        "  :\n",
        tiers("TierTests.test_gate_meta_forces_the_probes"),
    ),
    (
        "undecidable change read as no infrastructure change",
        LIB,
        "    *) return 0 ;;\n",
        "    *) return 1 ;;\n",
        tiers("TierTests.test_undecidable_change_counts_as_infrastructure_changed"),
    ),
]

probe_lib.main(MUTANTS, tests=AFFECTED_TESTS, command=(), timeout=900, baseline=False)
