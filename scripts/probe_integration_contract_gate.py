#!/usr/bin/env python3
"""Prove the downstream-contract freshness gate detects capability drift.

The one mutant removes a promised capability identifier; the gate
(`scripts/check-integration-contract.sh`) must go red, and must be green on
the untouched tree before and after.
"""

from pathlib import Path

import probe_lib

ROOT = Path(__file__).resolve().parent.parent
SOURCE = ROOT / "crates/contracts/common/base/src/capability_id.rs"
CHECK = ROOT / "scripts/check-integration-contract.sh"

MUTANTS = [
    ("promised capability identifier removed", SOURCE, "        MESH_BOOLEAN,\n", ""),
]

probe_lib.main(MUTANTS, tests=[[str(CHECK)]], command=(), timeout=1800, baseline=True)
