#!/usr/bin/env python3
"""Mutation probes for the axiolid-field gates.

Each probe injects one semantic defect into the crate source and asserts that
the field test suites FAIL. A probe that survives means the suite cannot
detect that defect, so the gate is decorative. Sources are restored verbatim
from an in-memory snapshot, and the restoration is verified by hash.
"""

from pathlib import Path

import probe_lib


ROOT = Path(__file__).resolve().parents[1]
OPS_SRC = Path("crates/algorithms/sampled/field/src")
VALUE_SRC = Path("crates/representations/sampled/field/src")

# (name, file, old, new) -- each must change real behaviour, not a no-op.
PROBES = [
    (
        "coverage-fabricates-occupancy",
        OPS_SRC / "sample.rs",
        "LayeredCell::with_layers(hits, Vec::new())?",
        "LayeredCell::with_layers(hits.clone(), hits.first().map(|h| vec![axiolid_core::Interval::new(h.w(), h.w() + 1.0)]).unwrap_or_default())?",
    ),
    (
        "coverage-drops-stacked-layers",
        OPS_SRC / "sample.rs",
        "evidence.surface_hits += hits.len();",
        "hits.truncate(1);\n            evidence.surface_hits += hits.len();",
    ),
    (
        "cell-merges-touching-occupancy",
        VALUE_SRC / "cell.rs",
        ".any(|pair| pair[0].end >= pair[1].start)",
        ".any(|pair| pair[0].end > pair[1].start)",
    ),
    (
        "config-ignores-cell-budget",
        VALUE_SRC / "config.rs",
        "if cells > budget.max_cells {",
        "if false && cells > budget.max_cells {",
    ),
    (
        "config-accepts-left-handed-frame",
        VALUE_SRC / "config.rs",
        "&& frame.x.cross(frame.y).dot(frame.z) > 0.0",
        "&& frame.x.cross(frame.y).dot(frame.z) != 0.0",
    ),
    (
        "morphology-radius-becomes-cell-count",
        OPS_SRC / "morphology.rs",
        "let reach = radius / config.cell_size();",
        "let reach = radius;",
    ),
    (
        "clearance-includes-reference-surface",
        OPS_SRC / "clearance.rs",
        "*value > w + linear",
        "*value > w - 1.0e9",
    ),
    (
        "navigation-ignores-max-step",
        OPS_SRC / "navigate.rs",
        "if rise > envelope.max_step + linear {",
        "if false && rise > envelope.max_step + linear {",
    ),
    (
        "navigation-ignores-agent-height",
        OPS_SRC / "navigate.rs",
        "if report.distance + linear < envelope.agent_height {",
        "if false && report.distance + linear < envelope.agent_height {",
    ),
    (
        "occupancy-accepts-unbalanced-crossings",
        VALUE_SRC / "cell.rs",
        "if self.surfaces.len() % 2 != 0 {\n            return Err(LayeredFieldError::UnbalancedCrossings);\n        }",
        "if false {\n            return Err(LayeredFieldError::UnbalancedCrossings);\n        }",
    ),
]

TEST_CMD = [
    "cargo", "+1.88.0", "test",
    "-p", "axiolid-field", "-p", "axiolid-field-ops",
    "--all-features", "--quiet",
]

probe_lib.main(PROBES, tests=[TEST_CMD], command=(), timeout=1800, baseline=True)
