#!/usr/bin/env python3
"""Mutation probes for tessellation and interference."""
from pathlib import Path

import probe_lib

ROOT = Path(__file__).resolve().parent.parent

SUITES = [
    ["cargo", "test", "-p", "axiolid-reference", "--all-features", "--test", "tessellate"],
    ["cargo", "test", "-p", "axiolid-reference", "--all-features", "--test", "clash"],
]

PROBES = [
    (
        "tess-ignores-the-chord-budget",
        "crates/algorithms/reference/src/tessellate.rs",
        "Some(s) if s > budget.chord_tolerance => {}",
        "Some(s) if s > budget.chord_tolerance * 1.0e12 => {}",
    ),
    (
        "tess-cell-diagonal-choice-inverted",
        "crates/algorithms/reference/src/tessellate.rs",
        "if shorter_diagonal_is_ad(&positions, a, b, c, d) {",
        "if !shorter_diagonal_is_ad(&positions, a, b, c, d) {",
    ),
    (
        "tess-sagitta-measured-to-the-start-not-the-chord",
        "crates/algorithms/reference/src/tessellate.rs",
        "let deviation = point_to_segment(pm, pa, pb);",
        "let deviation = (pm - pa).length();",
    ),
    (
        "tess-refinement-never-doubles",
        "crates/algorithms/reference/src/tessellate.rs",
        "n = (2 * (n - 1) + 1).min(budget.max_samples_per_direction);",
        "n = (n + 1).min(budget.max_samples_per_direction);",
    ),
    (
        "tess-exhaustion-is-hidden",
        "crates/algorithms/reference/src/tessellate.rs",
        "nu >= budget.max_samples_per_direction || nv >= budget.max_samples_per_direction",
        "false",
    ),
    (
        "tess-accepts-a-single-sample",
        "crates/algorithms/reference/src/tessellate.rs",
        "if max_samples_per_direction < 2 {",
        "if max_samples_per_direction < 1 {",
    ),
    # `tess-projection-unclamped` was removed after measurement rather than
    # left as a permanent leak. The sagitta probe evaluates the MIDPOINT of a
    # chord; a midpoint's projection foot always lies within its own segment,
    # so the clamp cannot be reached from `tessellate_patch`. It is defensive
    # code for future callers that pass an arbitrary point, not a live branch,
    # and a probe that no reachable input can kill is noise in the report.
]

CLASH_PROBES = [
    (
        "clash-broad-phase-never-rejects",
        "crates/algorithms/reference/src/clash.rs",
        "if !box_a.intersects(box_b) {",
        "if false {",
    ),
    (
        "clash-proper-crossing-is-not-penetration",
        "crates/algorithms/reference/src/clash.rs",
        "report.penetrating_pairs.push((i, j));\n                    report.kind = Interference::Penetrating;",
        "report.penetrating_pairs.push((i, j));\n                    report.kind = Interference::Touching;",
    ),
    (
        "clash-containment-never-checked",
        "crates/algorithms/reference/src/clash.rs",
        "if report.kind != Interference::Penetrating && !a.indices.is_empty() && !b.indices.is_empty() {",
        "if false && !a.indices.is_empty() && !b.indices.is_empty() {",
    ),
    (
        "clash-boundary-counts-as-inside",
        "crates/algorithms/reference/src/clash.rs",
        "Some(w > 0.75)",
        "Some(w > 0.25)",
    ),
    (
        "clash-coplanar-always-overlaps",
        "crates/algorithms/reference/src/clash.rs",
        "if coplanar_pair_overlaps(ta, tb) {",
        "if true {",
    ),
    (
        "clash-coplanar-never-overlaps",
        "crates/algorithms/reference/src/clash.rs",
        "if coplanar_pair_overlaps(ta, tb) {",
        "if false {",
    ),
    (
        "clash-touching-reported-as-clear",
        "crates/algorithms/reference/src/clash.rs",
        "if report.kind == Interference::Clear {\n                        report.kind = Interference::Touching;\n                    }\n                }\n                // `Coplanar` short-circuits",
        "if false {\n                        report.kind = Interference::Touching;\n                    }\n                }\n                // `Coplanar` short-circuits",
    ),
]

PROBES = PROBES + CLASH_PROBES

probe_lib.main(PROBES, tests=SUITES, command=(), timeout=1800, baseline=True)
