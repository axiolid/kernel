#!/usr/bin/env python3
"""Mutation probes for the oracle and conformance suite (ADR 0017 sections 5-6).

The oracle is the thing everything else is judged against, so it needs the
strictest proof that its own tests can fail. Each probe injects one semantic
defect and asserts the suites FAIL. A surviving probe means the gate is blind.

Sources are restored byte-identically afterwards.
"""
import os
from pathlib import Path

import probe_lib

ROOT = Path(
    os.environ.get("PROBE_ROOT", probe_lib.ROOT)
)

SUITES = [
    ["cargo", "test", "-p", "axiolid-reference", "--all-features", "--test", "oracle"],
    ["cargo", "test", "-p", "axiolid-mesh-boolean-boolmesh", "--all-features", "--test", "conformance"],
]

PROBES = [
    # --- oracle correctness (ADR 0017 section 5) ---
    (
        "oracle-union-of-nested-returns-inner",
        "crates/algorithms/reference/src/boolean.rs",
        "(BooleanOperator::Union, Arrangement::ToolInsideSubject) => subject.clone(),",
        "(BooleanOperator::Union, Arrangement::ToolInsideSubject) => tool.clone(),",
    ),
    (
        "oracle-difference-forgets-the-cavity",
        "crates/algorithms/reference/src/boolean.rs",
        """            (
                BooleanOperator::Difference | BooleanOperator::SymmetricDifference,
                Arrangement::ToolInsideSubject,
            ) => concatenate(subject, &reversed(tool)),""",
        """            (
                BooleanOperator::Difference | BooleanOperator::SymmetricDifference,
                Arrangement::ToolInsideSubject,
            ) => subject.clone(),""",
    ),
    (
        "oracle-cavity-is-not-reversed",
        "crates/algorithms/reference/src/boolean.rs",
        "    for triangle in indices.as_chunks_mut::<3>().0 {\n        triangle.swap(0, 1);\n    }",
        "    for _triangle in indices.as_chunks_mut::<3>().0 {}",
    ),
    (
        "oracle-self-difference-is-not-empty",
        "crates/algorithms/reference/src/boolean.rs",
        """            (
                BooleanOperator::Difference | BooleanOperator::SymmetricDifference,
                Arrangement::Identical,
            ) => empty(),""",
        """            (
                BooleanOperator::Difference | BooleanOperator::SymmetricDifference,
                Arrangement::Identical,
            ) => subject.clone(),""",
    ),
    (
        "oracle-disjoint-intersection-is-not-empty",
        "crates/algorithms/reference/src/boolean.rs",
        "(BooleanOperator::Intersection, Arrangement::Disjoint) => empty(),",
        "(BooleanOperator::Intersection, Arrangement::Disjoint) => subject.clone(),",
    ),
    # --- oracle honesty: it must refuse what it cannot answer ---
    (
        "oracle-guesses-at-interpenetration",
        "crates/algorithms/reference/src/boolean.rs",
        "    if surfaces_intersect(subject, tool, options)? {",
        "    if false && surfaces_intersect(subject, tool, options)? {",
    ),
    (
        "oracle-crossing-test-reopens-the-shared-edge-hole",
        "crates/algorithms/reference/src/boolean.rs",
        "        if !(positive && negative) {\n            return true;\n        }",
        "        if !(positive && negative)\n            && !signs.contains(&Sign::Zero)\n        {\n            return true;\n        }",
    ),
    (
        "oracle-accepts-empty-operands",
        "crates/algorithms/reference/src/boolean.rs",
        "        if mesh.positions.is_empty() || mesh.indices.is_empty() {",
        "        if false {",
    ),
    # --- containment classification ---
    (
        "oracle-point-on-surface-counts-as-inside",
        "crates/algorithms/reference/src/boolean.rs",
        "        if side_origin == Sign::Zero {\n            // The point is ON the surface: neither inside nor outside.\n            return Some(false);\n        }",
        "        if side_origin == Sign::Zero {\n            return Some(true);\n        }",
    ),
    (
        "oracle-parity-inverted",
        "crates/algorithms/reference/src/boolean.rs",
        "    Some(crossings % 2 == 1)",
        "    Some(crossings % 2 == 0)",
    ),
    # --- conformance suite must detect defects (ADR 0017 section 6) ---
    (
        "suite-treats-failure-as-conformant",
        "crates/contracts/operations/mesh-boolean/src/conformance.rs",
        "        !self\n            .checks\n            .iter()\n            .any(|check| matches!(check.outcome, Outcome::Failed { .. }))",
        "        true",
    ),
    (
        "suite-counts-skips-as-passes",
        "crates/contracts/operations/mesh-boolean/src/conformance.rs",
        "            .filter(|check| check.outcome == Outcome::Passed)\n            .count()",
        "            .filter(|check| check.outcome != Outcome::Passed || true)\n            .count()",
    ),
    (
        "suite-skips-the-disjoint-algebra-check",
        "crates/contracts/operations/mesh-boolean/src/conformance.rs",
        "    check_disjoint_algebra(provider, &mut report);",
        "    let _ = check_disjoint_algebra;",
    ),
    (
        "suite-ignores-a-wrong-union-volume",
        "crates/contracts/operations/mesh-boolean/src/conformance.rs",
        "        if (measured - 2.0).abs() < 1e-9 {",
        "        if true {",
    ),
    (
        "suite-accepts-an-unrefused-empty-operand",
        "crates/contracts/operations/mesh-boolean/src/conformance.rs",
        "            Ok(_) => Outcome::Failed {\n                detail: \"an empty mesh is not a solid and must be refused\".into(),\n            },",
        "            Ok(_) => Outcome::Passed,",
    ),
    # --- registration gate ---
    (
        "registration-admits-non-conformant-providers",
        "crates/execution/dispatch/src/boolean.rs",
        "        if !report.is_conformant() {\n            return Err(Box::new(report));\n        }",
        "        if false {\n            return Err(Box::new(report));\n        }",
    ),
]

probe_lib.main(PROBES, tests=SUITES, root=ROOT, cwd=ROOT, command=(), timeout=1800, baseline=True)
