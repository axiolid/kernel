#!/usr/bin/env python3
"""Mutation probes for the boolean contract (ADR 0017 sections 1-4).

Each probe injects one semantic defect and asserts the contract suites FAIL.
A surviving probe means the gate is blind to that defect and the tests are
decorative. Sources are restored byte-identically afterwards.
"""
from pathlib import Path

import probe_lib

ROOT = Path(__file__).resolve().parents[1]

SUITES = [
    ["cargo", "test", "-p", "axiolid-dispatch", "--all-features", "--test", "mesh_boolean_contract"],
    ["cargo", "test", "-p", "axiolid-mesh-boolean-boolmesh", "--all-features", "--test", "symmetric_difference"],
    # `csg_deferral` was deleted when the deferral ended (ADR 0017 fully
    # landed); the conformance suite superseded it.
    ["cargo", "test", "-p", "axiolid-mesh-boolean-boolmesh", "--all-features", "--test", "conformance"],
]

# (name, relative path, find, replace)
PROBES = [
    (
        "ops-drop-symmetric-difference",
        "crates/foundation/core/src/operation.rs",
        "    pub const ALL: [Self; 4] = [\n        Self::Union,\n        Self::Intersection,\n        Self::Difference,\n        Self::SymmetricDifference,\n    ];",
        "    pub const ALL: [Self; 3] = [Self::Union, Self::Intersection, Self::Difference];",
    ),
    (
        "ops-difference-claims-commutative",
        "crates/foundation/core/src/operation.rs",
        "        !matches!(self, Self::Difference)",
        "        true",
    ),
    (
        "registry-skips-precondition-validation",
        "crates/execution/dispatch/src/boolean.rs",
        "        SolidRequirements::Oriented.validate_operands(subject, &[tool])?;",
        "        let _ = SolidRequirements::Oriented;",
    ),
    (
        "solid-accepts-inside-out",
        "crates/contracts/common/mesh/src/solid.rs",
        "        if six_volume < 0.0 {",
        "        if false {",
    ),
    (
        "solid-accepts-zero-volume",
        "crates/contracts/common/mesh/src/solid.rs",
        "        if six_volume == 0.0 {",
        "        if false {",
    ),
    (
        "cancellation-never-fires",
        "crates/contracts/common/base/src/cancel.rs",
        "        if self.is_cancelled() {\n            return Err(GeomError::Cancelled);\n        }",
        "        if false {\n            return Err(GeomError::Cancelled);\n        }",
    ),
    (
        "cancellation-token-clones-dont-share",
        "crates/contracts/common/base/src/cancel.rs",
        "        Arc::ptr_eq(&self.flag, &other.flag)",
        "        self.is_cancelled() == other.is_cancelled()",
    ),
    (
        "evidence-absorb-loses-sub-operations",
        "crates/contracts/operations/mesh-boolean/src/evidence.rs",
        "        self.sub_operations += other.sub_operations;",
        "        self.sub_operations = other.sub_operations;",
    ),
    (
        "evidence-absorb-overwrites-input-counts",
        "crates/contracts/operations/mesh-boolean/src/evidence.rs",
        "        self.output_triangles = other.output_triangles;",
        "        self.subject_triangles = other.subject_triangles;\n        self.output_triangles = other.output_triangles;",
    ),
    (
        "xor-reports-one-pass",
        "crates/contracts/operations/mesh-boolean/src/contract.rs",
        "    evidence.sub_operations = 3;",
        "    evidence.sub_operations = 1;",
    ),
    (
        "xor-skips-the-final-difference",
        "crates/contracts/operations/mesh-boolean/src/contract.rs",
        "    if intersection.mesh.indices.is_empty() {",
        "    if true {",
    ),
]

probe_lib.main(PROBES, tests=SUITES, command=(), timeout=1800, baseline=True)
