"""Mutation probe for the exact extruder's profile check: each mutant must be caught.

The check of a union's regions (`refuse_invalid_region(&region.outer, ..)`)
has no mutant: a union of members returns pieces that are valid regions
whenever its members are, so no input reaches it; it stays as a guard and is
exercised through `refuse_invalid_region` by the unit tests.
"""
import pathlib

import probe_lib

ROOT = pathlib.Path(__file__).resolve().parents[1]
SRC = ROOT / "crates/algorithms/construction/construct/src"

MUTANTS = [
    ("profile check dropped", "extrude_exact.rs",
     "    // region builds a solid that is not a valid manifold (#401).\n"
     "    refuse_invalid_region(&lowered[0], &lowered[1..])?;\n",
     "    // region builds a solid that is not a valid manifold (#401).\n"),
    ("holes left out of the profile check", "extrude_exact.rs",
     "    refuse_invalid_region(&lowered[0], &lowered[1..])?;\n    let mut rings",
     "    refuse_invalid_region(&lowered[0], &[])?;\n    let mut rings"),
    ("composite members not checked", "extrude_exact.rs",
     "            refuse_invalid_region(&lowered[0], &lowered[1..])?;\n        }\n    }\n    let regions",
     "        }\n    }\n    let regions"),
    ("composite member holes not checked", "extrude_exact.rs",
     "            refuse_invalid_region(&lowered[0], &lowered[1..])?;\n        }\n    }\n    let regions",
     "            refuse_invalid_region(&lowered[0], &[])?;\n        }\n    }\n    let regions"),
    ("rings along one line not refused as flat", "extrude_exact.rs",
     "    for ring in lowered.iter().filter(|ring| is_flat(ring)) {",
     "    for ring in lowered.iter().filter(|_| false) {"),
    ("defect mapped to Ok", "extrude_exact.rs",
     "    if !definite.is_empty() {",
     "    if false {"),
    ("defect refused as unsupported", "extrude_exact.rs",
     "        return Err(GeomError::InvalidInput(format!(\n            \"exact extrusion profile is not a valid region: {definite:?}\"\n        )));",
     "        return Err(unsupported(\"exact extrusion profile is not a valid region\"));"),
    ("defect not named", "extrude_exact.rs",
     "is not a valid region: {definite:?}\"",
     "is not a valid region: {}\", definite.len()"),
    ("undecided nesting mapped to Ok", "extrude_exact.rs",
     "    if defects\n        .iter()\n        .any(|defect| !matches!(defect, Defect::Ring { .. }))\n    {",
     "    if false {"),
    ("undecided nesting refused as invalid", "extrude_exact.rs",
     "        return Err(unsupported(\n            \"exact extrusion profile with a contact undecided at the tolerance\",\n        ));",
     "        return Err(GeomError::InvalidInput(\"undecided\".into()));"),
    ("ring degeneracies refused as region defects", "extrude_exact.rs",
     "                    | Defect::DisconnectedInterior { .. }\n",
     "                    | Defect::DisconnectedInterior { .. }\n                    | Defect::Ring { .. }\n"),
    ("disconnected interior accepted", "extrude_exact.rs",
     "                    | Defect::DisconnectedInterior { .. }\n",
     ""),
]

probe_lib.main(MUTANTS, tests=[["-p", "axiolid-construct", "--test", "extrude_exact_validate"],
                               ["-p", "axiolid-construct", "--lib", "extrude_exact"]],
               root=SRC, timeout=900)
