"""Mutation probe for exact unions and intersections of placed operands and
their measured deviation (#376).

Each mutant reintroduces a way the operator handling of
`compile_boolean` failed or could fail -- a union or an intersection
refused again, the operator the exact boolean runs not the node's, the
tool restriction (#326) dropped with the operator gate, the emptied result
named for the wrong operator -- and must turn a test red.
`exact_placed_boolean.rs` checks the volumes and refusals of the exact
unions and intersections; `boolean_deviation.rs` measures their meshes
against the exact result. `--only NAME` runs only the mutants whose names
contain NAME.

Equivalent mutants, deliberately not listed: the bare-prism path before
the gate (it takes every operator on its own and is probed by
`probe_boolean_deviation_mutants.py` and the prism tests).
"""
import probe_lib

C = "crates/execution/compile/src/exact/boolean.rs"
PLACED = ["-p", "axiolid-mesh-compile", "--test", "exact_placed_boolean"]
MEASURED = ["-p", "axiolid-mesh-compile", "--test", "boolean_deviation"]

REFUSE = (
    "        let tool = self.placed(right)?;\n",
    "        let tool = self.placed(right)?;\n"
    "        if operator == BooleanOperator::{op} {{\n"
    "            return Err(unsupported(\n"
    '                "exact union or intersection of placed operands",\n'
    "            ));\n"
    "        }}\n",
)
RUN = "boolean_with_report(&subject, &tool, operator, self.options.tolerance())"

MUTANTS = [
    ("union of placed operands refused again", C,
     REFUSE[0], REFUSE[1].format(op="Union"), [PLACED, MEASURED]),
    ("intersection of placed operands refused again", C,
     REFUSE[0], REFUSE[1].format(op="Intersection"), [PLACED, MEASURED]),
    ("every placed boolean run as a difference", C,
     RUN, "boolean_with_report(&subject, &tool, BooleanOperator::Difference, self.options.tolerance())",
     [PLACED, MEASURED]),
    ("a union run as an intersection", C,
     RUN, "boolean_with_report(&subject, &tool, match operator { BooleanOperator::Union => BooleanOperator::Intersection, other => other }, self.options.tolerance())",
     [PLACED, MEASURED]),
    ("an intersection run as a union", C,
     RUN, "boolean_with_report(&subject, &tool, match operator { BooleanOperator::Intersection => BooleanOperator::Union, other => other }, self.options.tolerance())",
     [PLACED, MEASURED]),
    ("tool restriction dropped with the operator gate", C,
     "        if !matches!(tool, Placed::Extrusion | Placed::Mesh) {",
     "        if false && !matches!(tool, Placed::Extrusion | Placed::Mesh) {",
     [PLACED, MEASURED]),
    ("an emptied result named for a difference whatever the operator", C,
     "BooleanError::EmptyResult => emptied(operator),",
     "BooleanError::EmptyResult => emptied(BooleanOperator::Difference),",
     [PLACED]),
]

probe_lib.main(MUTANTS, timeout=1800, reject_uncompiled=True)
