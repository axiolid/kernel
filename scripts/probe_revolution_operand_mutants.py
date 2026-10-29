"""Mutation probe for revolutions, and any compiled solid, as operands of an
exact boolean and of its measured deviation (#378).

Each mutant reintroduces a way a revolution operand, or the compiled-solid
reading of operands, failed or could fail -- the revolution refused again
as subject or as tool (bare or under a placement), a boolean refused again
as tool, a symmetric difference let through, the search settling at 10%
above its lower bound across the budget -- and must turn a test red.
`exact_placed_boolean.rs` checks the volumes of a torus cut through its axis
(as subject and as tool) and of a ring groove turned from a rectangle;
`boolean_deviation.rs` measures the box and bent tube of the issue under
every operator against the exact result. `--only NAME` runs only the
mutants whose names contain NAME.

Not listed: the `touching or crossing the axis` refusal of a full-turn
solid cylinder and the tapered revolution's: they are `compile`'s own and
are probed by `probe_exact_mutants.py`.
"""
import probe_lib

C = "crates/execution/compile/src/exact/boolean.rs"
S = "crates/execution/compile/src/certify.rs"
PLACED = ["-p", "axiolid-mesh-compile", "--test", "exact_placed_boolean"]
MEASURED = ["-p", "axiolid-mesh-compile", "--test", "boolean_deviation"]

SUBJECT = "        let subject = self.compile(left)?;\n"
TOOL = "        let tool = self.compile(right)?;\n"


def refuse(operand, kinds):
    """Refuse the operand, read through its placements, when it is `kinds`."""
    return (
        f"        let mut at = {operand};\n"
        "        while let Some(GeometryNode::Instance(i)) = self.graph.get(at) {\n"
        "            at = i.source;\n"
        "        }\n"
        f"        if matches!(self.graph.get(at), Some(GeometryNode::SolidOperation({kinds}))) {{\n"
        '            return Err(unsupported("exact boolean operand that is not an extrusion"));\n'
        "        }\n"
    )


REVOLUTION = "SolidOperation::Revolution { .. }"
BOOLEAN = "SolidOperation::Boolean { .. }"

MUTANTS = [
    ("a revolution tool refused again", C,
     TOOL, refuse("right", REVOLUTION) + TOOL, [PLACED, MEASURED]),
    ("a revolution subject refused again", C,
     SUBJECT, refuse("left", REVOLUTION) + SUBJECT, [PLACED, MEASURED]),
    ("a revolution subject or tool refused again", C,
     SUBJECT, refuse("left", REVOLUTION) + refuse("right", REVOLUTION) + SUBJECT,
     [PLACED, MEASURED]),
    ("a boolean tool refused again", C,
     TOOL, refuse("right", BOOLEAN) + TOOL, [PLACED, MEASURED]),
    ("the search settles across the budget", S,
     "            || (top.upper <= TIGHT * lower && !straddles)\n",
     "            || (top.upper <= TIGHT * lower)\n", [MEASURED]),
]

probe_lib.main(MUTANTS, timeout=1800, reject_uncompiled=True)
