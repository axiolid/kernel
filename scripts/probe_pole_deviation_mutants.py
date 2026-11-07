"""Mutation probe for the deviation oracle and tessellator at poles (#379).

A face reaching a sphere pole or a cone apex has a loop that leaves the pole
along its line, a gap in parameters that is no gap on the surface. Each
mutant weakens one decision of `crates/execution/compile/src/deviation/
boolean.rs` (the pole-gap rule, the cone's box ending at its apex) or of
`brep.rs` (the named refusal of a trim closing across a pole) and must turn a
test red: the partial-cone case of `tests/boolean_deviation.rs`, the `open`
unit test of the oracle, or `a_face_closing_across_an_apex_is_refused_by_name`.
A build error counts as killed.

Not listed, as equivalent at test resolution: the pole gap skipped at any
height, not only where the loop leaves the pole at the height it reached it.
A loop jumping between a sphere's two poles would close wrongly, and no
fixture builds one.
"""
import probe_lib

D = "crates/execution/compile/src/deviation/boolean.rs"
B = "crates/execution/compile/src/brep.rs"
TESTS = [
    ["-p", "axiolid-mesh-compile", "--test", "boolean_deviation"],
    ["-p", "axiolid-mesh-compile", "--lib"],
]

MUTANTS = [
    ("pole gap not skipped", D,
     "                if end.y == start.y && at_pole(surface, end) && at_pole(surface, start) {",
     "                if false {"),
    ("every gap skipped", D,
     "                if end.y == start.y && at_pole(surface, end) && at_pole(surface, start) {",
     "                if true {"),
    ("cone box not ended at its apex", D,
     "            if slope != 0.0 {\n                let apex",
     "            if false {\n                let apex"),
    ("cone box ended on the wrong side of its apex", D,
     "                if slope < 0.0 {\n                    domain.hi.y",
     "                if slope > 0.0 {\n                    domain.hi.y"),
    ("trim across a pole tessellated", B,
     "    end != start && pole(end) && pole(start)",
     "    end != start && pole(end) && pole(start) && false"),
]

probe_lib.main(MUTANTS, tests=TESTS, command=probe_lib.CARGO_TEST, timeout=1800)
