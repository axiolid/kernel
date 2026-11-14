"""Mutation probe for the deviation oracle and tessellator at poles (#379).

A face reaching a sphere pole or a cone apex has a loop that leaves the pole
along its line, a gap in parameters that is no gap on the surface. Each
mutant weakens one decision of `crates/execution/compile/src/deviation/
boolean.rs` (the pole-gap rule, the cone's box ending at its apex) or of
`brep/pole.rs` (the named refusal of a trim closing across a collapsed edge
of a surface without a pole kind, P1 of #463) and must turn a test red: the
partial-cone case of `tests/boolean_deviation.rs`, the `open` unit test of
the oracle, or `a_collapsed_edge_of_a_spline_surface_is_refused` in
`tests/pole_faces.rs`. A build error counts as killed.

Faces that do reach a sphere's or a cone's pole are tessellated since #463;
`scripts/probe_pole_faces_mutants.py` probes that.

Not listed, as equivalent at test resolution: the pole gap skipped at any
height, not only where the loop leaves the pole at the height it reached it.
A loop jumping between a sphere's two poles would close wrongly, and no
fixture builds one.
"""
import probe_lib

D = "crates/execution/compile/src/deviation/boolean.rs"
P = "crates/execution/compile/src/brep/pole.rs"
TESTS = [
    ["-p", "axiolid-mesh-compile", "--test", "boolean_deviation"],
    ["-p", "axiolid-mesh-compile", "--lib"],
]
POLE_FACES = [["-p", "axiolid-mesh-compile", "--test", "pole_faces"]]

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
    ("sphere box not ended at its poles", D,
     "        if matches!(surface, Surface::Sphere(_)) {",
     "        if false {"),
    ("trim across a collapsed edge tessellated", P,
     "    end != start && pole(end) && pole(start)",
     "    end != start && pole(end) && pole(start) && false", POLE_FACES),
]

probe_lib.main(MUTANTS, tests=TESTS, command=probe_lib.CARGO_TEST, timeout=1800)
