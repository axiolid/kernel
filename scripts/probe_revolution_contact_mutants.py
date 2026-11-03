"""Mutation probe for a revolution's contact with a face (#381, #382).

#381: a box less a torus whose tube leaves through the box's top face. The
face touches the latitude circles where the tube's quarter faces meet, the
traced section's double root splits in two roots 2e-8 apart, and they must be
one cut. #382: a column united with a revolved bend whose start cap lies in
the column's foot. The plane through the torus axis meets the torus in the
cap's circles, which the cap imprints in closed form; a traced piece along an
edge must take the edge's curve, and the plane through the axis meets the
torus in its two meridian circles in closed form. Each mutant breaks one
step and must turn a test red.

Equivalent mutant, deliberately not listed: the `SectionAlongEdge` residue
gate in `along_edge_curve` set to `true`. A piece along an edge was cut at
that edge's ends within `IncidentPoint`, so its ends never lie further off
the edge than the tolerance in any constructible case; the gate is a
fail-closed guard, equivalent under every test.
"""
import probe_lib

S = "crates/algorithms/construction/brep-boolean/src/section.rs"
T = "crates/algorithms/parametric/nurbs/src/torus_section.rs"
TOP = ["-p", "axiolid-brep-boolean", "--test", "curved", "a_box_less_a_torus"]
BEND = ["-p", "axiolid-brep-boolean", "--test", "curved", "bend"]
MERIDIAN = ["-p", "axiolid-nurbs", "--test", "exact_intersection", "meridian_circles"]

MUTANTS = [
    ("a traced section's cuts kept apart", S,
     "matches!(curve, Curve3::Line(_) | Curve3::ImplicitSection(_))",
     "matches!(curve, Curve3::Line(_))",
     TOP),
    ("a traced piece along an edge keeps its trace", S,
     "    if !matches!(curve, Curve3::ImplicitSection(_) | Curve3::PairSection(_)) {\n"
     "        return Ok((curve, span));",
     "    if true {\n        return Ok((curve, span));",
     BEND),
    ("the adopted stretch runs the other way round", S,
     "            if tm_up < te_up {",
     "            if tm_up > te_up {",
     BEND),
    ("no meridian closed form", T,
     "        if is_zero(&c) {\n            let n = p.frame.z;",
     "        if false && is_zero(&c) {\n            let n = p.frame.z;",
     MERIDIAN),
    ("meridian circles on one side", T,
     "vec![meridian(along), meridian(-along)]",
     "vec![meridian(along), meridian(along)]",
     MERIDIAN),
]

probe_lib.main(MUTANTS, tests=TOP, timeout=1800, reject_uncompiled=True)
