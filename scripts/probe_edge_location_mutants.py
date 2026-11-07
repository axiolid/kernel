"""Mutation probe for locating shared edge samples on a pcurve (#413).

Each mutant weakens one decision of `crates/execution/compile/src/brep/
locate.rs`, `edges.rs` or the trim terms of `brep.rs`: the residual a located
sample may have and the bound it joins, the strictly monotone parameters,
the bracket of the search, the end limit of an exact boolean's reading and
the refinement of the sampled edge. `tests/shared_edge_location.rs` or the
modules' unit tests must fail on it. A build error also counts as killed.

Not listed, as not observable (no fixture, not equivalent):
- the strictly monotone check after the search, as a guard by construction:
  a bracket that starts at the sample before cannot return an earlier
  parameter, so the check only repeats the refusal at the bracket's start.
- the fast path that keeps a sample on its uniform parameter when the
  pcurve produced it: dropping it only lets the search find the same
  parameter to its tolerance.
- the iteration cap of the search: no fixture needs more than a few dozen.
"""
import probe_lib

L = "crates/execution/compile/src/brep/locate.rs"
B = "crates/execution/compile/src/brep.rs"
E = "crates/execution/compile/src/brep/edges.rs"
TARGETS = [
    ["-p", "axiolid-mesh-compile", "--test", "shared_edge_location"],
    ["-p", "axiolid-mesh-compile", "--lib", "brep::locate"],
    ["-p", "axiolid-mesh-compile", "--lib", "brep::tests"],
    ["-p", "axiolid-mesh-compile", "--lib", "brep::exact_tests"],
]

MUTANTS = [
    ("the residual of a sample is not limited", L,
     "        if !miss.is_finite() || miss > limit {\n",
     "        if !miss.is_finite() {\n"),
    ("the bracket's lower end is not refused: a minimum at the sample before is kept", L,
     "        if found - lo <= margin {\n            break;\n        }\n",
     ""),
    ("the bracket spans the whole domain", L,
     "    let mut width = 2.0 * step;\n",
     "    let mut width = 1.0;\n"),
    ("the bracket is never widened", L,
     "        width *= 2.0;\n",
     "        break;\n"),
    ("the ends take the interior limit", L,
     "        let limit = if i == 0 || i == n {\n            limits.end\n        } else {\n            limits.interior\n        };\n",
     "        let limit = limits.interior;\n"),
    ("the guess ignores the pace of the pcurve", L,
     "        uniform * to_walk / walked\n",
     "        uniform\n"),
    ("a located parameter keeps the search's rounding of a uniform one", L,
     "(located - t).abs() <= SNAP * (located - before).abs()",
     "false"),
    ("the corner of a hairpin is measured from the line, not the chord", B,
     "                ((p - p0).dot(along) / length).clamp(0.0, 1.0)\n",
     "                (p - p0).dot(along) / length\n"),
    ("the offset of the samples is added to the trim bound again", B,
     "            bound.trim(&pcurve, &ts, &params);\n",
     "            bound.trim(&pcurve, &ts, &params);\n"
     "            for (&index, p) in shared.iter().zip(&params) {\n"
     "                let exact = surface.evaluate(p.x, p.y)?;\n"
     "                bound.trim += (state.mesh.positions[index as usize] - exact).length();\n"
     "            }\n"),
    ("the trim is bounded over the uniform parameters, not the located ones", B,
     "            bound.trim(&pcurve, &ts, &params);\n",
     "            bound.trim(&pcurve, &uniform_times(&pcurve, shared.len() - 1), &params);\n"),
    ("a span across a corner is bounded without the corner's distance", B,
     "            off = off.max((p - (p0 + along * t)).length());\n",
     "            off = off.max(0.0 * (p - (p0 + along * t)).length());\n"),
    ("the pieces of a span across a corner are not bounded", B,
     "        worst_piece = worst_piece.max(axiolid_reference::bound::chord_bound2(curve, from, cut)?);\n",
     "        let _ = axiolid_reference::bound::chord_bound2(curve, from, cut)?;\n"),
    ("the sampled edge is not refined for the other uses", E,
     "        if n >= MAX_CURVED_EDGE_SEGMENTS || located_within_budget(ctx, plan, source, &points)? {\n",
     "        if true {\n"),
    ("a located span over the budget does not refine the edge", E,
     "            if fits == Some(false) {\n",
     "            if fits == Some(true) && false {\n"),
    ("a vertex is welded only within the chord budget", B,
     "        if !disagreement.is_finite() || disagreement > state.endpoint_limit {\n",
     "        if !disagreement.is_finite() || disagreement > state.chord_error {\n"),
    ("the located ends are limited by the chord budget", B,
     "                        end: state.endpoint_limit,\n",
     "                        end: ctx.chord_error,\n"),
]

probe_lib.main(MUTANTS, tests=TARGETS, command=probe_lib.CARGO_TEST + ("--all-features",), timeout=1800)
