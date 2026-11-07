"""Mutation probe for the exit check of a boolean's result (#445).

`axiolid-brep-boolean` never lets a window's crossing replace a feature cut
(`section.rs`), and refuses a result whose vertex lies off the curve of an
edge it bounds (`consistency.rs`, `BooleanError::VertexOffCurve`). Each
mutant lets a window end replace a cut, loosens or removes one bound of the
check (the cap by the caller's tolerance, the bound of a non-exact report,
the NaN refusal), measures nothing on the edges or the pcurves, drops the
check from `run` or from the recentred chain, or reads the retried tolerance
as the caller's, and must turn a test red: the far corbel sweep, the
sub-ulp refusal, the consistency unit tests or the report's.

Each mutant runs in its own process group under a timeout that kills the
whole group; the mutated file is always restored. A missing or ambiguous
anchor aborts the probe.
"""
import probe_lib

D = "crates/algorithms/construction/brep-boolean/src/"
S = D + "section.rs"
C = D + "consistency.rs"
L = D + "lib.rs"
R = D + "report.rs"
BOOLEAN = [
    ["-p", "axiolid-brep-boolean", "--test", "far_origin"],
    ["-p", "axiolid-brep-boolean", "--lib", "consistency::"],
    ["-p", "axiolid-brep-boolean", "--lib", "section::tests::a_window_end"],
    ["-p", "axiolid-brep-boolean", "--test", "tangent_contact"],
    ["-p", "axiolid-brep-boolean", "--lib", "report::"],
]

MUTANTS = [
    ('a window end replaces a cut', S,
     '        .filter(|&c| !cuts.iter().any(|&f| same_cut(c, f)))\n',
     '        .filter(|&c| cuts.is_empty() || c.is_finite())\n', BOOLEAN),
    ('an exact report bounded by the floor alone', C,
     '        (true, true) => floor.min(eps),',
     '        (true, true) => floor,', BOOLEAN),
    ('a report with readings bounded above eps', C,
     '        (false, _) => eps,',
     '        (false, _) => eps.max(floor),', BOOLEAN),
    ('a residue that is not a number passes', C,
     '    residue.is_nan() || residue > bound',
     '    residue > bound', BOOLEAN),
    ('an edge end is not measured', C,
     '        worst = worst.max(distance(point, vertex)?);\n    }\n    Ok(worst)\n}\n\n/// How far `start` and `end` lie from `pcurve`',
     '        worst = worst.max(distance(point, vertex)? * 0.0);\n    }\n    Ok(worst)\n}\n\n/// How far `start` and `end` lie from `pcurve`', BOOLEAN),
    ('a use end is not measured', C,
     '        let point = surface::evaluate(surface, uv.x, uv.y).map_err(|_| BooleanError::Evaluation)?;\n        worst = worst.max(distance(point, vertex)?);',
     '        let point = surface::evaluate(surface, uv.x, uv.y).map_err(|_| BooleanError::Evaluation)?;\n        worst = worst.max(distance(point, vertex)? * 0.0);', BOOLEAN),
    ('the check removed from the cut', L,
     '    consistency::check(&result, tolerance)?;\n',
     '', BOOLEAN),
    ('the check removed after the move back', L,
     '            .and_then(|result| consistency::check(&result, tolerance).map(|()| result))\n',
     '', BOOLEAN),
    ('the caller tolerance read from the retried one', R,
     '.map(|s| s.caller))',
     '.map(|s| s.tolerance))', BOOLEAN),
]

probe_lib.main(MUTANTS, timeout=1200)
