"""Mutation probe for open surfaces compiled through compile_mesh (#311).

Each mutant refuses loose faces or reports a trimmed sheet as a solid,
under-reports or overruns a trimmed patch's bound, certifies a cell across
a crease, turns a trimmed sheet the wrong way, misreads or fails to refuse
its parameter range, and must turn a test red.
"""
import probe_lib

T = "crates/execution/compile/src/trimmed.rs"
BR = "crates/execution/compile/src/brep.rs"
CO = "crates/execution/compile/src/compiler.rs"
TESTS = [
    ["-p", "axiolid-mesh-compile", "--test", "open_surfaces"],
    ["-p", "axiolid-mesh-compile", "--test", "surface_models"],
]

MUTANTS = [
    ('loose faces refused', BR, '                (vec![(&loose, false)], MeshClosure::Surface)', '                return Err(GeomError::InvalidInput("loose faces".into()));', TESTS),
    ('trimmed sheet reported a solid', CO, '                Built::with_closure(mesh, axiolid_mesh_compile_contract::MeshClosure::Surface)\n                    .with_deviation(crate::deviation::Deviation::one(\n                        crate::deviation::DeviationPath::TrimmedSurface,', '                Built::with_closure(mesh, axiolid_mesh_compile_contract::MeshClosure::Solid)\n                    .with_deviation(crate::deviation::Deviation::one(\n                        crate::deviation::DeviationPath::TrimmedSurface,', TESTS),
    ('loose faces reported a solid', BR, '                (vec![(&loose, false)], MeshClosure::Surface)', '                (vec![(&loose, false)], MeshClosure::Solid)', TESTS),
    ('range longer than a period accepted', T, '    if period.is_some_and(|period| hi - lo > period) {', '    if false && period.is_some_and(|period| hi - lo > period) {', TESTS),
    ('trimmed bound under-reported', T, '                    let bound = 0.125 * (', '                    let bound = 0.03125 * (', TESTS),
    ('too few cells for the budget', T, '    let count = (0.5 * length * (k / chord).sqrt()).ceil();', '    let count = (0.25 * length * (k / chord).sqrt()).ceil();', TESTS),
    ('sense ignored for the facing', T, '    let flip = trim.u_sense != trim.v_sense;', '    let flip = false;', TESTS),
    ('seam not wrapped forward', T, '        (true, Some(period)) => (first, second + period),', '        (true, Some(period)) => (second, first),', TESTS),
    ('sense against a non-periodic parameter accepted', T, '        (false, _) if second < first => (second, first),', '        (false, _) => (first.min(second), first.max(second)),', TESTS),
    ('cell across a crease certified', T, '                Some(cell) if cell.smooth => {', '                Some(cell) => {', TESTS),
    ('pole triangles kept degenerate', T, '                if pp == qp || qp == rp || rp == pp {', '                if false {', TESTS),
]

probe_lib.main(MUTANTS, timeout=1200)
