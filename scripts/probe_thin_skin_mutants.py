"""Mutation probe for skins thinner than the tolerance (#276).

A difference whose tool stops a rounding error short of its host's face
(or reaches that far past it) leaves a skin or sliver that thin. The mesh
compiler moves the operands onto each other's faces within the linear
tolerance before the mesh boolean (`crates/execution/compile/src/snap.rs`)
and reports the move; the exact boolean cuts an imprinted edge at its own
end where the other operand's boundary crosses it within the tolerance
(`imprint_cuts` in `axiolid-brep-boolean`'s `section.rs`). Each mutant
breaks the resolution, its report, or the tolerance limit, and must turn
the regression tests red.

Equivalent mutants, deliberately not listed:

- Not recording the imprint merge (`report::near` replaced by the bare
  comparison): every configuration that merges a cut there also reads the
  door's end face as coincident with the wall's face, by the same
  distance, so the report already bounds the move.
- The snap's own zero-tolerance guard (`eps <= 0.0`): with `eps = 0` no
  residue passes `> eps`, so nothing moves anyway; the floor mutant below
  is the one that would break zero tolerance.
- The cheap residue prefilter in `target` (`residue.abs() > eps`): the
  triangle distance it precedes is never smaller than the residue.
"""
import pathlib, subprocess, sys

ROOT = pathlib.Path(__file__).resolve().parents[1]
SNAP = "crates/execution/compile/src/snap.rs"
COMPILER = "crates/execution/compile/src/compiler.rs"
SECTION = "crates/algorithms/construction/brep-boolean/src/section.rs"
MESH = ["-p", "axiolid-mesh-compile", "--lib", "--test", "thin_skin"]
EXACT = ["-p", "axiolid-brep-boolean", "--test", "openings"]
BOTH = [MESH, EXACT, ["-p", "axiolid-mesh-compile", "--test", "thin_skin"]]

MUTANTS = [
    # Resolution.
    ('mesh boolean cuts the given operands', COMPILER,
     '            settled.subject.as_ref().unwrap_or(&subject.mesh),\n            settled.tool.as_ref().unwrap_or(tool),',
     '            &subject.mesh,\n            tool,', [MESH]),
    ('tool not moved onto the subject', SNAP,
     '    if let Some((moved_tool, by)) = onto(tool, subject, eps) {',
     '    if let Some((moved_tool, by)) = onto(tool, subject, eps).filter(|_| false) {', [MESH]),
    ('subject not moved onto the tool', SNAP,
     '    if let Some((moved_subject, by)) = onto(subject, fixed, eps) {',
     '    if let Some((moved_subject, by)) = onto(subject, fixed, eps).filter(|_| false) {', [MESH]),
    ('a vertex near an edge lands on one plane only', SNAP,
     '    let step = smallest_move(&planes)?;',
     '    let step = smallest_move(&planes[..1])?;', [MESH]),
    ('a vertex between parallel planes is moved anyway', SNAP,
     '                if same > PARALLEL * (1.0 + p.abs().max_element()) {\n                    return None;\n                }',
     '                let _ = same;', [MESH]),
    ('a move that flattens an operand is made', SNAP,
     '        if before != after && after.dot(before) <= 0.0 {',
     '        if before != after && after.dot(before) < 0.0 {', [MESH]),
    ('exact imprint keeps cuts next to the edge\'s ends', SECTION,
     '        let cuts = imprint_cuts(&curve, span, cuts, tolerance)?;',
     '        let _ = imprint_cuts;', [EXACT]),
    # Report.
    ('the snap is not reported', COMPILER,
     '        if settled.moved > 0.0 {\n            built.deviation.add(',
     '        if settled.moved > 0.0 && false {\n            built.deviation.add(', [MESH]),
    ('the snap reported smaller than the move', COMPILER,
     '                crate::deviation::DeviationBound::Certified(settled.moved),',
     '                crate::deviation::DeviationBound::Certified(settled.moved * 0.25),', [MESH]),
    ('an inner snap dropped by the outer boolean', COMPILER,
     '        built.deviation.carry_snaps(&subject.deviation);',
     '        let _ = &subject.deviation;', [MESH]),
    ('an inner snap dropped by the measured deviation', COMPILER,
     '                        deviation.carry_snaps(&built.deviation);',
     '                        let _ = &built.deviation;', [MESH]),
    # Tolerance limit.
    ('the snap reaches ten tolerances', SNAP,
     '    let eps = tolerance.linear();\n    let mut out = Settled {',
     '    let eps = 20.0 * tolerance.linear();\n    let mut out = Settled {', [MESH]),
    ('the snap has a floor at zero tolerance', SNAP,
     '    let eps = tolerance.linear();\n    let mut out = Settled {',
     '    let eps = tolerance.linear().max(1e-9);\n    let mut out = Settled {', [MESH]),
    ('a multi-plane move beyond the tolerance is made', SNAP,
     '    if step.length() > eps {',
     '    if step.length() > 1e6 * eps {', [MESH]),
    ('exact imprint merges ends a thousand tolerances off', SECTION,
     '        if !report::near(ToleranceDecisionKind::MergedPoints, gap, tolerance) {',
     '        if !report::near(ToleranceDecisionKind::MergedPoints, gap / 1e3, tolerance) {', BOTH),
]


def run(targets):
    for target in targets:
        code = subprocess.run(
            ["cargo", "test", "-q", *target],
            cwd=ROOT, capture_output=True, text=True, timeout=3600,
        ).returncode
        if code != 0:
            return code
    return 0


# Optional substrings select mutants by name.
if sys.argv[1:]:
    MUTANTS = [m for m in MUTANTS if any(a in m[0] for a in sys.argv[1:])]

survivors = []
for name, rel, old, new, targets in MUTANTS:
    path = ROOT / rel
    original = path.read_text()
    assert original.count(old) == 1, f"anchor for '{name}' not unique/found"
    path.write_text(original.replace(old, new))
    try:
        try:
            code = run(targets)
        except subprocess.TimeoutExpired:
            code = -1
    finally:
        path.write_text(original)
    status = "killed" if code != 0 else "SURVIVED"
    print(f"{status:8} {name}", flush=True)
    if code == 0:
        survivors.append(name)
print(f"{len(MUTANTS) - len(survivors)}/{len(MUTANTS)} killed")
sys.exit(1 if survivors else 0)
