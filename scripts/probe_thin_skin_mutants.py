"""Mutation probe for skins at the rounding scale (#276, #291).

A difference whose tool stops a rounding error short of its host's face
(or reaches that far past it) leaves a skin or sliver that thin. The mesh
compiler sets the operands' coordinates onto each other's axis-aligned
faces before the mesh boolean, by at most a rounding residue of their
coordinates, exactly, and reports the move; where the snapped boolean is
refused it cuts the operands as given
(`crates/execution/compile/src/snap.rs`, `build_boolean` in `compiler.rs`).
The exact boolean cuts an imprinted edge at its own end where the other
operand's boundary crosses it within the tolerance (`imprint_cuts` in
`axiolid-brep-boolean`'s `section.rs`), and retries a boolean that merge
left refused with ends merged only within the rounding floor (`lib.rs`,
`report.rs`). Each mutant breaks the resolution, its reach, its exact
landing, the never-worse guard, or the report, and must turn the
regression tests red.

Equivalent mutants, deliberately not listed:

- Not recording the imprint merge (`report::near` replaced by the bare
  comparison): every configuration that merges a cut there also reads the
  door's end face as coincident with the wall's face, by the same
  distance, so the report already bounds the move.
- The snap's own zero-tolerance guard (`eps <= 0.0`): with `eps = 0` the
  reach is `min(0, ...) = 0` anyway, and nothing moves; the floor mutant
  below is the one that would break zero tolerance.
- The cheap residue prefilter in `onto` (`(face.at - p[k]).abs() > reach`):
  the triangle distance it precedes is never smaller than the residue.
- Which refusal is returned when the snapped and the given operands are
  both refused: the tests' refusing boolean refuses both alike.
- The exact end merge reaching a thousand tolerances (`gap / 1e3`):
  killed before #291 by the door ten tolerances past the face, which it
  left unsewn; now masked by the retry, which cuts that door with ends
  merged only within the rounding floor, as before #276. The retry itself
  is pinned by the three `exact ...` mutants below.
- The pinch check on the snapped result alone: `cut` checks both results
  the same way, so a pinch is a refusal on either path (the pinch itself
  is pinned by `tests/tangent_void.rs`).
"""
import pathlib, subprocess, sys

ROOT = pathlib.Path(__file__).resolve().parents[1]
SNAP = "crates/execution/compile/src/snap.rs"
COMPILER = "crates/execution/compile/src/compiler.rs"
SECTION = "crates/algorithms/construction/brep-boolean/src/section.rs"
BB_LIB = "crates/algorithms/construction/brep-boolean/src/lib.rs"
BB_REPORT = "crates/algorithms/construction/brep-boolean/src/report.rs"
MESH = ["-p", "axiolid-mesh-compile", "--lib", "--test", "thin_skin"]
FLUSH = ["-p", "axiolid-mesh-compile", "--test", "flush_through_cut"]
EXACT = ["-p", "axiolid-brep-boolean", "--test", "openings"]

MUTANTS = [
    # Resolution.
    ('mesh boolean cuts the given operands', COMPILER,
     '                settled.subject.as_ref().unwrap_or(&subject.mesh),\n                settled.tool.as_ref().unwrap_or(tool),',
     '                &subject.mesh,\n                tool,', [MESH]),
    ('tool not moved onto the subject', SNAP,
     '    if let Some((moved_tool, by)) = onto(tool, subject, reach) {',
     '    if let Some((moved_tool, by)) = onto(tool, subject, reach).filter(|_| false) {', [MESH]),
    ('subject not moved onto the tool', SNAP,
     '    if let Some((moved_subject, by)) = onto(subject, fixed, reach) {',
     '    if let Some((moved_subject, by)) = onto(subject, fixed, reach).filter(|_| false) {', [MESH]),
    ('a vertex near an edge lands on one plane only', SNAP,
     '                to[k] = *at;',
     '                if to == *p {\n                    to[k] = *at;\n                }', [MESH]),
    ('a value between parallel planes is moved anyway', SNAP,
     '                    if *at != Some(face.at) {\n                        *at = None;\n                    }',
     '                    *at = Some(face.at);', [MESH]),
    ('a move that flattens an operand is made', SNAP,
     '        if before != after && after.dot(before) <= 0.0 {',
     '        if before != after && after.dot(before) < 0.0 {', [MESH]),
    ('exact imprint keeps cuts next to the edge\'s ends', SECTION,
     '        let cuts = imprint_cuts(&curve, span, cuts, tolerance)?;',
     '        let _ = imprint_cuts;', [EXACT]),
    # Reach: the rounding scale, not the tolerance.
    ('the reach is the tolerance', SNAP,
     '    (REACH_EPSILONS * Scalar::EPSILON * scale).min(eps / 3.0_f64.sqrt())',
     '    (eps / 3.0_f64.sqrt())', [MESH, FLUSH]),
    ('the reach is a thousand times wider', SNAP,
     'pub(crate) const REACH_EPSILONS: Scalar = 16.0;',
     'pub(crate) const REACH_EPSILONS: Scalar = 16000.0;', [MESH, FLUSH]),
    ('the reach is a quarter as wide', SNAP,
     'pub(crate) const REACH_EPSILONS: Scalar = 16.0;',
     'pub(crate) const REACH_EPSILONS: Scalar = 4.0;', [MESH, FLUSH]),
    ('the reach ignores the coordinates\' magnitude', SNAP,
     '        .fold(0.0, Scalar::max);\n    (REACH_EPSILONS',
     '        .fold(1.0, Scalar::min);\n    (REACH_EPSILONS', [MESH, FLUSH]),
    ('the reach is not capped by the tolerance', SNAP,
     '    (REACH_EPSILONS * Scalar::EPSILON * scale).min(eps / 3.0_f64.sqrt())',
     '    REACH_EPSILONS * Scalar::EPSILON * scale', [MESH]),
    ('the snap has a floor at zero tolerance', SNAP,
     '    let eps = tolerance.linear();\n    if eps <= 0.0',
     '    let eps = tolerance.linear().max(1e-9);\n    if eps <= 0.0', [MESH]),
    # Landing: exactly, on exactly axis-aligned faces, a face whole.
    ('a face axis-aligned only within rounding is landed on', SNAP,
     '    (0..3).find(|&k| corners[0][k] == corners[1][k] && corners[1][k] == corners[2][k])',
     '    (0..3).find(|&k| {\n        (corners[0][k] - corners[1][k]).abs() <= 1e-12\n            && (corners[1][k] - corners[2][k]).abs() <= 1e-12\n    })',
     [MESH, FLUSH]),
    ('only vertices inside the other operand\'s box move', SNAP,
     '            if let Some(Some(at)) = held.get(&key(p[k])) {',
     '            if let Some(Some(at)) = held.get(&key(p[k])).filter(|_| {\n                bounds(&fixed.positions)\n                    .is_some_and(|(lo, hi)| p.cmpge(lo).all() && p.cmple(hi).all())\n            }) {',
     [MESH, FLUSH]),
    # Never worse.
    ('a refused snapped boolean is not cut as given', COMPILER,
     '            )\n            .ok()\n        } else {',
     '            )\n            .map(Some)?\n        } else {', [MESH]),
    ('the fallback cuts the snapped operands again', COMPILER,
     '            None => (self.cut(&subject.mesh, tool, operator, options)?, 0.0),',
     '            None => (\n                self.cut(\n                    settled.subject.as_ref().unwrap_or(&subject.mesh),\n                    settled.tool.as_ref().unwrap_or(tool),\n                    operator,\n                    options,\n                )?,\n                0.0,\n            ),',
     [MESH]),
    ('the fallback reports the snap it did not use', COMPILER,
     '            None => (self.cut(&subject.mesh, tool, operator, options)?, 0.0),',
     '            None => (self.cut(&subject.mesh, tool, operator, options)?, settled.moved),',
     [MESH]),
    ('exact boolean not retried after an end merge', BB_LIB,
     '        Err(_) if session.merged_an_end() => {',
     '        Err(_) if session.merged_an_end() && false => {', [EXACT]),
    ('exact retry merges ends within the tolerance again', BB_REPORT,
     '    if within {\n        tolerance',
     '    if within || true {\n        tolerance', [EXACT]),
    ('exact end merge not noted', SECTION,
     '            report::merged_an_end();',
     '            {}', [EXACT]),
    # Report.
    ('the snap is not reported', COMPILER,
     '        if moved > 0.0 {\n            built.deviation.add(',
     '        if moved > 0.0 && false {\n            built.deviation.add(', [MESH]),
    ('the snap reported smaller than the move', COMPILER,
     '                crate::deviation::DeviationBound::Certified(moved),',
     '                crate::deviation::DeviationBound::Certified(moved * 0.25),', [MESH]),
    ('an inner snap dropped by the outer boolean', COMPILER,
     '        built.deviation.carry_snaps(&subject.deviation);',
     '        let _ = &subject.deviation;', [MESH]),
    ('an inner snap dropped by the measured deviation', COMPILER,
     '                        deviation.carry_snaps(&built.deviation);',
     '                        let _ = &built.deviation;', [MESH]),
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
