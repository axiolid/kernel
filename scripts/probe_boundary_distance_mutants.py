"""Mutation probe for the certified boundary distance (#125, C18).

Each mutant weakens a guarantee -- a lower bound that could overshoot the
true distance, a witness that is not on the boundary, a pruning rule that
drops the nearest pair -- or the refinement that closes the interval, and
must turn a test red.
"""
import pathlib, subprocess, sys

ROOT = pathlib.Path(__file__).resolve().parents[1]
D = "crates/algorithms/query/measure/src/exact_distance.rs"
M = "crates/algorithms/query/measure/src/exact_domain.rs"
SOLIDS = ["-p", "axiolid-construct", "--test", "boundary_distance"]
UNIT = ["-p", "axiolid-measure", "--features", "exact", "--lib"]
BOTH = (SOLIDS, UNIT)

MUTANTS = [
    ('patch sphere halved', D, '    let radius = 0.5 * ((hi.x - lo.x).abs() * lu + (hi.y - lo.y).abs() * lv);', '    let radius = 0.25 * ((hi.x - lo.x).abs() * lu + (hi.y - lo.y).abs() * lv);', BOTH),
    ('cylinder u-rate ignores radius', D, '        Surface::Cylinder(c) => (c.radius.abs() * frame_scale(&c.frame), c.frame.z.length()),', '        Surface::Cylinder(c) => (frame_scale(&c.frame), c.frame.z.length()),', BOTH),
    ('trig range misses the peak', D, '    if reaches(peak) {\n        hi = amplitude;\n    }', '', BOTH),
    ('projection gap one-sided', D, '        best = best.max(b_lo - a_hi).max(a_lo - b_hi);', '        best = best.max(b_lo - a_hi);', BOTH),
    ('normal cone ignores spread', D, '    let angle = spread + aperture + 1e-9;', '    let angle = aperture + 1e-9;', BOTH),
    ('cone apex pruned', D, "            if !apex.is_finite() || (apex >= lo.y.min(hi.y) - 1e-9 && apex <= lo.y.max(hi.y) + 1e-9) {\n                return None;\n            }", '', UNIT),
    ('witness without certification', D, '        let witness = if inside {', '        let witness = if true {', BOTH),
    ('outside patch kept as inside', D, '                        Some(false) => return Ok(None),', '                        Some(false) => inside = true,', BOTH),
    ('pole winding dropped', M, '            if pole > q.y {\n                total += winding;\n            }', '', UNIT),
    ('touch test trusts boxes', M, '        if inside(arc.a) || inside(arc.b) || depth >= 40 {', '        if true {', SOLIDS),
    ('clearance rounds to a verdict', D, '        } else if self.lower > limit {', '        } else if self.upper >= limit {', SOLIDS),
]

def run(targets):
    if not isinstance(targets, tuple):
        targets = (targets,)
    for target in targets:
        code = subprocess.run(
            ["cargo", "test", "-q", *target],
            cwd=ROOT, capture_output=True, text=True, timeout=1200,
        ).returncode
        if code != 0:
            return code
    return 0

survivors = []
for name, rel, old, new, target in MUTANTS:
    path = ROOT / rel
    original = path.read_text()
    assert original.count(old) == 1, f"anchor for '{name}' not unique/found"
    path.write_text(original.replace(old, new))
    try:
        try:
            code = run(target)
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
