"""Mutation probe for Hausdorff distance and the certified Frechet decision
(#147, ledger row H4), both in `src/polyline_hausdorff.rs`.

Each mutant drops a family of tie candidates, weakens the evaluation to
something other than the true minimum, skips an endpoint or a refusal,
breaks the two-sided distance's asymmetry, or corrupts the certified
decision's margin or branching -- and must turn a test red.
"""
import pathlib, subprocess, sys

ROOT = pathlib.Path(__file__).resolve().parents[1]
F = "crates/algorithms/query/measure/src/polyline_hausdorff.rs"
TESTS = ["-p", "axiolid-measure", "--test", "hausdorff", "--test", "frechet_certified"]

MUTANTS = [
    (
        'vertex-vertex tie candidates dropped',
        F,
        '                if slope != 0.0 {',
        '                if false && slope != 0.0 {',
        TESTS,
    ),
    (
        'vertex-vs-line tie candidates dropped',
        F,
        '                    if (0.0..=1.0).contains(&t) && foot_active(p0, d, t, u, w) {',
        '                    if false {',
        TESTS,
    ),
    (
        'line-vs-line tie candidates dropped',
        F,
        '                    if (0.0..=1.0).contains(&t)\n                        && foot_active(p0, d, t, u, w)\n                        && foot_active(p0, d, t, u2, w2)\n                    {',
        '                    if false {',
        TESTS,
    ),
    (
        'evaluation takes the farthest feature instead of the nearest',
        F,
        '        .map(|s| point_segment_distance(p, s[0], s[1]))\n        .fold(f64::INFINITY, f64::min)',
        '        .map(|s| point_segment_distance(p, s[0], s[1]))\n        .fold(f64::INFINITY, f64::max)',
        TESTS,
    ),
    (
        "a segment's own endpoints are not among the candidates",
        F,
        '    let mut candidates: Vec<f64> = vec![0.0, 1.0];',
        '    let mut candidates: Vec<f64> = vec![];',
        TESTS,
    ),
    (
        'the two-sided distance takes the smaller direction',
        F,
        '    Ok(a_to_b.max(b_to_a))',
        '    Ok(a_to_b.min(b_to_a))',
        TESTS,
    ),
    (
        'one-sided Hausdorff refusals skipped',
        F,
        'pub fn one_sided_polyline_hausdorff_distance(a: &[Point3], b: &[Point3]) -> Result<f64, FrechetError> {\n    check(a)?;\n    check(b)?;',
        'pub fn one_sided_polyline_hausdorff_distance(a: &[Point3], b: &[Point3]) -> Result<f64, FrechetError> {\n    let _ = check(a);\n    let _ = check(b);',
        TESTS,
    ),
    (
        'certified decision refusals skipped',
        F,
        '    eps: f64,\n) -> Result<FrechetDecision, FrechetError> {\n    check(a)?;\n    check(b)?;',
        '    eps: f64,\n) -> Result<FrechetDecision, FrechetError> {\n    let _ = check(a);\n    let _ = check(b);',
        TESTS,
    ),
    (
        'the AtMost branch grows the leash instead of shrinking it',
        F,
        '    if eps >= margin && decide(a, b, eps - margin) {',
        '    if eps >= margin && decide(a, b, eps + margin) {',
        TESTS,
    ),
    (
        "the MoreThan branch's condition is inverted",
        F,
        '    if !decide(a, b, eps + margin) {',
        '    if decide(a, b, eps + margin) {',
        TESTS,
    ),
    (
        'Undecided is never returned',
        F,
        '    Ok(FrechetDecision::Undecided)',
        '    Ok(FrechetDecision::AtMost)',
        TESTS,
    ),
    (
        'the error margin is always zero',
        F,
        '    K * (r + eps.abs() + 1.0) * f64::EPSILON',
        '    0.0',
        TESTS,
    ),
    (
        'AtMost and MoreThan are confused',
        F,
        '        return Ok(FrechetDecision::AtMost);',
        '        return Ok(FrechetDecision::MoreThan);',
        TESTS,
    ),
]


def run(target):
    return subprocess.run(
        ["cargo", "test", "-q", *target],
        cwd=ROOT, capture_output=True, text=True, timeout=1200,
    ).returncode


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
