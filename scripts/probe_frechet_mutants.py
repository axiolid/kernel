"""Mutation probe for Fréchet distance between polylines (#147).

Each mutant drops a family of critical values, lets the walk run backwards,
skips a refusal or a coupling constraint, or answers with a bound instead of
the distance, and must turn a test red.
"""
import pathlib, subprocess, sys

ROOT = pathlib.Path(__file__).resolve().parents[1]
F = "crates/algorithms/query/measure/src/frechet.rs"
TESTS = ["-p", "axiolid-measure", "--test", "frechet"]
UNIT = ["-p", "axiolid-measure", "--lib", "frechet"]

MUTANTS = [
    ('equidistant-point critical values dropped', F, '                    if along == 0.0 {', '                    if true || along == 0.0 {', TESTS),
    ('vertex-to-segment critical values dropped', F, '                    if (0.0..=1.0).contains(&t) {\n                        visit(h);', '                    if false && (0.0..=1.0).contains(&t) {\n                        visit(h);', TESTS),
    ('walk may run backwards across a cell (right edge)', F, '                    lo: l.lo.max(right.lo),', '                    lo: right.lo,', TESTS),
    ('walk may run backwards across a cell (top edge)', F, '                    lo: bottom.lo.max(top.lo),', '                    lo: top.lo,', TESTS),
    ('a point need only reach one vertex', F, '        return other.iter().all(|&v| distance(point, v) <= eps);', '        return other.iter().any(|&v| distance(point, v) <= eps);', TESTS),
    ('discrete walk may skip the diagonal predecessor', F, '                _ => above.min(diagonal).min(row[j - 1]),', '                _ => above.min(row[j - 1]),', TESTS),
    ('discrete walk ignores the current pair', F, '            row[j] = here.max(best);', '            row[j] = best.max(here * 0.5);', TESTS),
    ('empty polylines accepted', F, '    if points.is_empty() {\n        return Err(FrechetError::EmptyPolyline);', '    if false {\n        return Err(FrechetError::EmptyPolyline);', TESTS),
    ('non-finite points accepted', F, '    if points.iter().any(|p| !p.is_finite()) {', '    if false && points.iter().any(|p| !p.is_finite()) {', TESTS),
    ('invalid leash accepted', F, '    if !eps.is_finite() || eps < 0.0 {', '    if eps < 0.0 {', TESTS),
    ('answer is the discrete upper bound', F, '    candidates.get(first).copied().unwrap_or(hi)', '    hi', TESTS),
    ('bisection moves the wrong bound', F, '        if decide(a, b, mid) {\n            hi = mid;', '        if !decide(a, b, mid) {\n            hi = mid;', UNIT),
    ('the first row need not start at the corner', F, '        left.push(if open && span.lo == 0.0 { span } else { Span::EMPTY });', '        left.push(if !span.is_empty() { span } else { Span::EMPTY });', TESTS),
    ('the first column need not start at the corner', F, '        let mut bottom = if bottom_open && span.lo == 0.0 {', '        let mut bottom = if !span.is_empty() {', TESTS),
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
