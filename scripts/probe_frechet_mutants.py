"""Mutation probe for Fréchet distance between polylines (#147).

Each mutant drops a family of critical values, lets the walk run backwards,
skips a refusal or a coupling constraint, or answers with a bound instead of
the distance, and must turn a test red.
"""
import pathlib

import probe_lib

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
    ('the first row need not start at the corner', F, '        left.push(if open && span.lo == 0.0 {', '        left.push(if !span.is_empty() {', TESTS),
    ('the first column need not start at the corner', F, '        let mut bottom = if bottom_open && span.lo == 0.0 {', '        let mut bottom = if !span.is_empty() {', TESTS),
]

probe_lib.main(MUTANTS, timeout=1200)
