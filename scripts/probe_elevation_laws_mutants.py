"""Mutation probe for vertical circular arcs and intrinsic (clothoid)
elevation laws (#238).

Each mutant reads the wrong circle, loses the stable form, mis-rebases a
piece, reads an intrinsic profile on a branch where plan distance runs
backwards, drops a term of the inversion, or understates a chord bound,
and must turn a test red.

Not probed: dropping the first-order residual correction of the inverted
height. The residual is at most 1e-12 * max(1, d), so the correction moves
a height by less than the tests can resolve; it is an equivalent mutant at
double precision, kept because it is free and states the contract.
"""
import pathlib, subprocess, sys

ROOT = pathlib.Path(__file__).resolve().parents[1]
C = "crates/representations/analytic/curve/src/elevation.rs"
E = "crates/algorithms/parametric/evaluate/src/elevation.rs"
A = "crates/algorithms/parametric/evaluate/src/arc_length.rs"
B = "crates/algorithms/parametric/evaluate/src/banked.rs"
CURVE = ["-p", "axiolid-curve", "--test", "elevation"]
EVAL = ["-p", "axiolid-evaluate", "--test", "elevation_laws"]

MUTANTS = [
    ('arc height by the cancelling form', C,
     'Some(height + distance * (arc.sin + arc.sin0) / (arc.cos0 + arc.cos))',
     'Some(height + radius * (arc.cos0 - arc.cos))', CURVE),
    ('arc turns the wrong way', C,
     '        let sin = sin0 + distance / radius;',
     '        let sin = sin0 - distance / radius;', CURVE),
    ('arc grade from the start angle', C,
     'Some(arc.sin / arc.cos).filter(|g| g.is_finite())',
     'Some(arc.sin0 / arc.cos).filter(|g| g.is_finite())', CURVE),
    ('zero radius well formed', C,
     '} => height.is_finite() && grade.is_finite() && radius.is_finite() && *radius != 0.0,',
     '} => height.is_finite() && grade.is_finite() && radius.is_finite(),', CURVE),
    ('piece_at not rebased', C,
     '                law.piece_at(local)\n',
     '                law.piece_at(distance)\n', CURVE),
    ('intrinsic height_at guesses the start height', C,
     '            // No closed form: an evaluator integrates it (module docs).\n            Self::Intrinsic { .. } => None,\n        }\n    }\n\n    /// Grade',
     '            // No closed form: an evaluator integrates it (module docs).\n            Self::Intrinsic { height, .. } => Some(*height),\n        }\n    }\n\n    /// Grade', CURVE),
    ('inversion drops the start angle', E,
     '        let theta = theta0 + heading(&curve, s)?;',
     '        let theta = heading(&curve, s)?;', EVAL),
    ('below-vertical certificate disabled', E,
     '    let vertical = || invalid(',
     '    if lo <= hi {\n        return Ok(1.0);\n    }\n    let vertical = || invalid(', EVAL),
    ('negative plan distance accepted', E,
     '    if d < 0.0 {\n        return Err(invalid(\n            "an intrinsic profile is defined forward from its start only",',
     '    if d < -1e300 {\n        return Err(invalid(\n            "an intrinsic profile is defined forward from its start only",', EVAL),
    ('intrinsic chord bound drops 1/cos^3', E,
     'Some(kappa / (floor * floor * floor)).filter(|s| s.is_finite())',
     'Some(kappa).filter(|s| s.is_finite())', EVAL),
    ('arc chord bound drops cos^2', E,
     'Some(1.0 / (radius.abs() * cos * cos * cos)).filter(|s| s.is_finite())',
     'Some(1.0 / (radius.abs() * cos)).filter(|s| s.is_finite())', EVAL),
    ('chord bound across a seam', E,
     '            if lo < lower || hi > upper {',
     '            if lo.is_nan() && lower < upper {', EVAL),
    ('clothoid curvature sup from the start only', E,
     '            start.abs().max((start + rate * span).abs())',
     '            start.abs()', EVAL),
    ('polynomial chord bound skips the quadratic term', E,
     '                    .skip(2)',
     '                    .skip(3)', EVAL),
    ('elevated point reads only closed forms', A,
     '    let height = crate::elevation::elevation_height(&curve.elevation, d)?;',
     '    let height = curve\n        .elevation\n        .height_at(d)\n        .ok_or_else(|| invalid("no height"))?;', EVAL),
    ('banked tangent reads only closed-form grades', B,
     '    let grade = crate::elevation::elevation_grade(&curve.base.elevation, d)?;',
     '    let grade = curve\n        .base\n        .elevation\n        .grade_at(d)\n        .ok_or_else(|| invalid("no grade"))?;', EVAL),
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
