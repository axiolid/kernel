"""Mutation probe for sweeps along elevated and banked directrices (#252).

Each mutant drops a half of the composed 3D chord bound, a term of a
derivative bound (the plan's `k^2`, the arc's factor 3, the polynomial's
falling factorial, the pivot's length scaling), a part of the second
derivative, a named break or a grade corner, misreads the plan's domain,
or lets the compiler sweep an elevated curve as one segment, as smooth
across a grade break, or without its end tangents, and must turn a test
red. `tests/elevated.rs` checks the bounds against dense sampling and
finite differences; `tests/elevated_directrix.rs` samples the exact tubes
against the compiled meshes.

Not probed: dropping `elevated_chord_bound`'s own refusal across a break.
It is equivalent: the profile's `elevation_chord_bound` and a chain's
`chain_plan_bounds` each refuse a span across their seams, and
`chord_bound3` refuses across `continuity_breaks3` before it dispatches;
the check states the contract where it is read.
"""
import pathlib, subprocess, sys

ROOT = pathlib.Path(__file__).resolve().parents[1]
EL = "crates/algorithms/parametric/evaluate/src/elevated.rs"
EV = "crates/algorithms/parametric/evaluate/src/elevation.rs"
AL = "crates/algorithms/parametric/evaluate/src/arc_length.rs"
CH = "crates/algorithms/parametric/evaluate/src/chain.rs"
BK = "crates/algorithms/parametric/evaluate/src/banked.rs"
DX = "crates/execution/compile/src/directrix/elevated.rs"
EVAL = ["-p", "axiolid-evaluate", "--test", "elevated"]
SWEEP = ["-p", "axiolid-mesh-compile", "--test", "elevated_directrix"]

MUTANTS = [
    ('chord bound drops the plan', EL,
     '    let bound = plan.hypot(profile) * (1.0 + ROUNDING);',
     '    let bound = profile * (1.0 + ROUNDING);', [EVAL]),
    ('chord bound drops the profile', EL,
     '    let bound = plan.hypot(profile) * (1.0 + ROUNDING);',
     '    let bound = plan * (1.0 + ROUNDING);', [EVAL]),
    ('banked chord bound drops the pivot', EL,
     '    let bound = plan.hypot(profile + pivot) * (1.0 + ROUNDING);',
     '    let bound = plan.hypot(profile) * (1.0 + ROUNDING);', [EVAL]),
    ('plan third derivative drops k^2', EL,
     '            [stretch, stretch * k, stretch * k * k]',
     '            [stretch, stretch * k, 0.0]', [EVAL]),
    ('clothoid third derivative drops k^2', EL,
     '            [stretch, stretch * k, stretch * (rate + k * k)]',
     '            [stretch, stretch * k, stretch * rate]', [EVAL]),
    ('arc profile third derivative without its 3', EV,
     '                3.0 * reach / (r * r * cos.powi(5)),',
     '                reach / (r * r * cos.powi(5)),', [EVAL]),
    ('polynomial bound without the falling factorial', EV,
     '                        let falling: Scalar = (0..k).map(|j| (i - j) as Scalar).product();\n                        falling * c.abs() * m.powi((i - k) as i32)\n                    })\n                    .sum()\n            };\n            [term(1), term(2), term(3)]\n        }\n        ElevationLaw::Piecewise',
     '                        let falling: Scalar = 1.0;\n                        falling * c.abs() * m.powi((i - k) as i32)\n                    })\n                    .sum()\n            };\n            [term(1), term(2), term(3)]\n        }\n        ElevationLaw::Piecewise', [EVAL]),
    ('pivot bend not scaled by L^2', BK,
     '        in_xi[1] / (length * length),',
     '        in_xi[1] / length,', [EVAL]),
    ('second derivative drops the plan', EL,
     '    let out = Vec3::new(plan.x, plan.y, bend);',
     '    let out = Vec3::new(0.0, 0.0, bend);', [EVAL]),
    ('arc profile z\'\' without cos^3', EV,
     '            1.0 / (radius * cos * cos * cos)',
     '            1.0 / radius', [EVAL]),
    ('clothoid plan normal flipped', AL,
     '                Vec2::new(-heading.sin(), heading.cos()) * kappa,',
     '                Vec2::new(heading.sin(), -heading.cos()) * kappa,', [EVAL]),
    ('chain parametric curvature not divided by speed^2', CH,
     '            (c2 - t * t.dot(c2)) / speed2',
     '            c2 - t * t.dot(c2)', [EVAL]),
    ('plan curvature seams not named', EL,
     '            Curve2::Intrinsic(i) => out.extend(i.curvature.seams_within(i.length)),',
     '            Curve2::Intrinsic(_) => {}', [EVAL]),
    ('profile seams not named', EL,
     '    out.extend(crate::elevation::elevation_seams(&curve.elevation));\n    if k >= 2 {',
     '    if k >= 2 {', [EVAL]),
    ('grade corners read one side', EL,
     '        let mut before = crate::elevation::elevation_grade_before(&elevated.elevation, seam)?;',
     '        let mut before = crate::elevation::elevation_grade(&elevated.elevation, seam)?;', [EVAL]),
    ('pivot rate ignored at a corner', EL,
     '            before += rate(true)?;',
     '            before += rate(false)?;', [EVAL]),
    ('line plan given a unit domain', EL,
     '        Curve2::Line(_) => return Some(Scalar::INFINITY),',
     '        Curve2::Line(_) => return Some(1.0),', [EVAL]),
    ('parabola read as a segment', DX,
     '    if coefficients.iter().skip(2).any(|c| *c != 0.0) {',
     '    if coefficients.iter().skip(3).any(|c| *c != 0.0) {', [SWEEP]),
    ('grade break taken as smooth', DX,
     '        .map_or(true, |corners| corners.iter().any(|&c| c > lo && c < hi))',
     '        .map_or(true, |_| false)', [SWEEP]),
    ('elevated end tangents not reported', DX,
     '    Some([unit(start)?, unit(end)?])',
     '    Some([unit(start)?, unit(end)?]).filter(|_| false)', [SWEEP]),
    ('unbounded line plan swept', DX,
     '    if span.start.is_finite() && span.end.is_finite() {',
     '    if true {', [SWEEP]),
]


def run(target):
    return subprocess.run(
        ["cargo", "test", "-q", *target],
        cwd=ROOT, capture_output=True, text=True, timeout=2400,
    ).returncode


survivors = []
for name, rel, old, new, targets in MUTANTS:
    path = ROOT / rel
    original = path.read_text()
    assert original.count(old) == 1, f"anchor for '{name}' not unique/found"
    path.write_text(original.replace(old, new))
    code = 0
    try:
        for target in targets:
            try:
                code = run(target)
            except subprocess.TimeoutExpired:
                code = -1
            if code != 0:
                break
    finally:
        path.write_text(original)
    status = "killed" if code != 0 else "SURVIVED"
    print(f"{status:8} {name}", flush=True)
    if code == 0:
        survivors.append(name)
print(f"{len(MUTANTS) - len(survivors)}/{len(MUTANTS)} killed")
sys.exit(1 if survivors else 0)
