"""Mutation probe for arc-length trims and arc-length chains (#239).

Each mutant breaks one step an arc-length position or a chain's chord
bound rests on: the quadrature's acceptance test, a line's speed, the
inverse's end-of-curve refusal, the trim's sense, a piece's rigid
placement or local-frame check, the curvature and chord bounds of a
piece, the joins a flattener splits at, and the model's validation of the
selector. Each must turn a test red: `arc_length_chain.rs` checks
positions and tangents against independent references (the cubic
parabola's binomial series, closed forms, hand-placed chain ends) and
samples every chord bound it asserts; `arc_length_trim.rs` sweeps along
arc-length trims; `trim_arc_length.rs` validates the selector.

Not probed, because no test can tell them apart: the Newton step's sign
(the bracketed secant and bisection still converge), splitting the
quadrature at knots and vertices (adaptive bisection still reaches the
tolerance, only slower), and which piece owns a join (both evaluate the
same point within the local-frame tolerance).
"""
import pathlib, subprocess, sys

ROOT = pathlib.Path(__file__).resolve().parents[1]
Q = "crates/algorithms/parametric/evaluate/src/arc_parameter.rs"
C = "crates/algorithms/parametric/evaluate/src/chain.rs"
B = "crates/algorithms/parametric/evaluate/src/bound.rs"
D = "crates/execution/compile/src/directrix.rs"
V = "crates/representations/modeling/graph/src/validation.rs"
CHAIN = ["-p", "axiolid-evaluate", "--test", "arc_length_chain"]
TRIM = ["-p", "axiolid-mesh-compile", "--test", "arc_length_trim"]
MODEL = ["-p", "axiolid-model", "--test", "trim_arc_length"]

MUTANTS = [
    ('quadrature accepts a single panel', Q,
     '        if (left + right - whole).abs() <= share || !(mid > lo && mid < hi) {',
     '        if true || !(mid > lo && mid < hi) {', [CHAIN]),
    ('backward arc length reported positive', Q,
     '            Ok(if to < from { -magnitude } else { magnitude })',
     '            Ok(magnitude)', [CHAIN]),
    ('line read at unit speed', Q,
     '        Family::Proportional(rate) => return Ok(from + length / rate),',
     '        Family::Proportional(_) => return Ok(from + length),', [CHAIN, TRIM]),
    ('a length past the curve end is clamped, not refused', Q,
     '            if available < target - tolerance {',
     '            if false {', [CHAIN, TRIM]),
    ('placement turns the wrong way', C,
     '            self.x.x * v.x - self.x.y * v.y,',
     '            self.x.x * v.x + self.x.y * v.y,', [CHAIN]),
    ('pieces not moved to the previous end', C,
     '            origin: placement.place(end),',
     '            origin: end,', [CHAIN]),
    ('piece start off its local origin accepted', C,
     '    if offset > PIECE_FRAME_TOLERANCE * length.max(1.0) || angle > PIECE_FRAME_TOLERANCE {',
     '    if angle > PIECE_FRAME_TOLERANCE {', [CHAIN]),
    ('piece start tangent off +x accepted', C,
     '    if offset > PIECE_FRAME_TOLERANCE * length.max(1.0) || angle > PIECE_FRAME_TOLERANCE {',
     '    if offset > PIECE_FRAME_TOLERANCE * length.max(1.0) {', [CHAIN]),
    ('intrinsic chord factor halved', C,
     '            (local_hi - local_lo).powi(2) * 0.125 * k',
     '            (local_hi - local_lo).powi(2) * 0.0625 * k', [CHAIN]),
    ('curvature bound read at the span start', C,
     '    let reach = lo.abs().max(hi.abs());',
     '    let reach = lo.abs();', [CHAIN]),
    ('chord bound spans a join', C,
     '    if local_hi > length + ARC_LENGTH_TOLERANCE * length.max(1.0) {\n        return None;\n    }',
     '    if false {\n        return None;\n    }', [CHAIN]),
    ('a cornered piece certified', C,
     '                    && crate::bound::continuity_breaks2(curve, 1).is_empty()',
     '', [CHAIN]),
    ('chain joins not named as breaks', B,
     '        Curve2::Chain(c) if k >= 1 => c.joins().unwrap_or_default(),',
     '', [CHAIN]),
    ('trim sense ignored', D,
     '            let signed = if sense { *length } else { -*length };',
     '            let signed = *length;', [TRIM]),
    ('arc length on a relation basis read as a parameter', D,
     '                parameter_kind(selector, None, true, label)',
     '                match selector {\n                    TrimSelector::ArcLength(v) => Ok(Some(*v)),\n'
     '                    _ => parameter_kind(selector, None, true, label),\n                }',
     [TRIM]),
    ('non-finite arc length validated', V,
     '        TrimSelector::ArcLength(value) => value.is_finite(),',
     '        TrimSelector::ArcLength(_) => true,', [MODEL]),
    ('arc length not a parameter-kind selector', V,
     '                TrimSelector::Parameter(_) | TrimSelector::ArcLength(_)\n            )\n        }),',
     '                TrimSelector::Parameter(_)\n            )\n        }),', [MODEL]),
]


def run(target):
    return subprocess.run(
        ["cargo", "test", "-q", *target],
        cwd=ROOT, capture_output=True, text=True, timeout=1800,
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
