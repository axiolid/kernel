"""Mutation probe for offset curves as station bases (#289, ADR 0082
amendment).

Each mutant breaks one step a station along an offset rests on: in
`axiolid-evaluate`'s offset module the displacement each law gives (the
planar side, the 3D `V x T` direction, the linear interpolation, the plan
frame), the closed forms (taken at all, the offset circle's radius, the
own measure's inverse), the numerical reading (the plan-measured speed and
rate, the one-sided differences at a piece's ends), the refusals (cusp,
collapse, self-crossing), the banked frame, and the splitting of a base at
its seams with the law carried across pieces; in the composite the offset
piece's convention, exactness, measure exactness, its conversion to a path
piece and the named corner refusal; in `axiolid-curve` an offset of a line
being a line and the inner placement's exactness; and in the compiler's
flattening the planar and 3D laws, the stations' order, the plan frame and
the refused trim. Each must turn `evaluate/tests/station_offset.rs`,
`axiolid-curve`'s unit tests or `compile/tests/station_offset.rs` red.

Not probed, because no test can tell them apart: the planar frame built
directly where the offset runs level (the reference-up frame equals it
there, rounding aside), the Newton step inside a panel (a bisection
fallback reaches the same root), and the measure tables a composite caches
(an offset piece read without one builds the same table).
"""
import pathlib, subprocess, sys

ROOT = pathlib.Path(__file__).resolve().parents[1]
O = "crates/algorithms/parametric/evaluate/src/station/offset.rs"
E = "crates/algorithms/parametric/evaluate/src/station/composite.rs"
V = "crates/representations/analytic/curve/src/path.rs"
B = "crates/execution/compile/src/station/basis.rs"
EVAL = ["-p", "axiolid-evaluate", "--test", "station_offset"]
CURVE = ["-p", "axiolid-curve", "--lib"]
COMP = ["-p", "axiolid-mesh-compile", "--test", "station_offset"]

MUTANTS = [
    # The laws.
    ('a planar offset lies to the right', O,
     'OffsetLaw::Planar { distance } => Ok(distance * section.lateral),',
     'OffsetLaw::Planar { distance } => Ok(-distance * section.lateral),', [EVAL]),
    ('a 3D offset lies along T x V', O,
     '                let normal = reference.cross(section.tangent);',
     '                let normal = section.tangent.cross(reference);', [EVAL, COMP]),
    ('an offset by distances is not interpolated', O,
     '                let at = start.lerp(end, (v / length).clamp(0.0, 1.0));',
     '                let at = start;', [EVAL]),
    ('the plan frame is the section frame', O,
     '                    OffsetFrame::Plan => section.plan()?,',
     '                    OffsetFrame::Plan => *section,', [EVAL]),
    # Closed forms.
    ('no closed form', O,
     '        if !self.law.is_constant() {\n            return Ok(None);\n        }',
     '        if true {\n            return Ok(None);\n        }', [EVAL]),
    ('an offset circle grows by the displacement towards its centre', O,
     '        let scale = (radius - self.displacement(&section, 0.0)?.dot(inward)) / radius;',
     '        let scale = (radius + self.displacement(&section, 0.0)?.dot(inward)) / radius;',
     [EVAL]),
    ('a closed form inverts its scale the wrong way', O,
     '            return Ok((m / scale).clamp(0.0, length));',
     '            return Ok((m * scale).clamp(0.0, length));', [EVAL]),
    # The numerical reading.
    ('a plan-measured offset is measured in 3D', O,
     '        let (along, forward) = if self.plan_measured() {',
     '        let (along, forward) = if false {', [EVAL]),
    ('a plan-measured base runs at unit rate', O,
     '            1.0 / horizontal', '            1.0', [EVAL]),
    ('central differences read past a piece\'s start', O,
     '        let derivative = if v - 2.0 * h < 0.0 {',
     '        let derivative = if false {', [EVAL]),
    # Refusals.
    ('no cusp refused', O,
     '        if !(forward.is_finite() && forward > CUSP_TOLERANCE) {',
     '        if !forward.is_finite() {', [EVAL]),
    ('no collapse refused', O,
     '        if !(scale.is_finite() && scale > CUSP_TOLERANCE) {',
     '        if !scale.is_finite() {', [EVAL]),
    ('no self-crossing refused', O,
     '        self.refuse_self_crossing(&table)?;\n', '', [EVAL]),
    # The banked frame.
    ('a banked base\'s roll is dropped', O,
     '    if banked {\n', '    if false {\n', [EVAL]),
    # Splitting a base.
    ('a base is not split at its seams', O,
     '        cuts.extend(piece.seams()?.into_iter().map(|seam| seam.distance));\n', '',
     [EVAL]),
    ('a law restarts on every base piece', O,
     '                    start: start.lerp(end, (run + a) / total),\n'
     '                    end: start.lerp(end, (run + b) / total),',
     '                    start: start.lerp(end, a / total),\n'
     '                    end: start.lerp(end, b / total),', [EVAL]),
    # The composite.
    ('an offset is measured in arc length whatever its base', E,
     '            Self::Offset(offset) => offset.convention(),',
     '            Self::Offset(_) => DistanceConvention::ArcLength3d,', [EVAL]),
    ('no offset frame is exact', E,
     '            Self::Offset(offset) => offset.frame_is_exact(),',
     '            Self::Offset(_) => false,', [EVAL, COMP]),
    ('every offset\'s measure is exact', E,
     '            Self::Offset(offset) => offset.measure_is_exact(),',
     '            Self::Offset(_) => true,', [COMP]),
    ('an offset piece hands out its base reversed', E,
     '                PathPiece::from(offset.base()),',
     '                PathPiece::from(offset.base()).reversed(),', [EVAL]),
    ('an offset across a corner is a plain gap', E,
     '    if [before, after]\n        .iter()\n'
     '        .any(|piece| matches!(piece.curve, StationCurve::Offset(_)))',
     '    if false', [EVAL, COMP]),
    # The neutral path.
    ('an offset of a line is no line', V,
     '            Self::Offset(offset) => offset.law.is_constant() && offset.base.curve.is_line(),',
     '            Self::Offset(_) => false,', [CURVE]),
    ('an offset\'s base placement is not read', V,
     '                offset.base.placement_exact && offset.base.curve.inner_placement_exact()',
     '                offset.base.curve.inner_placement_exact()', [CURVE]),
    # The compiler's flattening.
    ('a planar offset flattened to the right', B,
     '                None => OffsetLaw::Planar {\n                    distance: *distance,',
     '                None => OffsetLaw::Planar {\n                    distance: -*distance,',
     [COMP]),
    ('a 3D offset flattened along -V', B,
     '                    reference_direction: *direction,',
     '                    reference_direction: -*direction,', [COMP]),
    ('stations flattened from the later one', B,
     '            start: offsets(&pair[0]),', '            start: offsets(&pair[1]),', [COMP]),
    ('a plan run flattened in the section frame', B,
     '        StationFrame::Plan => OffsetFrame::Plan,',
     '        StationFrame::Plan => OffsetFrame::Section,', [COMP]),
    ('a trim of an offset read by its own length', B,
     '                    CurveRelation::Offset { .. } | CurveRelation::OffsetByStations { .. }',
     '                    CurveRelation::ParameterCurve { .. }', [COMP]),
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
