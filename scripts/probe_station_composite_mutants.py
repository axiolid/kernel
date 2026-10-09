"""Mutation probe for curve relations as station bases (#285, ADR 0082
amendment).

Each mutant breaks one step a composite basis rests on: which piece a
joint reads on each side and the tolerance that puts a station on it,
reading a piece from inside at its own ends, a reversed piece (read from
its end, tangent and lateral negated, seam sides swapped), the joint
checks (a gap, an undeclared reversed piece) and the convention check, the
exactness claim (every piece up to the one read, placements included),
the seams (joints, pieces' own seams offset, inexact joints), a trim of a
composite, a closed conic's span across its parameter seam, a line's
measure, the frame a tilted placement carries; and in the compiler the
flattening (segment sense, trim sense, trim parameters, placed pieces
carried). Each must turn a test red: `evaluate/tests/station_composite.rs`
checks composites against closed forms and stations on the pieces;
`compile/tests/station_composite.rs` resolves stations and placements on
composite, nested, trimmed and segmented graph curves against stations on
the pieces and placements composed by hand; `compile/tests/station_placement.rs`
checks tilted placements.

Not probed, because no test can tell them apart: the nesting depth limit,
the 2D point-selector inversion of a trimmed basis (every trim tested is
by parameter), the periodic reading of a closed conic's trim whose sense
disagrees (a quarter arc reads the same either way when its sense agrees),
and the surface curve's 3D curve as a basis (the same flattening as its
curve).
"""
import pathlib, subprocess, sys

ROOT = pathlib.Path(__file__).resolve().parents[1]
E = "crates/algorithms/parametric/evaluate/src/station.rs"
P = "crates/algorithms/parametric/evaluate/src/station/composite.rs"
B = "crates/execution/compile/src/station/basis.rs"
C = "crates/execution/compile/src/station.rs"
EVAL = ["-p", "axiolid-evaluate", "--test", "station_composite"]
COMP = ["-p", "axiolid-mesh-compile", "--test", "station_composite"]
PLACE = ["-p", "axiolid-mesh-compile", "--test", "station_placement"]

MUTANTS = [
    ('joint incoming reads the outgoing piece', P,
     '                    SeamSide::Incoming => (joint - 1, self.pieces[joint - 1].length()),',
     '                    SeamSide::Incoming => (joint, 0.0),', [EVAL, COMP]),
    ('a station on a joint not recognised', P,
     '            if (at - distance).abs() <= slack(distance) {',
     '            if at == distance {', [EVAL]),
    ('a piece read past its start', P,
     '        let side = if u <= slack(u) {',
     '        let side = if false {', [EVAL]),
    ('a piece read past its end', P,
     '        } else if u >= length - slack(length) {',
     '        } else if false {', [EVAL]),
    ('reversed tangent kept', P,
     '            section.tangent = -section.tangent;\n', '', [EVAL, COMP]),
    ('reversed lateral kept', P,
     '            section.lateral = -section.lateral;\n', '', [EVAL, COMP]),
    ('reversed piece read from its start', P,
     '            (self.end - u, swapped)',
     '            (self.start + u, swapped)', [EVAL, COMP]),
    ('reversed seam sides not swapped', P,
     '                SeamSide::Incoming => SeamSide::Outgoing,\n                _ => SeamSide::Incoming,',
     '                SeamSide::Incoming => SeamSide::Incoming,\n                _ => SeamSide::Outgoing,',
     [EVAL]),
    ('a gap accepted', P,
     '    if meets(end, start) {',
     '    if true {', [EVAL, COMP]),
    ('a reversed piece reported as a gap', P,
     '    if reversed {\n        return Err(invalid(format!(',
     '    if false {\n        return Err(invalid(format!(', [EVAL, COMP]),
    ('mixed conventions accepted', P,
     '            if piece.convention() != convention {',
     '            if false {', [EVAL]),
    ('exactness ignores the pieces before', P,
     '        self.pieces[..=piece]',
     '        self.pieces[piece..=piece]', [EVAL]),
    ('joints not listed as seams', P,
     '                out.push(StationSeam::new(at, at, false, exact_before));\n', '',
     [EVAL, COMP]),
    ('a joint after an ellipse claimed exact', P,
     '            exact_before &= piece.curve.measure_is_exact();\n', '', [EVAL]),
    ('a piece\'s seams not offset', P,
     '                let distance = at + seam.distance;',
     '                let distance = seam.distance;', [EVAL]),
    ('a trim of a reversed piece clipped forwards', P,
     '                clipped.start = piece.end - b;\n                clipped.end = piece.end - a;',
     '                clipped.start = piece.start + a;\n                clipped.end = piece.start + b;',
     [EVAL]),
    ('a closed conic not folded past its seam', P,
     '        Ok(if m > turn { m - turn } else { m })',
     '        Ok(m)', [EVAL]),
    ('a 2D line measured by its parameter', P,
     '            Self::Two(Curve2::Line(line)) => t * line.direction.length(),',
     '            Self::Two(Curve2::Line(_)) => t,', [COMP]),
    ('a tilted placement carried as moved', E,
     '        if (rigid.transform_vector3(Vec3::Z) - Vec3::Z).length() <= KEEPS_UP_TOLERANCE {',
     '        if true {', [EVAL, PLACE]),
    ('a tilted plan-measured source accepted', E,
     '        if measure == DistanceConvention::PlanDistance {',
     '        if false {', [EVAL, PLACE]),
    ('a placed atomic basis not carried', B,
     '                    Some(rigid) => section.carried(*rigid, curve.convention())?,',
     '                    Some(_) => section,', [PLACE]),
    ('a placed atomic basis claimed exact', B,
     '                    exact: exact && placement.exact,',
     '                    exact,', [COMP]),
    ('placed pieces not carried', B,
     '                        .map(|piece| piece.placed(placement.transform, placement.exact))',
     '                        .map(|piece| piece)', [COMP]),
    ('placed pieces claimed exact', B,
     '                        .map(|piece| piece.placed(placement.transform, placement.exact))',
     '                        .map(|piece| piece.placed(placement.transform, true))', [COMP]),
    ('a segment\'s sense ignored', B,
     '                    false => reversed(pieces),',
     '                    false => pieces,', [COMP]),
    ('a trim\'s sense ignored', B,
     '            Ok(Flat::Pieces(if *sense_agreement {',
     '            Ok(Flat::Pieces(if true {', [COMP]),
    ('a trim of an atomic curve read whole', B,
     '                    let piece = StationPiece::between_parameters(curve, lo, hi)?;',
     '                    let piece = StationPiece::whole(curve)?;', [COMP]),
    ('a composite frame always claimed exact', C,
     '        basis.frame_is_exact(station.distance, *seam),',
     '        true,', [COMP]),
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
