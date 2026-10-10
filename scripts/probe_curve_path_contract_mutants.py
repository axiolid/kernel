"""Mutation probe for curve paths on the curve-evaluation contract (#290,
ADR 0082 amendment).

Each mutant breaks one step the `path_*` queries rest on: the contract's
defaults (the plain queries are the outgoing reading; a provider without
paths refuses every query by a typed `UnsupportedInput` naming itself and
`CURVE_PATH_UNSUPPORTED`, reports no convention and claims nothing
exact), the conformance checks that tell a joint read from the wrong side
and a path refused as bad input from a supported reading; in the
reference provider the side passed through to `CompositeBasis::section_on`,
the refused native parameter, the reference up a path frame is built
against and the banked piece's `+Z` check, the exactness claim and the
convention; in `axiolid-evaluate` the conversion of a path into a
composite (sense, placement, placement exactness, the span checked
against its curve, and back); in `axiolid-curve` the path's own
composition (reversal order and sense, placement order, exactness); and
the compiler's `curve_path`. Each must turn
`evaluate/tests/curve_evaluate_path.rs`, `axiolid-curve`'s unit tests or
`compile/tests/station_curve_path.rs` red.

Not probed, because no test can tell them apart: the refusal of an
unknown `SeamSide` or `CurveMeasure` variant (neither enum has one to
pass), and the refusal of an unknown `PathCurve` kind (it has none yet;
#289 adds one).
"""
import pathlib, subprocess, sys

ROOT = pathlib.Path(__file__).resolve().parents[1]
C = "crates/contracts/operations/curve-evaluate/src/contract.rs"
K = "crates/contracts/operations/curve-evaluate/src/conformance.rs"
P = "crates/algorithms/parametric/evaluate/src/provider.rs"
E = "crates/algorithms/parametric/evaluate/src/station/composite.rs"
V = "crates/representations/analytic/curve/src/path.rs"
B = "crates/execution/compile/src/station/basis.rs"
PATH = ["-p", "axiolid-evaluate", "--test", "curve_evaluate_path"]
CURVE = ["-p", "axiolid-curve", "--lib"]
COMP = ["-p", "axiolid-mesh-compile", "--test", "station_curve_path"]

MUTANTS = [
    # The contract's defaults.
    ('plain point reads the incoming side', C,
     '        self.path_point_at_on(path, at, SeamSide::Outgoing)',
     '        self.path_point_at_on(path, at, SeamSide::Incoming)', [PATH]),
    ('plain tangent reads the incoming side', C,
     '        self.path_tangent_at_on(path, at, SeamSide::Outgoing)',
     '        self.path_tangent_at_on(path, at, SeamSide::Incoming)', [PATH]),
    ('plain frame reads the incoming side', C,
     '        self.path_frame_at_on(path, at, SeamSide::Outgoing)',
     '        self.path_frame_at_on(path, at, SeamSide::Incoming)', [PATH]),
    ('default path refusal names no provider', C,
     '        backend: provider.descriptor().id,\n'
     '        operation: Operation::CurveEvaluation,\n'
     '        input: CURVE_PATH_UNSUPPORTED,',
     '        backend: axiolid_contracts::BackendId::new("unnamed"),\n'
     '        operation: Operation::CurveEvaluation,\n'
     '        input: CURVE_PATH_UNSUPPORTED,', [PATH]),
    ('default path refusal is untyped', C,
     '    GeomError::UnsupportedInput {\n'
     '        backend: provider.descriptor().id,\n'
     '        operation: Operation::CurveEvaluation,\n'
     '        input: CURVE_PATH_UNSUPPORTED,\n'
     '    }',
     '    {\n'
     '        let _ = provider;\n'
     '        GeomError::InvalidInput(CURVE_PATH_UNSUPPORTED.into())\n'
     '    }', [PATH]),
    ('default claims every path frame exact', C,
     '        let _ = (path, at, side);\n        false\n',
     '        let _ = (path, at, side);\n        true\n', [PATH]),
    ('default reports a path convention', C,
     '        let _ = path;\n        DistanceConvention::Unsupported',
     '        let _ = path;\n        DistanceConvention::ArcLength3d', [PATH]),
    # The conformance suite.
    ('conformance skips paths', K,
     '    out.extend(check_paths(provider));\n', '', [PATH]),
    ('conformance takes any refusal as unsupported', K,
     '        (Err(a), Err(b), Err(c)) if side_refused(&a) && side_refused(&b) && side_refused(&c) => {',
     '        (Err(_), Err(_), Err(_)) => {', [PATH]),
    ('conformance expects the incoming side to read the outgoing piece', K,
     '                    (SeamSide::Incoming, place.incoming, place.past_arc.0),',
     '                    (SeamSide::Incoming, place.outgoing, place.past_arc.0),', [PATH]),
    # The reference provider.
    ('provider reads every joint outgoing', P,
     '        let section = basis.section_on(distance, side)?;',
     '        let section = basis.section_on(distance, SeamSide::Outgoing)?;', [PATH]),
    ('provider reads a native parameter as a distance', P,
     '            CurveMeasure::Parameter(_) => Err(GeomError::UnsupportedInput {',
     '            CurveMeasure::Parameter(p) if p.is_finite() => Ok(p),\n'
     '            CurveMeasure::Parameter(_) => Err(GeomError::UnsupportedInput {', [PATH]),
    ('provider frames a path against +Z whatever its up', P,
     '        let (basis, distance, section) = Self::path_section(path, at, side)?;\n'
     '        if self.up == Vec3::Z {',
     '        let (basis, distance, section) = Self::path_section(path, at, side)?;\n'
     '        if true {', [PATH]),
    ('provider frames a banked piece against another up', P,
     '        if let StationCurve::Three(Curve3::Banked(_)) = basis.piece_read(distance, side)?.curve {\n'
     '            self.banked_up()?;',
     '        if let StationCurve::Three(Curve3::Banked(_)) = basis.piece_read(distance, side)?.curve {\n'
     '            let _ = self.banked_up();', [PATH]),
    ('provider exactness ignores the side', P,
     '.is_ok_and(|basis| basis.frame_is_exact_at(distance, side))',
     '.is_ok_and(|basis| basis.frame_is_exact_at(distance, SeamSide::Outgoing))', [PATH]),
    ('provider reports arc length for every path', P,
     '            .map_or(DistanceConvention::Unsupported, |basis| basis.convention())',
     '            .map_or(DistanceConvention::Unsupported, |_| DistanceConvention::ArcLength3d)',
     [PATH]),
    # A path read as a composite.
    ('a path piece loses its sense', E,
     '            reversed: piece.reversed,\n            placement: piece.placement,',
     '            reversed: false,\n            placement: piece.placement,', [PATH]),
    ('a path piece loses its placement', E,
     '            reversed: piece.reversed,\n            placement: piece.placement,',
     '            reversed: piece.reversed,\n            placement: None,', [PATH, COMP]),
    ('a path piece\'s placement claimed exact', E,
     '            placement_exact: piece.placement_exact,\n            ..Self::between',
     '            placement_exact: true,\n            ..Self::between', [PATH]),
    ('a path piece\'s span not checked', E,
     '            ..Self::between(curve, piece.start, piece.end)?',
     '            ..Self {\n'
     '                curve,\n'
     '                start: piece.start,\n'
     '                end: piece.end,\n'
     '                reversed: false,\n'
     '                placement: None,\n'
     '                placement_exact: true,\n'
     '            }', [PATH]),
    ('a station piece hands out no sense', E,
     '        out.reversed = piece.reversed;\n', '', [PATH, COMP]),
    ('a station piece hands out no placement', E,
     '        out.placement = piece.placement;\n', '', [PATH, COMP]),
    # The path's own composition.
    ('a reversed path keeps its order', V,
     '            .rev()\n            .map(PathPiece::reversed)',
     '            .map(PathPiece::reversed)', [CURVE, PATH]),
    ('a reversed path keeps its pieces\' sense', V,
     '            .rev()\n            .map(PathPiece::reversed)',
     '            .rev()', [CURVE, PATH]),
    ('placements compose innermost last', V,
     '            Some(inner) => rigid * inner,',
     '            Some(inner) => inner * rigid,', [CURVE]),
    ('a placed line claimed exact whatever its placement', V,
     '        self.curve.is_line() && self.placement_exact',
     '        self.curve.is_line()', [CURVE, PATH]),
    ('a 3D line not exact', V,
     '            Self::Two(Curve2::Line(_)) | Self::Three(Curve3::Line(_))',
     '            Self::Two(Curve2::Line(_))', [CURVE]),
    # The compiler's flattening, handed out.
    ('curve_path hands out each piece reversed', B,
     '        .map(PathPiece::from)\n        .collect())',
     '        .map(|piece| PathPiece::from(piece.reversed()))\n        .collect())', [COMP]),
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
