"""Mutation probe for the seam side on the curve-evaluation contract
(#286, ADR 0082 amendment) and the axes of `frame_at` (#242).

Each mutant breaks one step the sided queries rest on: the contract's
defaults (every side but `Outgoing` refused by a typed `UnsupportedInput`
naming the provider, `Outgoing` delegated to the side-less query), the
conformance checks that tell a silent outgoing frame and a mis-typed
refusal from a supported side, and in the reference provider the side
passed through to `station_section3_on`, the seam predicate it asks, the
measures it locates (an elevated or banked curve's native parameter, a
polyline's refused), the reference up a sided frame is built against, the
banked curve's `+Z` check, and the frame layout `x` tangent, `y` up, `z`
right. Each must turn `evaluate/tests/curve_evaluate_seam_side.rs` or
`evaluate/tests/curve_evaluate_contract.rs` red.

Not probed, because no test can tell them apart: the refusal of an
unknown `SeamSide` (the enum has no third variant to pass), and the
skipped seam lookup on a family without sided seams (an optimisation:
the lookup returns no seam there either).
"""
import pathlib, subprocess, sys

ROOT = pathlib.Path(__file__).resolve().parents[1]
C = "crates/contracts/operations/curve-evaluate/src/contract.rs"
K = "crates/contracts/operations/curve-evaluate/src/conformance.rs"
P = "crates/algorithms/parametric/evaluate/src/provider.rs"
SIDE = ["-p", "axiolid-evaluate", "--test", "curve_evaluate_seam_side"]
CONTRACT = ["-p", "axiolid-evaluate", "--test", "curve_evaluate_contract"]

MUTANTS = [
    ('default point answers Incoming with the outgoing point', C,
     '            SeamSide::Outgoing => self.point_at(curve, at),\n'
     '            _ => Err(side_unsupported(self)),',
     '            SeamSide::Outgoing => self.point_at(curve, at),\n'
     '            _ => self.point_at(curve, at),', [SIDE]),
    ('default tangent answers Incoming with the outgoing tangent', C,
     '            SeamSide::Outgoing => self.tangent_at(curve, at),\n'
     '            _ => Err(side_unsupported(self)),',
     '            SeamSide::Outgoing => self.tangent_at(curve, at),\n'
     '            _ => self.tangent_at(curve, at),', [SIDE]),
    ('default frame answers Incoming with the outgoing frame', C,
     '            SeamSide::Outgoing => self.frame_at(curve, at),\n'
     '            _ => Err(side_unsupported(self)),',
     '            SeamSide::Outgoing => self.frame_at(curve, at),\n'
     '            _ => self.frame_at(curve, at),', [SIDE]),
    ('default refuses Outgoing too', C,
     '            SeamSide::Outgoing => self.frame_at(curve, at),',
     '            SeamSide::Outgoing => Err(side_unsupported(self)),', [SIDE]),
    ('default refusal names no provider', C,
     '        backend: provider.descriptor().id,',
     '        backend: axiolid_contracts::BackendId::new("unnamed"),', [SIDE]),
    ('default refusal is untyped', C,
     '    GeomError::UnsupportedInput {\n'
     '        backend: provider.descriptor().id,\n'
     '        operation: Operation::CurveEvaluation,\n'
     '        input: SEAM_SIDE_UNSUPPORTED,\n'
     '    }',
     '    {\n'
     '        let _ = provider;\n'
     '        GeomError::InvalidInput(SEAM_SIDE_UNSUPPORTED.into())\n'
     '    }', [SIDE]),
    ('conformance skips the sides at seams', K,
     '    out.extend(check_sides_at_seams(provider));\n', '', [SIDE]),
    ('conformance takes any refusal as unsupported', K,
     '    matches!(\n'
     '        error,\n'
     '        GeomError::Unsupported { .. } | GeomError::UnsupportedInput { .. }\n'
     '    )',
     '    let _ = error;\n    true', [SIDE]),
    ('provider reads the outgoing side always', P,
     '        station_section3_on(curve, distance, side).map(Some)',
     '        station_section3_on(curve, distance, SeamSide::Outgoing).map(Some)', [SIDE]),
    ('provider never finds a seam', P,
     '        if !distance.is_finite() || !station_on_seam3(curve, distance)? {',
     '        if !distance.is_finite() || true {', [SIDE]),
    ('provider reads every measure as a station', P,
     '        if !distance.is_finite() || !station_on_seam3(curve, distance)? {',
     '        if !distance.is_finite() {', [SIDE]),
    ('an elevated or banked parameter not located', P,
     '                if matches!(curve, Curve3::Elevated(_) | Curve3::Banked(_)) =>',
     '                if false =>', [SIDE]),
    ('a polyline parameter answered on the incoming side', P,
     '                if side == SeamSide::Incoming',
     '                if false', [SIDE]),
    ('a sided frame ignores the reference up', P,
     '        if self.up == Vec3::Z {\n            return Ok(section.frame());',
     '        if true {\n            return Ok(section.frame());', [SIDE]),
    ('a sided banked frame under another up', P,
     '            self.banked_up()?;\n            return Ok(section.frame());',
     '            return Ok(section.frame());', [SIDE]),
    ('frame axes laid out z up (#242)', P,
     '            y: up,\n            z: right,',
     '            y: -right,\n            z: up,', [SIDE, CONTRACT]),
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
