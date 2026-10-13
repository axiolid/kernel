"""Mutation probe for open curve-path ends (#298, ADR 0082 amendment).

Each mutant breaks one step an open piece end rests on: in `axiolid-curve`
the constructors that leave the end open (`to_end`, `whole`); in
`axiolid-evaluate` its resolution to the curve's end by the evaluator's
own length (the line's parameter domain, an offset's own length rather
than its base's, the estimate itself rather than a number near it) and the
certification of the length it ends, for an offset and for a B-spline; and
the conformance checks that tell a provider reading an open end wrongly,
refusing it as bad input, or certifying its joint elsewhere. Each must turn
`evaluate/tests/station_certify.rs`, `evaluate/tests/curve_evaluate_path.rs`,
`axiolid-curve`'s unit tests or `compile/tests/station_offset.rs` red.

Not probed, because no test can tell them apart: the conformance check
that a distance past an open end is refused (no provider under test reads
past the end of a path it resolves; the reference provider's own refusal
is covered by `station_certify.rs`).
"""
import pathlib

import probe_lib

ROOT = pathlib.Path(__file__).resolve().parents[1]
V = "crates/representations/analytic/curve/src/path.rs"
E = "crates/algorithms/parametric/evaluate/src/station/composite.rs"
C = "crates/algorithms/parametric/evaluate/src/station/certify.rs"
K = "crates/contracts/operations/curve-evaluate/src/conformance.rs"
CERT = ["-p", "axiolid-evaluate", "--test", "station_certify"]
PATH = ["-p", "axiolid-evaluate", "--test", "curve_evaluate_path"]
CURVE = ["-p", "axiolid-curve", "--lib"]
COMP = ["-p", "axiolid-mesh-compile", "--test", "station_offset"]

MUTANTS = [
    # The constructors.
    ('to_end leaves the end stated', V,
     '            open_end: true,\n',
     '            open_end: false,\n', [CURVE, CERT]),
    ('whole starts past the curve\'s start', V,
     '        Self::to_end(curve, 0.0)',
     '        Self::to_end(curve, 1e-3)', [CURVE, CERT, COMP]),
    # The resolution.
    ('an open end is read as stated', E,
     '        let end = if piece.open_end {',
     '        let end = if false && piece.open_end {', [CERT, PATH, COMP]),
    ('an open end stops short of the estimate', E,
     '        let end = if piece.open_end {\n                    curve.end()?',
     '        let end = if piece.open_end {\n                    curve.end()? * (1.0 - 1e-12)',
     [CERT, COMP]),
    ('an open offset ends at its base\'s length', E,
     '    fn end(self) -> GeomResult<Scalar> {\n',
     '    fn end(self) -> GeomResult<Scalar> {\n'
     '        if let Self::Offset(offset) = self {\n'
     '            return Ok(offset.base_length());\n'
     '        }\n', [CERT, PATH, COMP]),
    ('an open line ends at one unit of its measure', E,
     '            Some(speed) => Ok(speed),',
     '            Some(_) => Ok(1.0),', [PATH, CERT]),
    # The certification of the length an open end closes.
    ('an offset\'s certified length stops short of its end', C,
     '            let b = offset.base_measure(piece.end, table)?;',
     '            let b = offset.base_measure(piece.end, table)? * (1.0 - 1e-6);',
     [CERT, PATH]),
    ('a B-spline\'s certified length stops short of its end', C,
     '        curve => curve_span(curve, piece.start, piece.end, tolerance),',
     '        curve => curve_span(curve, piece.start, piece.end * (1.0 - 1e-7), tolerance),',
     [CERT]),
    # The conformance suite.
    ('conformance skips open paths', K,
     '    out.extend(check_open_paths(provider));\n', '', [PATH]),
    ('conformance takes any reading of an open end', K,
     '                if off > 1e-9 {\n'
     '                    out.push(fail(\n'
     '                        "an open end resolves to the curve\'s end",',
     '                if off > 1e9 {\n'
     '                    out.push(fail(\n'
     '                        "an open end resolves to the curve\'s end",', [PATH]),
    ('conformance takes any refusal of an open end as unsupported', K,
     '            Err(error) if side_refused(&error) => {}\n'
     '            Err(error) => out.push(fail(\n'
     '                "an open path is answered or refused as unsupported",',
     '            Err(_) => {}\n'
     '            #[allow(unreachable_patterns)]\n'
     '            Err(error) => out.push(fail(\n'
     '                "an open path is answered or refused as unsupported",', [PATH]),
    ('conformance takes any joint of an open path', K,
     '            } else if let Err(detail) = certifies(intervals[0], slant, JOINT_TOLERANCE) {',
     '            } else if let Err(detail) = certifies(intervals[0], intervals[0].lower, 1.0) {',
     [PATH]),
]

probe_lib.main(MUTANTS, timeout=1800)
