"""Mutation probe for curve-path piece bounds at a parameter or a point
(#304, ADR 0082 amendment).

Each mutant breaks one step a bounded piece rests on: in `axiolid-curve`
the constructors that set a bound and `is_stated`; in `axiolid-evaluate`
the resolution of a bound to a distance (a start bound, an end bound
after a stated start, the closed conic's fold and its wrap round the seam,
the refusal of an end both open and bounded, two parameters read as
`between_parameters`) and the inversion of a point (the joint tolerance
either way, the closed-form check, Newton on a B-spline, the refusal of
several places and the merging of one place reached twice, on a polyline
and on a B-spline); and the conformance checks that tell a provider
reading a bound elsewhere, refusing it as bad input, or reading a point it
must refuse. Each must turn `evaluate/tests/curve_path_bounds.rs`,
`axiolid-curve`'s unit tests or `compile/tests/station_curve_path.rs` red.

Not probed, because no test can tell them apart: adding a B-spline's knots
to its samples (the tests' splines have no interior knot); the conformance
check's joint certification, which the shifted provider already fails by
its readings.
"""
import probe_lib

V = "crates/representations/analytic/curve/src/path.rs"
E = "crates/algorithms/parametric/evaluate/src/station/composite.rs"
B = "crates/algorithms/parametric/evaluate/src/station/path_bound.rs"
K = "crates/contracts/operations/curve-evaluate/src/conformance.rs"
BOUNDS = ["-p", "axiolid-evaluate", "--test", "curve_path_bounds"]
CURVE = ["-p", "axiolid-curve", "--lib"]
COMP = ["-p", "axiolid-mesh-compile", "--test", "station_curve_path"]

MUTANTS = [
    # The representation.
    ('with_start sets no bound', V,
     '        self.start_bound = Some(bound);',
     '        self.start_bound = None;', [CURVE, BOUNDS]),
    ('with_end keeps an open end', V,
     '        self.end = Scalar::INFINITY;\n        self.open_end = false;',
     '        self.end = Scalar::INFINITY;', [CURVE, BOUNDS]),
    ('a bounded piece counts as stated', V,
     '        !self.open_end && self.start_bound.is_none() && self.end_bound.is_none()',
     '        !self.open_end', [CURVE, BOUNDS]),
    # The resolution.
    ('a start bound is read as the stated start', E,
     '                    Some(t) => curve.measure_at(curve.turned(t))?,',
     '                    Some(_) => piece.start,', [BOUNDS, COMP]),
    ('a start angle is not folded into the first turn', E,
     '                    Some(t) => curve.measure_at(curve.turned(t))?,',
     '                    Some(t) => curve.measure_at(t)?,', [BOUNDS]),
    ('an end after a stated start does not cross the seam', E,
     '        if self.is_closed_conic() && m <= start {',
     '        if false && m <= start {', [BOUNDS]),
    ('two parameters are read as two measures', E,
     '                Self::between_parameters(curve, t0, t1)?',
     '                Self::between(curve, curve.measure_at(t0)?, curve.measure_at(t1)?)?',
     [BOUNDS, COMP]),
    ('a closed conic does not run round', E,
     '        if !self.is_closed_conic() {\n            return (t0, t1);',
     '        if true {\n            return (t0, t1);', [BOUNDS, COMP]),
    ('a whole number of turns spans nothing', E,
     '                0.0 => turn,',
     '                0.0 => 0.0,', [BOUNDS]),
    ('an end both open and bounded is read', E,
     '        if piece.open_end && piece.end_bound.is_some() {',
     '        if false && piece.open_end && piece.end_bound.is_some() {', [BOUNDS]),
    # The inversion of a point.
    ('the point tolerance is a thousand times looser', B,
     '    JOINT_TOLERANCE * point.abs().max_element().max(1.0)',
     '    1e3 * JOINT_TOLERANCE * point.abs().max_element().max(1.0)', [BOUNDS]),
    ('the point tolerance is a thousand times tighter', B,
     '    JOINT_TOLERANCE * point.abs().max_element().max(1.0)',
     '    1e-3 * JOINT_TOLERANCE * point.abs().max_element().max(1.0)', [BOUNDS]),
    ('a closed-form parameter is not checked against the curve', B,
     '                Ok(if self.distance(t, point)? <= tolerance {',
     '                Ok(if true {', [BOUNDS]),
    ('a B-spline candidate is not refined', B,
     '            let refined = self.refine(at[i], a, b, point)?;',
     '            let refined = at[i];', [BOUNDS]),
    ('several places resolve to the first', B,
     '        _ => Err(refused(PATH_POINT_AMBIGUOUS)),',
     '        [t, ..] => Ok(*t),', [BOUNDS]),
    ('one place reached twice is two', B,
     '            Some(place) if joins => {',
     '            Some(place) if false && joins => {', [BOUNDS]),
    ('a polyline\'s places always merge', B,
     '        Ok((first..last)',
     '        Ok(true || (first..last)', [BOUNDS]),
    ('a B-spline\'s places always merge', B,
     '            Ok(between && self.distance(0.5 * (ta + tb), point)? <= tolerance)',
     '            Ok(true || between)', [BOUNDS]),
    # The conformance suite.
    ('conformance skips bounded paths', K,
     '    out.extend(check_bounded_paths(provider));\n', '', [BOUNDS]),
    ('conformance takes any reading of a bound', K,
     '                if off > 1e-9 {\n'
     '                    out.push(fail(\n'
     '                        "a bound resolves to the distance of its parameter or point",',
     '                if off > 1e9 {\n'
     '                    out.push(fail(\n'
     '                        "a bound resolves to the distance of its parameter or point",',
     [BOUNDS]),
    ('conformance takes any refusal of a bound as unsupported', K,
     '            Err(error) if side_refused(&error) => {}\n'
     '            Err(error) => out.push(fail(\n'
     '                "a bounded path is answered or refused as unsupported",',
     '            Err(_) => {}\n'
     '            #[allow(unreachable_patterns)]\n'
     '            Err(error) => out.push(fail(\n'
     '                "a bounded path is answered or refused as unsupported",', [BOUNDS]),
    ('conformance reads points it must refuse', K,
     '        if let Ok(point) = provider.path_point_at(path, CurveMeasure::Distance(0.5)) {',
     '        if let (false, Ok(point)) =\n'
     '            (true, provider.path_point_at(path, CurveMeasure::Distance(0.5)))\n'
     '        {', [BOUNDS]),
]

probe_lib.main(MUTANTS, timeout=1800)
