"""Mutation probe for joints near a window (#305, ADR 0082 amendment).

Each mutant breaks one step `path_joint_distances_within` rests on: in
`axiolid-evaluate`'s composite the window restriction itself, the chord
that carries a joint past the window without certifying its piece, the
plan chord of a plan-measured piece, the prefix certified before the
window (its tolerance shares, its sum, an uncertifiable piece within reach
read as its estimate), and the joint's index; the shared interval test
(`DistanceInterval::meets`) and the window's validation; the reference
provider's override; the contract's default (its filter); and the
conformance checks that tell a provider answering joints outside the
window, missing one inside it, certifying one wrongly, disagreeing with
its own whole-path answer or answering a window that is not finite. Each
must turn `evaluate/tests/station_certify.rs` or
`evaluate/tests/curve_evaluate_path.rs` red.

Not probed, because no test can tell them apart: the rounding the chord is
reduced by (`CHORD_ROUNDING`, a hundred times the points' rounding, moves
the floor by about 1e-12 of the coordinates) and the rounding-down of the
floor's sum; the order check of the conformance suite (the reference
provider answers in order, and a provider out of order is also caught by
the certificate check of the joint it misnames).

Every mutant runs in its own process group with a timeout that kills the
whole group, and the file is restored whatever happens.
"""
import pathlib

import probe_lib

ROOT = pathlib.Path(__file__).resolve().parents[1]
TIMEOUT = 1800
K = "crates/algorithms/parametric/evaluate/src/station/composite.rs"
P = "crates/algorithms/parametric/evaluate/src/provider.rs"
I = "crates/contracts/operations/curve-evaluate/src/interval.rs"
T = "crates/contracts/operations/curve-evaluate/src/contract.rs"
F = "crates/contracts/operations/curve-evaluate/src/conformance.rs"

CERT = ["-p", "axiolid-evaluate", "--all-features", "--test", "station_certify"]
PATH = ["-p", "axiolid-evaluate", "--all-features", "--test", "curve_evaluate_path"]

MUTANTS = [
    # The window restriction and what lies past it.
    ('every joint is returned whatever the window', K,
     '            if span.is_none_or(|span| interval.meets(&span)) {',
     '            if span.is_none_or(|span| interval.meets(&span) || true) {', [CERT, PATH]),
    ('no chord carries a joint past the window', K,
     '                if floor > span.upper {\n                    break;',
     '                if false && floor > span.upper {\n                    break;', [CERT, PATH]),
    ('a piece\'s estimate is taken for its floor', K,
     '        Ok((chord - CHORD_ROUNDING * scale).max(0.0))',
     '        Ok(self.length() + 0.0 * (chord - CHORD_ROUNDING * scale))', [CERT, PATH]),
    ('a plan-measured piece is floored by its 3D chord', K,
     '            DistanceConvention::PlanDistance => d.x.hypot(d.y),',
     '            DistanceConvention::PlanDistance => d.length(),', [CERT]),
    # The prefix before the window.
    ('the window shares the tolerance among the pieces before it only', K,
     '        let shares = self.shares(tolerance);\n        let joints',
     '        let shares = match span {\n'
     '            Some(span) => self.shares_until(tolerance, span.upper),\n'
     '            None => self.shares(tolerance),\n'
     '        };\n        let joints', [CERT]),
    ('the prefix is summed from estimates under a window', K,
     '            run = run.plus(piece_length(piece, table.as_ref(), share)?);\n'
     '            let interval = run.interval();',
     '            run = run.plus(if span.is_some() {\n'
     '                DistanceInterval::point(piece.length())\n'
     '            } else {\n'
     '                piece_length(piece, table.as_ref(), share)?\n'
     '            });\n'
     '            let interval = run.interval();', [CERT, PATH]),
    ('an uncertifiable piece within reach is read as its estimate', K,
     '            run = run.plus(piece_length(piece, table.as_ref(), share)?);\n'
     '            let interval = run.interval();',
     '            run = run.plus(piece_length(piece, table.as_ref(), share).or_else(|error| {\n'
     '                if span.is_some() {\n'
     '                    Ok(DistanceInterval::point(piece.length()))\n'
     '                } else {\n'
     '                    Err(error)\n'
     '                }\n'
     '            })?);\n'
     '            let interval = run.interval();', [CERT, PATH]),
    ('a joint is named by the next one\'s index', K,
     '                out.push((k, interval));', '                out.push((k + 1, interval));',
     [CERT, PATH]),
    # The interval test and the window's validation.
    ('meeting a window means lying inside it', I,
     '        self.lower <= other.upper && other.lower <= self.upper',
     '        other.lower <= self.lower && self.upper <= other.upper', [CERT, PATH]),
    ('a window centre that is not finite is accepted', T,
     '    if !around.is_finite() {', '    if false {', [CERT]),
    ('a window that is negative or not finite is accepted', T,
     '    if !(window.is_finite() && window >= 0.0) {', '    if false {', [CERT, PATH]),
    # The provider and the contract.
    ('the reference provider certifies the whole path', P,
     '        CompositeBasis::from_path(path)?.joint_intervals_within(around, window, tolerance)',
     '        let span = DistanceInterval::new(around - window, around + window);\n'
     '        Ok(CompositeBasis::from_path(path)?\n'
     '            .joint_intervals(tolerance)?\n'
     '            .into_iter()\n'
     '            .enumerate()\n'
     '            .filter(|(_, interval)| interval.meets(&span))\n'
     '            .collect())', [CERT, PATH]),
    ('the default keeps every joint', T,
     '            .filter(|(_, interval)| interval.meets(&span))',
     '            .filter(|(_, interval)| interval.meets(&span) || true)', [PATH]),
    # The conformance suite.
    ('the conformance suite accepts a joint outside the window', F,
     '        if !interval.meets(&span) {\n            return Err((',
     '        if false && !interval.meets(&span) {\n            return Err((', [PATH]),
    ('the conformance suite accepts a missing joint', F,
     '        if span.contains(*truth) && !found.iter().any(|(j, _)| *j == joint) {',
     '        if false && span.contains(*truth) && !found.iter().any(|(j, _)| *j == joint) {',
     [PATH]),
    ('the conformance suite accepts an uncertified joint near a window', F,
     '        certifies(*interval, *truth, tolerance).map_err(|detail| {',
     '        Ok::<(), String>(()).map_err(|detail: String| {', [PATH]),
    ('the conformance suite skips agreement with the whole path', F,
     '                        if all.get(*joint).is_some_and(|whole| !whole.meets(interval)) {',
     '                        if all.get(*joint).is_some_and(|whole| false && !whole.meets(interval)) {',
     [PATH]),
    ('the conformance suite accepts a window that is not finite', F,
     '        if provider\n            .path_joint_distances_within(&turns, around, window, tolerance)',
     '        if false\n            && provider\n            .path_joint_distances_within(&turns, around, window, tolerance)',
     [PATH]),
]

probe_lib.main(MUTANTS, timeout=TIMEOUT, reject_uncompiled=True)
