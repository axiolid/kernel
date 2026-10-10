"""Mutation probe for planar region detection (#131).

Each mutant weakens the certificate or the growth -- keeps a region whose
bound exceeds the tolerance, reports a bound below the actual distances,
drops the angle test or the rounding allowance, or calls noisy regions
coplanar -- and must turn a test red.

Equivalent mutant, deliberately not listed: rounding the final square
root or quotient to nearest instead of up. It moves the bound by an ulp,
which no test distance can resolve; it is kept for the proof.

Also not listed: growing without the distance test. Peeling enforces the
certificate afterwards regardless, so results stay correct; the test
during growth keeps regions from overgrowing and peeling back, which is
quality and cost, not soundness.
"""
import pathlib

import probe_lib

ROOT = pathlib.Path(__file__).resolve().parents[1]
P = "crates/algorithms/query/inspect/src/planes.rs"
TESTS = ["-p", "axiolid-inspect", "--test", "planes"]

MUTANTS = [
    ('regions never peeled', P, '            if bound <= distance + fit.slack || members.len() == 1 {', '            if true {', TESTS),
    ('bound from one corner only', P, '        t.iter()\n            .map(|v| {', '        t.iter()\n            .take(1)\n            .map(|v| {', TESTS),
    ('angle not tested', P, '                if unit.dot(fit.normal) < cos_limit {', '                if false && unit.dot(fit.normal) < cos_limit {', TESTS),
    ('no rounding allowance', P, '            slack: 32.0 * f64::EPSILON * far,', '            slack: 0.0,', TESTS),
    ('every region called coplanar', P, '                let coplanar = coplanar(members.iter().flat_map(|&t| tri[t]));', '                let coplanar = true;', TESTS),
]

probe_lib.main(MUTANTS, timeout=1200)
