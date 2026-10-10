"""Mutation probe for the certified boundary distance (#125, C18).

Each mutant weakens a guarantee -- a lower bound that could overshoot the
true distance, a witness that is not on the boundary, a pruning rule that
drops the nearest pair -- or the refinement that closes the interval, and
must turn a test red.
"""
import pathlib

import probe_lib

ROOT = pathlib.Path(__file__).resolve().parents[1]
D = "crates/algorithms/query/measure/src/exact_distance.rs"
M = "crates/algorithms/query/measure/src/exact_domain.rs"
T = "crates/algorithms/query/measure/src/exact_distance/touch.rs"
SOLIDS = ["-p", "axiolid-construct", "--test", "boundary_distance"]
UNIT = ["-p", "axiolid-measure", "--features", "exact", "--lib"]
BOTH = (SOLIDS, UNIT)

MUTANTS = [
    ('patch sphere halved', D, '    let radius = 0.5 * ((hi.x - lo.x).abs() * lu + (hi.y - lo.y).abs() * lv);', '    let radius = 0.25 * ((hi.x - lo.x).abs() * lu + (hi.y - lo.y).abs() * lv);', BOTH),
    ('cylinder u-rate ignores radius', D, '        Surface::Cylinder(c) => (c.radius.abs() * frame_scale(&c.frame), c.frame.z.length()),', '        Surface::Cylinder(c) => (frame_scale(&c.frame), c.frame.z.length()),', BOTH),
    ('trig range misses the peak', D, '    if reaches(peak) {\n        hi = amplitude;\n    }', '', BOTH),
    ('projection gap one-sided', D, '        best = best.max(b_lo - a_hi).max(a_lo - b_hi);', '        best = best.max(b_lo - a_hi);', BOTH),
    ('normal cone ignores spread', D, '    let angle = spread + aperture + 1e-9;', '    let angle = aperture + 1e-9;', BOTH),
    ('cone apex pruned', D, "            if !apex.is_finite() || (apex >= lo.y.min(hi.y) - 1e-9 && apex <= lo.y.max(hi.y) + 1e-9)\n            {\n                return None;\n            }", '', UNIT),
    ('witness without certification', D, '        let witness = if inside {', '        let witness = if true {', BOTH),
    ('outside patch kept as inside', D, '                        Some(false) => return Ok(None),', '                        Some(false) => inside = true,', BOTH),
    ('pole winding dropped', M, '            if pole > q.y {\n                total += winding;\n            }', '', UNIT),
    ('touch test trusts boxes', M, '        if inside(arc.a) || inside(arc.b) || depth >= 40 {', '        if true {', SOLIDS),
    ('touch trusts a nanometre gap', T, '    d <= 1e-12 * (p.length() + q.length())', '    d <= 1e-6 * (p.length() + q.length())', UNIT),
    ('touch off the face trim', T, '        let on_face = inside\n            || match', '        let on_face = true\n            || match', UNIT),
    ('touch misses the crossing', T, '    let candidates = [crossing, Some(t0)', '    let candidates = [None, Some(t0)', UNIT),
    ('touch ignores edge pairs', T, '        (Shape::Edge { .. }, Shape::Edge { .. }) => edge_edge(side_a, a, side_b, b)?,', '        (Shape::Edge { .. }, Shape::Edge { .. }) => None,', UNIT),
    ('clearance rounds to a verdict', D, '        } else if self.lower > limit {', '        } else if self.upper >= limit {', SOLIDS),
]

probe_lib.main(MUTANTS, timeout=1200)
