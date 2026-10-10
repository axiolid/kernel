"""Mutation probe for bounding volumes (#118).

Each mutant decides containment wrongly, drops a support level, stops
rounding outward, understates the error bound, or skips the candidate
orientations that make an oriented box tight, and must turn a test red.
"""
import pathlib

import probe_lib

ROOT = pathlib.Path(__file__).resolve().parents[1]
C = "crates/algorithms/planar/overlay/src/circle.rs"
S = "crates/algorithms/construction/construct/src/bounding.rs"
CIRCLE = ["-p", "axiolid-overlay", "--test", "circle"]
BOUNDS = ["-p", "axiolid-construct", "--test", "bounding_volumes"]

UNROUNDED = [
    ('    let mut radius = (radius_high + centre_error).next_up();', '    let mut radius = radius_low;'),
    ('        radius = radius.max(distance_bounds(p, &at).1);', '        radius = radius.max(0.0);'),
]

MUTANTS = [
    ('circle: incircle sign read without the orientation', C,
     [('            side == Sign::Zero || side == turn', '            side == Sign::Zero || side == Sign::Positive')], CIRCLE),
    ('circle: three-point support never formed', C,
     [('                    support = vec![pi, pj, pk];', '                    support = vec![pi, pj];')], CIRCLE),
    ('circle: radius neither rounded up nor checked', C, UNROUNDED, CIRCLE),
    ('circle: error bound claimed zero', C,
     [('            error: (radius - radius_low).next_up(),', '            error: 0.0,')], CIRCLE),
    ('sphere: circumsphere side read without the denominator', S,
     [('            side == Sign::Zero || den == Sign::Zero || side != den', '            side == Sign::Zero || side == Sign::Negative')], BOUNDS),
    ('sphere: four-point support never formed', S,
     [('                        support = vec![pi, pj, pk, pl];', '                        support = vec![pi, pj, pk];')], BOUNDS),
    ('sphere: diametral test inverted', S,
     [('            }) != Sign::Positive', '            }) == Sign::Positive')], BOUNDS),
    ('sphere: triangle centre off by a factor two', S,
     [('    (n, ww.add(&ww))', '    (n, ww)')], BOUNDS),
    ('sphere: radius neither rounded up nor checked', S, UNROUNDED, BOUNDS),
    ('box: candidates never ranked, axis-aligned returned', S,
     [('        if volume < best_volume {', '        if false && volume < best_volume {')], BOUNDS),
    ('box: hull face normals not tried', S,
     [('                normals.push((b - a).cross(c - a));', '                let _ = (b - a).cross(c - a);')], BOUNDS),
    ('box: extents rounded to nearest, not up', S,
     [('    while e(x).sub(d).sign() == Some(Sign::Negative) {', '    while false && e(x).sub(d).sign() == Some(Sign::Negative) {')], BOUNDS),
    ('box: upper extremes not tracked', S,
     [('            } else if is_less(&high[k], &s) {', '            } else if false && is_less(&high[k], &s) {')], BOUNDS),
]

probe_lib.main(MUTANTS, timeout=1200)
