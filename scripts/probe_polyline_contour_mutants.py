"""Mutation probe for lowering polyline contour segments (#389).

Each mutant skips the interior vertices a polyline segment passes, evaluates
a stored point by the interpolation formula instead of copying it, drops the
domain check, or counts segments instead of ring vertices for the minimum. Others round a
backward start down, leave a seam index unwrapped, or keep a repeated closing
point.
Each must turn a test in `tests/exact_extrusion.rs` or
`tests/profile_validate.rs` red.
"""
import pathlib

import probe_lib

ROOT = pathlib.Path(__file__).resolve().parents[1]
L = "crates/algorithms/construction/construct/src/contour_lower.rs"
EXTRUSION = ["-p", "axiolid-construct", "--test", "exact_extrusion"]
VALIDATE = ["-p", "axiolid-construct", "--test", "profile_validate"]

MUTANTS = [
    ('interior vertices skipped traversing forward', L, '        while k < to {', '        while k < to && false {', EXTRUSION),
    ('interior vertices skipped traversing backward', L, '        while k > to {', '        while k > to && false {', EXTRUSION),
    ('stored point at an integer parameter recomputed', L, '        let i = t.floor();\n        if i == t {\n            return points[i as usize % count];\n        }', '        let mut i = t.floor();\n        if i == t && i > 0.0 {\n            i -= 1.0;\n        } else if i == t {\n            return points[0];\n        }', EXTRUSION),
    ('domain outside the points accepted (extruder)', L, '    if !inside(from) || !inside(to) {', '    if false {', EXTRUSION),
    ('domain outside the points accepted (validator)', L, '    if !inside(from) || !inside(to) {', '    if false {', VALIDATE),
    ('minimum counted in segments (validator)', L, '    if vertices.len() < 2 {', '    if contour.segments.len() < 2 {', VALIDATE),
    ('minimum counted in segments (extruder)', L, '    if vertices.len() < 2 {', '    if contour.segments.len() < 2 {', EXTRUSION),
    ('empty range accepted', L, '    if from == to {\n        return Err(degenerate(', '    if false {\n        return Err(degenerate(', VALIDATE),
    ('backward start rounded down instead of up', L, '        let mut k = from.ceil() - 1.0;', '        let mut k = from.floor() - 1.0;', EXTRUSION),
    ('seam index not wrapped', L, '            return points[i as usize % count];', '            return points[i as usize];', EXTRUSION),
    ('repeated closing point kept', L, '            vertices.pop();', '            let _ = &mut vertices;', EXTRUSION),
    ('repeated point always pushed', L, '        if point.x.to_bits() != last.x.to_bits() || point.y.to_bits() != last.y.to_bits() {', '        if true || point.x.to_bits() != last.x.to_bits() || point.y.to_bits() != last.y.to_bits() {', EXTRUSION),
]

probe_lib.main(MUTANTS, timeout=1200)
