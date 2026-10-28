"""Mutation probe for fixed-reference sweeps along chains of segments and
arcs (#122).

Each mutant lets the frame twist unpaid or unreported, loses the re-seed
check at a joint, makes a loop of a chain that merely touches itself,
turns the profile half a turn about the tangent, lets the
reference come within the tangent's line, measures that clearance against
the direction of travel only, loosens the corner tolerance or drops the
corner refusal, drops the fold check or its slack, carries one piece's
frame into the next instead of restarting the law on it, leaves a loop
open, claims no wall term for a chain, or drops the fold check a lone
arc's swept disk makes, and must turn `fixed_reference_chain` red (the
plan's own wall term, the unit tests beside it). The test samples the
exact swept surfaces against the compiled meshes and checks the refusals
by name.
"""
import pathlib

import probe_lib

ROOT = pathlib.Path(__file__).resolve().parents[1]
FR = "crates/algorithms/construction/construct/src/pipe/fixed_reference.rs"
C = "crates/execution/compile/src/compiler.rs"
CHAIN = ["-p", "axiolid-mesh-compile", "--test", "fixed_reference_chain"]
# The compiler reads the walls term itself (`deviation/paths.rs`), so the
# chain's own is checked where it is planned.
PLAN = ["-p", "axiolid-construct", "--lib", "fixed_reference"]

MUTANTS = [
    ('twisting arcs not left to the sampled sweep', FR, '        if term > TWIST_SHARE * chord {', '        if false {', [CHAIN]),
    ('twist term dropped from the walls', FR, '        0.5 * chord + twist', '        0.5 * chord', [CHAIN, PLAN]),
    ('joint re-seed left unchecked', FR, '        if shift > JOINT_SHARE * chord {', '        if false {', [CHAIN, PLAN]),
    ('a corner closing an undeclared chain made a loop', FR, '            && turn(last.end_tangent(), first.start_tangent()) <= allowed);', '            );', [CHAIN]),
    ('seed negated', FR, '        let seed = Frame::from_reference(piece.start_point(), piece.start_tangent(), reference)?.x;', '        let seed = -Frame::from_reference(piece.start_point(), piece.start_tangent(), reference)?.x;', [CHAIN]),
    ('reference negated in the plan', FR, '        reference: r,', '        reference: -r,', [CHAIN]),
    ('clearance check disabled', FR, '        if clearance <= angular {', '        if false {', [CHAIN]),
    ('clearance against the direction of travel only', FR, 'r.cross(t).length().atan2(r.dot(t).abs())', 'r.cross(t).length().atan2(r.dot(t))', [CHAIN]),
    ('corner tolerance at a quarter of the reach', FR, '    let allowed = joint_tolerance(chord, reach);', '    let allowed = joint_tolerance(chord, 0.25 * reach);', [CHAIN]),
    ('corner refusal removed', FR, '        if theta > allowed {', '        if false {', [CHAIN]),
    ('fold check removed', FR, '            if depth >= bend {', '            if false {', [CHAIN]),
    ('fold check without its slack', FR, '                + share;\n', '                + 0.0;\n', [CHAIN]),
    ('frame carried across pieces', FR, '        let seed = Frame::from_reference(piece.start_point(), piece.start_tangent(), reference)?.x;', '        let first = &chain.pieces[0];\n        let seed = Frame::from_reference(first.start_point(), first.start_tangent(), reference)?.x;', [CHAIN]),
    ('loop lofted open', FR, '    let mut mesh = loft::loft(rings, &stations, chain.closed())?;', '    let mut mesh = loft::loft(rings, &stations, false)?;', [CHAIN]),
    ('no wall term for a chain', FR, '    } else {\n        0.5 * chord + twist\n    };', '    } else {\n        twist\n    };', [PLAN]),
    ('lone arc swept disk not checked for a fold', C, '                            axiolid_construct::pipe::check_disk_along(arc, *radius)?;', '                            let _ = (arc, radius);', [CHAIN]),
]

probe_lib.main(MUTANTS, timeout=1200)
