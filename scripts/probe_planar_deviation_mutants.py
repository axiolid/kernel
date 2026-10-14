"""Mutation probe for an all-planar boolean's deviation bound (#314).

Each mutant breaks one step the near-rounding bound rests on: clipping a
plane's boundary cells to the face (the inward side, the mitred shift,
the corner's turn, whether clipped pieces are used at all), sizing the
slack from the chord bounds the trim flattening reached, and the search
that runs to the sampled maximum when every face is a plane, settled by
the budget only once a deviation is shown (sampled, or a slab width the
report adds). Each must
turn a test red: the `deviation::boolean` unit tests check the clipped
pieces against dense samples of known domains, the slack against dense
samples of curved trims, and the bound against a mesh with a moved
corner; `boolean_deviation.rs` pins the wall less two doors near
rounding and checks every bound against samples of the exact result.
"""
import pathlib

import probe_lib

ROOT = pathlib.Path(__file__).resolve().parents[1]
CERTIFY = "crates/execution/compile/src/certify.rs"
BOOL = "crates/execution/compile/src/deviation/boolean.rs"
MEASURED = ["-p", "axiolid-mesh-compile", "--test", "boolean_deviation"]
UNIT = ["-p", "axiolid-mesh-compile", "--lib", "deviation::boolean"]

MUTANTS = [
    # Clipping.
    ('clip keeps the outer side of an edge', BOOL,
     '            clip_half(polygon, self.segments[k][0], self.inward[k], shift)',
     '            clip_half(polygon, self.segments[k][0], -self.inward[k], shift)',
     [UNIT, MEASURED]),
    ('clip moves the edge inward', BOOL,
     '        let shift = self.slack + self.rounding;',
     '        let shift = -(self.slack + self.rounding);', [UNIT, MEASURED]),
    ('clip ignores the slack', BOOL,
     '        let shift = self.slack + self.rounding;',
     '        let shift = 0.0;', [UNIT, MEASURED]),
    ('clip takes every corner as turning inward', BOOL,
     '                if turn > 1e-9 * (b - v).length() {',
     '                if true {', [UNIT, MEASURED]),
    ('clip takes every corner as turning outward', BOOL,
     '                if turn > 1e-9 * (b - v).length() {',
     '                if false {', [UNIT, MEASURED]),
    ('inward normal parity inverted', BOOL,
     '    if odd {\n        n\n    } else {\n        -n\n    }',
     '    if !odd {\n        n\n    } else {\n        -n\n    }', [UNIT, MEASURED]),
    ('clipped pieces never used', CERTIFY,
     '    if coverage == Coverage::Boundary && [a, b, c] == [0.0; 3] {',
     '    if coverage == Coverage::Boundary && [a, b, c] == [0.0; 3] && false {',
     [UNIT, MEASURED]),
    ('clipped cells bound nothing', CERTIFY,
     '                worst = worst.max(piece_bound(&corners, &candidates, index));',
     '                worst = worst.max(0.0 * piece_bound(&corners, &candidates, index));',
     [UNIT, MEASURED]),
    # Slack sizing.
    ('plane slack ignores the chord bounds', BOOL,
     '        let delta = if affine { chord } else { delta };',
     '        let delta = if affine { 0.0 } else { delta };', [UNIT, MEASURED]),
    ('plane slack is what the flattening allowed', BOOL,
     '        let delta = if affine { chord } else { delta };',
     '        let delta = if affine { delta.max(chord) } else { delta };', [UNIT, MEASURED]),
    ('flattening does not record its chord bounds', BOOL,
     '        *chord = chord.max(d);',
     '        *chord = chord.max(0.0 * d);', [UNIT, MEASURED]),
    # Stopping rule.
    ('all-planar search stops at the budget', BOOL,
     '        certify_affine(&patches, cells, mesh, enough, shown)',
     '        certify_within(&patches, cells, mesh, target, enough)', [UNIT, MEASURED]),
    ('affine search stops at a tenth of the request', CERTIFY,
     '        (index.rounding, index.rounding)',
     '        (1e-4, index.rounding)', [UNIT, MEASURED]),
    # The budget settles an all-planar search only once a deviation is
    # shown, sampled or added by the caller.
    ('budget settles a search at rounding', CERTIFY,
     'lower.max(search.shown) > DEVIATING * stop',
     'lower.max(search.shown) >= 0.0', [UNIT, MEASURED]),
    ('budget never settles an affine search', CERTIFY,
     '(search.affine && top.upper <= enough &&',
     '(false && top.upper <= enough &&', [UNIT, MEASURED]),
    ('a slab width is not shown', CERTIFY,
     'lower.max(search.shown) > DEVIATING * stop',
     'lower > DEVIATING * stop', [UNIT, MEASURED]),
    ('a sampled deviation is not shown', CERTIFY,
     'lower.max(search.shown) > DEVIATING * stop',
     'search.shown > DEVIATING * stop', [UNIT, MEASURED]),
]

probe_lib.main(MUTANTS, timeout=1800)
