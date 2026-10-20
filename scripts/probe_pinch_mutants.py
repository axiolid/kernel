"""Mutation probe for the "result touches itself" check (#194, #328, #327).

`crates/execution/compile/src/pinch.rs` welds a boolean's result by exact
position and finds an edge of more than two faces or a vertex of separate
fans; `cut` in `compiler.rs` refuses such a result by name. Each mutant
makes the check refuse too much (a false positive on a convex chain of
half-space intersections) or too little (a true contact returned), and
must turn a test red: the unit tests in `pinch.rs`,
`tests/half_space_chain_contact.rs` (convex chains through earlier cut
edges at their analytic volume, a seeded sweep never refused, the notched
wall of #328 refused along its corner line), `tests/tangent_void.rs`
(#194's exact pinch) and `tests/curved_wall_openings.rs` (#327's tangent
box, refused).

Not listed: dropping the filter of triangles a weld collapses. A boolean
returns no triangle with two corners at one position, so no input reaches
it.

Runs on `scripts/probe_lib.py`: each mutant in its own process group
under a timeout that kills the group, the source always restored, and a
missing or ambiguous anchor aborts the probe. A mutant that does not
compile counts as BROKEN, not killed.
"""
import pathlib

import probe_lib

ROOT = pathlib.Path(__file__).resolve().parents[1]
P = "crates/execution/compile/src/pinch.rs"
C = "crates/execution/compile/src/compiler.rs"
CRATE = ["-p", "axiolid-mesh-compile"]
UNIT = CRATE + ["--lib", "pinch"]
CHAIN = CRATE + ["--test", "half_space_chain_contact"]
VOID = CRATE + ["--test", "tangent_void"]
CURVED = CRATE + ["--test", "curved_wall_openings"]
ALL = [UNIT, CHAIN, VOID, CURVED]
TIMEOUT = 1800

MUTANTS = [
    # Too much refused.
    ('every shared edge read as a contact', P,
     'sorted.iter().find(|(_, faces)| faces.len() > 2)',
     'sorted.iter().find(|(_, faces)| faces.len() > 1)', ALL),
    ('every vertex read as separate fans', P,
     '        if (1..faces.len()).any(|i| root(&mut parent, i) != first) {',
     '        if (1..faces.len()).any(|i| root(&mut parent, i) == first) {', ALL),
    # Too little refused.
    ('an edge of four faces not read as a contact', P,
     'sorted.iter().find(|(_, faces)| faces.len() > 2)',
     'sorted.iter().find(|(_, faces)| faces.len() > 4)', ALL),
    ('separate fans at a vertex not read as a contact', P,
     '        if (1..faces.len()).any(|i| root(&mut parent, i) != first) {',
     '        if (1..faces.len()).any(|i| root(&mut parent, i) != first && false) {', ALL),
    ('positions not welded', P,
     '        .map(|t| [0, 1, 2].map(|k| welded[t[k] as usize]))',
     '        .map(|t| [t[0], t[1], t[2]])\n        .inspect(|_| { let _ = &welded; })', ALL),
    ('negative zero not welded to zero (#283)', P,
     '            (p.x + 0.0).to_bits(),\n            (p.y + 0.0).to_bits(),\n            (p.z + 0.0).to_bits(),',
     '            p.x.to_bits(),\n            p.y.to_bits(),\n            p.z.to_bits(),', ALL),
    ('a result touching itself returned', C,
     '        if let Some(pinch) = crate::pinch::find(&outcome.mesh) {',
     '        if let Some(pinch) = crate::pinch::find(&outcome.mesh).filter(|_| false) {', ALL),
]


probe_lib.main(MUTANTS, timeout=TIMEOUT, reject_uncompiled=True)
