"""Mutation probe for the tapered revolution's welding of axis vertices (#364).

Each mutant drops, weakens or misreads a step of `tapered_welds`: the
refusal of a vertex on the axis at one end only, the refusal of one that
slides along the axis, and the welding of a vertex that stays on it. Each
must turn a test in `tests/sweeps.rs` red.

Not mutated: the slices of the hole flags in `tapered_welds`. A hole cannot
touch the axis, so those flags are all false and any misreading of them
yields the same weld.
"""
import pathlib

import probe_lib

ROOT = pathlib.Path(__file__).resolve().parents[1]
S = "crates/algorithms/construction/construct/src/sweep.rs"
L = "crates/algorithms/construction/construct/src/loft.rs"
TESTS = ["-p", "axiolid-construct", "--test", "sweeps"]

MUTANTS = [
    ('vertex on the axis at one end only welded', S, '        if here != there {', '        if false {', TESTS),
    ('sliding vertex welded', S, '        if here && from[f] != to[f] {', '        if false {', TESTS),
    ('every vertex on the axis refused', S, '        if here && from[f] != to[f] {', '        if here {', TESTS),
    ('nothing welded', S, '        welded.push(here);', '        welded.push(false);', TESTS),
    ('everything welded', S, '        welded.push(here);', '        welded.push(true);', TESTS),
    ('end flags read from the start ring', S, '    let at_end = flags[outer..2 * outer]', '    let at_end = flags[..outer]', TESTS),
    ('one-end refusal with the wrong variant', S, '    let refuse = |input| GeomError::UnsupportedInput {', '    let refuse = |input| GeomError::InvalidInput(input.to_owned());\n    let _unused = |input| GeomError::UnsupportedInput {', TESTS),
    ('weld flags not handed to the loft', S, 'loft::loft_tapered_welded(&far, start, &stations, &welded)?', 'loft::loft_tapered(&far, start, &stations)?', TESTS),
    ('tapered caps do not share the axis vertices', L, '            let end = t.map(|f| id[last + f as usize]);', '            let end = t.map(|f| (last + f as usize) as u32);', TESTS),
]

probe_lib.main(MUTANTS, timeout=1200)
