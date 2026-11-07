"""Mutation probe for snap rounding on the binary64 grid (#410, ADR 0100).

Each mutant breaks one decision the snap owes its callers (which cells an
edge meets and in what order, which side owns a tie, how far a boundary
moved, how opposite traversals cancel, what is refused and what is charged)
and must turn a test of the exact arrangement red: a unit test of `snap.rs`,
the scenarios of `tests/snap.rs` or its seeded property test. The visible
face provider's own use of the snap is probed by
`probe_visible_faces_mutants.py`.
"""
import probe_lib

SNAP = "crates/algorithms/planar/exact-arrangement/src/snap.rs"
TRACE = "crates/algorithms/planar/exact-arrangement/src/trace.rs"
TESTS = ["-p", "axiolid-exact-arrangement"]

MUTANTS = [
    # Which cells an edge meets, and in what order.
    ('interior cells ignored', SNAP,
     '        route.extend(order.iter().map(|&i| found[i].0));',
     '        route.extend(std::iter::empty::<usize>());', TESTS),
    ('interior order reversed', SNAP,
     '        route.extend(order.iter().map(|&i| found[i].0));',
     '        route.extend(order.iter().rev().map(|&i| found[i].0));', TESTS),
    ('closed parity inverted', SNAP,
     '        closed[axis] = c.abs().to_bits() & 1 == 0;',
     '        closed[axis] = c.abs().to_bits() & 1 == 1;', TESTS),
    ('the lower gap is used for both sides', SNAP,
     '        (true, true) => (value - down, up - value),',
     '        (true, true) => (value - down, value - down),', TESTS),
    ('an open lower end does not win a tie', SNAP,
     '                        Ordering::Equal if candidate.open => lower = candidate,',
     '                        Ordering::Equal if false => lower = candidate,', TESTS),
    ('an open upper end does not win a tie', SNAP,
     '                        Ordering::Equal if candidate.open => upper = candidate,',
     '                        Ordering::Equal if false => upper = candidate,', TESTS),
    # How far a boundary moved.
    ('a flat end contributes its half-width', SNAP,
     '\n            if points.flat(id).is_none() {\n                displacement',
     '\n            if true {\n                displacement', TESTS),
    ('an end in one cell does not contribute', SNAP,
     '                if points.flat(id).is_none() {\n                    displacement',
     '                if false {\n                    displacement', TESTS),
    ('the half-width of an interior cell is ignored', SNAP,
     '            displacement = displacement.max(cells[index].half_width);',
     '            let _ = cells[index].half_width;', TESTS),
    # How opposite traversals cancel.
    ('cancellation disabled', SNAP,
     '                *net.entry(piece).or_default() += if forward { 1 } else { -1 };',
     '                *net.entry(piece).or_default() += 1;', TESTS),
    ('survives is always true', SNAP,
     '                    .any(|(piece, _)| net.get(piece).is_some_and(|&count| count != 0))',
     '                    .any(|(piece, _)| net.get(piece).is_some_and(|&count| count != 0)) || true', TESTS),
    # What is refused.
    ('the fold check is removed', SNAP,
     '                _ => return Err(refused("snap rounding folded a boundary onto itself")),',
     '                _ => {}', TESTS),
    ('the crossing check is removed', SNAP,
     '    if arrangement.crossings != 0 {',
     '    if false && arrangement.crossings != 0 {', TESTS),
    ('check_winding accepts everything', TRACE,
     '        if outers - holes != i64::from(ring.hole) {',
     '        if false && outers - holes != i64::from(ring.hole) {', TESTS),
    ('check_winding counts holes as outers', TRACE,
     '                if container.hole {\n                    holes += 1;',
     '                if false {\n                    holes += 1;', TESTS),
    # What is charged.
    ('hot cells are free', SNAP,
     'pub const HOT_POINT_SCRATCH_BYTES: usize = 256;',
     'pub const HOT_POINT_SCRATCH_BYTES: usize = 0;', TESTS),
    ('fragments are free', SNAP,
     'pub const FRAGMENT_SCRATCH_BYTES: usize = 64;',
     'pub const FRAGMENT_SCRATCH_BYTES: usize = 0;', TESTS),
]

probe_lib.main(MUTANTS, timeout=3600, reject_uncompiled=True)
