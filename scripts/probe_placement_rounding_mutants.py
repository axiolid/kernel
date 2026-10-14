"""Mutation probe for the rounding of placing an instance (#295).

Every bound of an `Instance` or `InstanceAtStation` carries the absolute
rounding of applying its `f64` transform (`4 eps M`), named under
`DeviationPath::Placement`; the source is meshed with room for it (and
once more where its own coordinates round by more than the room); a
budget at or below that rounding is refused by name. Each mutant drops
the term, the room, the retry or the refusal, or weakens the derived
constant, and must turn a test red: a torus at `(6e6, 6e6, 100)` over its
budget or missing the term, a refusal that is not by name, or the
constant's unit pin.

Each mutant runs in its own process group under a timeout that kills the
whole group; the mutated file is always restored. A missing or ambiguous
anchor aborts the probe.
"""
import probe_lib

C = "crates/execution/compile/src/compiler.rs"
D = "crates/execution/compile/src/deviation.rs"
P = "crates/execution/compile/src/polyline.rs"
TESTS = [
    ["-p", "axiolid-mesh-compile", "--all-features", "--test", "placement_rounding"],
    ["-p", "axiolid-mesh-compile", "--all-features", "--lib", "deviation::"],
    ["-p", "axiolid-mesh-compile", "--all-features", "--test", "polyline"],
]
TIMEOUT = 1200

MUTANTS = [
    ('placement rounding not added to the bounds', D,
     '    value * factor + rounding\n',
     '    value * factor\n', TESTS),
    ('placement rounding not named in the report', D,
     '        if rounding > 0.0 {\n            placed.add(',
     '        if false {\n            placed.add(', TESTS),
    ('derived constant halved', D,
     'pub const PLACEMENT_EPSILONS: Scalar = 4.0;',
     'pub const PLACEMENT_EPSILONS: Scalar = 2.0;', TESTS),
    ('translated rows read as copies', D,
     '        let copies = offset == 0.0\n            && nonzero <= 1',
     '        let copies = nonzero <= 1', TESTS),
    ('source budget leaves no room for the rounding', C,
     '    let room = placement_room(transform, budget.chord, room)?;',
     '    let room = { placement_room(transform, budget.chord, room)?; 0.0 };', TESTS),
    ('no retry with room for the rounding the source showed', C,
     '        Ok(Some((source_id, transform, room)))',
     '        let _ = (source_id, room);\n        Ok(None)', TESTS),
    ('budget below the rounding not refused', C,
     '    if rounding < chord {\n        return Ok(());',
     '    if true {\n        return Ok(());', TESTS),
    ('floor at the translation scale not checked', C,
     '    refuse_below_rounding(chord, floor)?;',
     '    let _ = refuse_below_rounding;', TESTS),
    ('placed polyline bound without the rounding', P,
     '        let grown = crate::deviation::grown(bound, stretch, rounding);',
     '        let grown = bound * stretch;', TESTS),
    # Dropping only the polyline's initial room is equivalent: the retry
    # restores it and the refusal of the placed points still fires. So the
    # room and the refusal are dropped together, which must be refused by
    # name no longer.
    ('placed polyline refuses nothing below the rounding', P,
     [('    let mut room = crate::compiler::placement_room(transform, chord, 0.0)?;',
       '    let mut room = 0.0;'),
      ('        crate::compiler::refuse_below_rounding(chord, rounding)?;\n',
       '')], TESTS),
    ('station placement without the rounding', C,
     '                let rounding = placed_rounding(&source.mesh, placement.transform, budget.chord)?;',
     '                let rounding = 0.0;', TESTS),
]

probe_lib.main(MUTANTS, timeout=TIMEOUT)
