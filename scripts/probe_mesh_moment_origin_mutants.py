"""Mutation probe for mesh second moments far from the world origin (#343).

Each mutant sums about the world origin again or drops a parallel-axis
term of the shift home, and must turn a test red: a box on a national
grid against its closed-form second moments, and the existing closed-form
and parallel-axis tests at the origin.
"""
import probe_lib

M = "crates/algorithms/query/measure/src/mesh.rs"
TESTS = [
    ["-p", "axiolid-measure", "--test", "conditioning"],
    ["-p", "axiolid-measure", "--test", "mass_properties"],
]
SHIFT = '    let moments = moments + 2.0 * base * first + base * base * volume;'

MUTANTS = [
    ('moments summed about the world origin', M,
     '        let [a, b, c] = mesh.triangle(i).map(|j| mesh.position(j as usize) - base);\n'
     '        // Signed volume of the tetrahedron on `base`, times 6.',
     '        let [a, b, c] = mesh.triangle(i).map(|j| mesh.position(j as usize) - base * 0.0);\n'
     '        // Signed volume of the tetrahedron on `base`, times 6.', TESTS),
    ('moments not shifted home', M, SHIFT, '    let moments = moments;', TESTS),
    ('shift misses the cross term', M, SHIFT,
     '    let moments = moments + base * base * volume;', TESTS),
    ('shift misses the volume term', M, SHIFT,
     '    let moments = moments + 2.0 * base * first;', TESTS),
    ('first moment weighted as a third', M,
     '        first += (a + b + c) * (six_v / 24.0);',
     '        first += (a + b + c) * (six_v / 18.0);', TESTS),
    ('volume weighted twice', M,
     '        volume += six_v / 6.0;', '        volume += six_v / 3.0;', TESTS),
]

probe_lib.main(MUTANTS, timeout=1200)
