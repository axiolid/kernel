"""Mutation probe for exact mass properties far from the world origin (#343).

Each mutant integrates about the world origin again, measures the solid in
place instead of the moved copy, or drops a term of the shift home, and
must turn a test red: a box and a cylinder on a national grid against their
closed forms, and a wall less flush openings under three roundings of its
length unit.
"""
import probe_lib

E = "crates/algorithms/query/measure/src/exact.rs"
TESTS = [
    ["-p", "axiolid-construct", "--test", "exact_measure"],
    ["-p", "axiolid-mesh-compile", "--test", "site_coordinate_measure"],
]
# The in-place path loses a curved face at site coordinates to the
# quadrature anyway; only the box tells whether its planar fans re-base.
BOX = [["-p", "axiolid-construct", "--test", "exact_measure", "--", "a_box_on_a_national_grid"]]
MOVED = '        brep.transformed(&Transform3::from_translation(-origin))\n            .ok()'

MUTANTS = [
    ('fields integrated about the world origin', E,
     '    let origin = reference_point(brep);', '    let origin = Point3::ZERO;', TESTS),
    ('reference point ignores the vertices', E,
     '    if let Some(first) = vertices.next() {',
     '    if let Some(first) = vertices.next().filter(|_| false) {', TESTS),
    ('solid measured in place, not moved to the reference point', E,
     MOVED, MOVED + '.filter(|_| false)', TESTS),
    ('planar fans in place without the reference point', E,
     [(MOVED, MOVED + '.filter(|_| false)'),
      ('                for point in &mut ring {\n                    *point -= about;\n                }\n', '')],
     BOX),
    ('centroid not shifted home', E,
     '        centroid: home + volume_weighted / volume,',
     '        centroid: volume_weighted / volume,', TESTS),
    ('second moments miss the cross term', E,
     'moments + 2.0 * home * volume_weighted + home * home * volume',
     'moments + home * home * volume', TESTS),
    ('second moments not shifted home', E,
     'moments + 2.0 * home * volume_weighted + home * home * volume',
     'moments', TESTS),
]

probe_lib.main(MUTANTS, timeout=2400)
