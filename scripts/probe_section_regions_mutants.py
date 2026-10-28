"""Mutation probe for in-plane section regions (#359, ADR 0091).

Each mutant breaks one decision of `ScalarSection::section_request`: the side
of an in-plane face, the dropped-axis orientation, what an edge in the plane
contributes, the overlap refusals, the fan walk that links boundary edges, the
outer/hole decision, ring orientation and canonical form, or the chain
normalisation. It must turn a test red: the hand-derived cases of
`tests/section_request.rs`, the region conformance suite, or the unit tests of
the region tracing.
"""
import probe_lib

R = "crates/algorithms/reference/src/section/regions.rs"
Q = "crates/algorithms/reference/src/section/request.rs"
C = "crates/algorithms/reference/src/section/chains.rs"
TESTS = [
    ["-p", "axiolid-reference", "--test", "section_request"],
    ["-p", "axiolid-reference", "--test", "section_conformance"],
    ["-p", "axiolid-reference", "--lib", "section"],
]

SEGMENT = (
    "insert_segment(\n"
    "                &mut segments,\n"
    "                Segment::new(NodeKey::Vertex(edge.0), NodeKey::Vertex(edge.1))?,\n"
    "                limits,\n"
    "            )?;"
)

MUTANTS = [
    (
        'in-plane triangle refused', Q,
        '        let zero_count = signs.iter().filter(|&&sign| sign == Sign::Zero).count();',
        '        let zero_count = signs.iter().filter(|&&sign| sign == Sign::Zero).count();\n'
        '        if zero_count == 3 {\n'
        '            return Err(GeomError::Degenerate("in-plane triangle refused".into()));\n'
        '        }',
        TESTS,
    ),
    (
        'dropped-axis sign ignored', R,
        'Ok(sign_times(projected, self.frame_sign))', 'Ok(projected)', TESTS,
    ),
    (
        'side mapping swapped', R,
        'if self.facing == Sign::Positive {', 'if self.facing != Sign::Positive {', TESTS,
    ),
    (
        'region edge emits a cut segment', Q,
        '[Incidence::In, Incidence::Off(_)] | [Incidence::Off(_), Incidence::In] => {}',
        '[Incidence::In, Incidence::Off(_)] | [Incidence::Off(_), Incidence::In] => {\n'
        '                ' + SEGMENT + '\n            }',
        TESTS,
    ),
    (
        'touching edge emits a cut segment', Q,
        '                    touching_edges += 1;',
        '                    touching_edges += 1;\n                    ' + SEGMENT,
        TESTS,
    ),
    (
        'touching edge not counted', Q,
        '                    touching_edges += 1;', '                    touching_edges += 0;', TESTS,
    ),
    (
        'overlap check skipped', R,
        'if interiors_overlap(&flat[i], &flat[j])? {',
        'if false && interiors_overlap(&flat[i], &flat[j])? {', TESTS,
    ),
    (
        'separating test strict', R,
        'if orient_sign(p, q, r)? == Sign::Positive {',
        'if orient_sign(p, q, r)? != Sign::Negative {', TESTS,
    ),
    (
        'overlap refusals swapped', R,
        'return Err(refuse(if first.facing == second.facing {',
        'return Err(refuse(if first.facing != second.facing {', TESTS,
    ),
    (
        'shared-edge pair counted', R,
        'if shared >= 2 && first.facing == second.facing {',
        'if false && shared >= 2 && first.facing == second.facing {', TESTS,
    ),
    (
        'pair budget unchecked', R,
        'if pairs > max_pairs {', 'if false && pairs > max_pairs {', TESTS,
    ),
    (
        'fan walk ignores component', R,
        'edge = leaving_edge(set, &boundary, edge.1, entering)?;',
        'let _ = entering;\n'
        '            edge = *boundary\n'
        '                .range((edge.1, 0)..)\n'
        '                .next()\n'
        '                .ok_or_else(internal_topology_error)?;',
        TESTS,
    ),
    (
        'outer and hole swapped', R,
        'Sign::Positive => Ok(true),', 'Sign::Positive => Ok(false),', TESTS,
    ),
    (
        'Above rings not reversed', R,
        'facing == Sign::Negative)?;', 'false)?;', TESTS,
    ),
    (
        'collinear region vertex kept', R,
        'if projection.orient(a, b, c)? != Sign::Zero {',
        'if true || projection.orient(a, b, c)? != Sign::Zero {',
        TESTS,
    ),
    (
        'cut-segment node simplified away', R,
        'if protected.contains(&ring[at]) {', 'if false && protected.contains(&ring[at]) {',
        TESTS,
    ),
    (
        'open chain direction unnormalised', C,
        'std::cmp::Ordering::Greater => points.reverse(),',
        'std::cmp::Ordering::Greater => {}', TESTS,
    ),
    (
        'saddle vertex accepted', C,
        '(_, 1) | (_, 2) => {}',
        '(_, 1) | (_, 2) | (NodeKey::Vertex(_), _) => {}', TESTS,
    ),
]

probe_lib.main(MUTANTS, timeout=1800, reject_uncompiled=True)
