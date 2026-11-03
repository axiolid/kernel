"""Mutation probe for the merge of in-plane faces (#366, ADR 0093).

Each mutant breaks one decision of the exact arrangement of the in-plane
faces or of what is built on it: the class of a side (`Both`), the interior
test of a piece in a triangle, the side a carrier covers, where the source
edges are cut (proper crossings, an end inside the other segment, equal
points joined), the half-edge walk and the fan order that link boundary
pieces into rings, which collinear vertices stay, the owner of a hole, and
how `Both` regions join the cut (their rings, the pinch refusal, the
cut-side check, the registry check and the carried evidence). The budgets
(coverage pairs, the memory charge) and the output position of a crossing are
mutated too. A mutant must turn a test red: the hand-derived cases of
`tests/section_request.rs`, the region conformance suite, the unit tests of
the arrangement and of the regions, or the registry tests of
`mesh_section_contract.rs`.

Not mutated, because no input reaches the difference: the conformance check
that every `Both` region is a cut region (the fixtures pin the exact cut
regions, which already contain them).
"""
import probe_lib

A = "crates/algorithms/reference/src/section/arrangement.rs"
X = "crates/algorithms/planar/exact-arrangement/src/arrangement.rs"
T = "crates/algorithms/planar/exact-arrangement/src/trace.rs"
R = "crates/algorithms/reference/src/section/regions.rs"
Q = "crates/algorithms/reference/src/section/request.rs"
D = "crates/execution/dispatch/src/section_request.rs"
TESTS = [
    ["-p", "axiolid-reference", "--test", "section_request"],
    ["-p", "axiolid-reference", "--test", "section_conformance"],
    ["-p", "axiolid-reference", "--lib", "section"],
    ["-p", "axiolid-dispatch", "--all-features", "--test", "mesh_section_contract"],
    ["-p", "axiolid-exact-arrangement"],
]

MUTANTS = [
    (
        'Both mapped to Below', A,
        '(true, true) => Self::Both,', '(true, true) => Self::Below,', TESTS,
    ),
    (
        'interior test strict', A,
        [
            ('        if start == Sign::Negative {', '        if start != Sign::Positive {'),
            ('        if end == Sign::Negative {', '        if end != Sign::Positive {'),
        ],
        TESTS,
    ),
    (
        'carrier side flipped', A,
        [
            ('Sign::Positive => cover(piece, true, false, triangles[triangle].facing),',
             'Sign::Positive => cover(piece, false, true, triangles[triangle].facing),'),
            ('Sign::Negative => cover(piece, false, true, triangles[triangle].facing),',
             'Sign::Negative => cover(piece, true, false, triangles[triangle].facing),'),
        ],
        TESTS,
    ),
    (
        'every pair of segments taken for one of a vertex', X,
        '.position(|&other| other == end)',
        '.position(|_| true)', TESTS,
    ),
    (
        'proper crossings skipped', X,
        '            if proper {', '            if false && proper {', TESTS,
    ),
    (
        'ends of the later segment never split the earlier', X,
        [
            ('if around[0] == Sign::Zero && points.strictly_between(a, c, b)? {',
             'if false && around[0] == Sign::Zero && points.strictly_between(a, c, b)? {'),
            ('if around[1] == Sign::Zero && points.strictly_between(a, d, b)? {',
             'if false && around[1] == Sign::Zero && points.strictly_between(a, d, b)? {'),
        ],
        TESTS,
    ),
    (
        'ends of the earlier segment never split the later', X,
        [
            ('if around[2] == Sign::Zero && points.strictly_between(c, a, d)? {',
             'if false && around[2] == Sign::Zero && points.strictly_between(c, a, d)? {'),
            ('if around[3] == Sign::Zero && points.strictly_between(c, b, d)? {',
             'if false && around[3] == Sign::Zero && points.strictly_between(c, b, d)? {'),
        ],
        TESTS,
    ),
    (
        'coincident source vertices not joined', X,
        '        if points.same(local.exact[pair[0]], local.exact[pair[1]])? {\n'
        '            local.union(pair[0], pair[1]);',
        '        if points.same(local.exact[pair[0]], local.exact[pair[1]])? {\n'
        '            let _ = pair;',
        TESTS,
    ),
    (
        'equal points along a segment not joined', X,
        '            if points.same(local.exact[pair[0]], local.exact[pair[1]])? {\n'
        '                local.union(pair[0], pair[1]);',
        '            if points.same(local.exact[pair[0]], local.exact[pair[1]])? {\n'
        '                let _ = pair;',
        TESTS,
    ),
    (
        'pieces starting inside a triangle never tested', A,
        '    for (&triangle, held) in triangle_order.iter().zip(&triangle_swept) {',
        '    for (&triangle, held) in triangle_order.iter().zip(&triangle_swept).take(0) {',
        TESTS,
    ),
    (
        'triangles starting inside a piece never tested', A,
        '    for (&piece, held) in piece_order.iter().zip(&piece_swept) {',
        '    for (&piece, held) in piece_order.iter().zip(&piece_swept).take(0) {',
        TESTS,
    ),
    (
        'coverage pairs not counted', A,
        '        budget.count()?;\n        let (from, to) = (',
        '        let (from, to) = (', TESTS,
    ),
    (
        'memory charge disabled', X,
        'Some(needed) if needed <= budget => {', 'Some(_) => {', TESTS,
    ),
    (
        'crossing output position ignores the parameter', A,
        'let t = exact.segment_parameter(ends)?;', 'let t = 0.5;', TESTS,
    ),
    (
        'counter-clockwise next instead of clockwise', T,
        'let to = to_of(&fan[(index + fan.len() - step) % fan.len()]);',
        'let to = to_of(&fan[(index + step) % fan.len()]);', TESTS,
    ),
    (
        'fan order reversed within a half', X,
        [
            ('Sign::Positive => Ok(std::cmp::Ordering::Less),',
             'Sign::Positive => Ok(std::cmp::Ordering::Greater),'),
            ('Sign::Negative => Ok(std::cmp::Ordering::Greater),',
             'Sign::Negative => Ok(std::cmp::Ordering::Less),'),
        ],
        TESTS,
    ),
    (
        'collinear vertex where other boundary meets dropped', R,
        '            || on_boundary[vertex.0] != 2\n    };',
        '            || false\n    };', TESTS,
    ),
    (
        'hole owned by the outermost container', T,
        'if !within(points, &rings[candidate], &rings[other])? {',
        'if !within(points, &rings[other], &rings[candidate])? {', TESTS,
    ),
    (
        'pinch of a Both ring never noticed', R,
        'both_pinched |= distinct.len() != ring.nodes.len();',
        'both_pinched |= distinct.is_empty();', TESTS,
    ),
    (
        'Both rings reversed for the cut', R,
        '                nodes: outer.nodes.clone(),',
        '                nodes: outer.nodes.iter().rev().copied().collect(),', TESTS,
    ),
    (
        'Both rings not passed to the cut', Q,
        '                both_rings,\n                cut_limits,',
        '                Vec::new(),\n                cut_limits,', TESTS,
    ),
    (
        'pinched Both ring not refused', Q,
        '            if both_pinched {', '            if false && both_pinched {', TESTS,
    ),
    (
        'covered cut side not refused', Q,
        '                    .any(|&class| class != Class::Empty)',
        '                    .any(|_| false)', TESTS,
    ),
    (
        'hole containment pairs not counted', T,
        '            budget.count()?;\n'
        '            if within(points, &rings[hole], &rings[candidate])? {',
        '            if within(points, &rings[hole], &rings[candidate])? {', TESTS,
    ),
    (
        'innermost container pairs not counted', T,
        '                budget.count()?;\n'
        '                if !within(points, &rings[candidate], &rings[other])? {',
        '                if !within(points, &rings[candidate], &rings[other])? {', TESTS,
    ),
    (
        'collapsed ring accepted', R,
        '            || (0..local.len()).any(|index| local[index] == local[(index + 1) % local.len()])',
        '            || false', TESTS,
    ),
    (
        'piece pruned by the outer box of its crossings', A,
        'let (core, held) = (&piece_cores[piece], &triangle_boxes[triangle]);',
        'let (core, held) = (&piece_boxes[piece], &triangle_boxes[triangle]);', TESTS,
    ),
    (
        'registry accepts a Both region that is not a cut region', D,
        '                if !is_cut {', '                if false && !is_cut {', TESTS,
    ),
    (
        'registry drops the in-plane crossings from the evidence', D,
        '        .with_in_plane_crossings(evidence.in_plane_crossings)\n', '', TESTS,
    ),
]

probe_lib.main(MUTANTS, timeout=1800, reject_uncompiled=True)
