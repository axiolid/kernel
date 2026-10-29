"""Mutation probe for the cut regions of a section request (#372, ADR 0092).

Each mutant breaks one decision that turns the cut of a solid into closed
regions: which region edges bound the cut and in which direction, how the
rings are linked and simplified, the pinch refusal, the exact area sign, the
nesting (innermost container, containment, crossings, collinear overlaps), the
budgets, and the registry checks of the answer. A mutant must turn a test red:
the hand-derived cases of `tests/section_request.rs`, the region conformance
suite, the unit tests of the ring tracing and nesting, or the registry tests
of `mesh_section_contract.rs`.
"""
import probe_lib

Q = "crates/algorithms/reference/src/section/request.rs"
K = "crates/algorithms/reference/src/section/cut_regions.rs"
N = "crates/algorithms/reference/src/section/ring_nesting.rs"
D = "crates/execution/dispatch/src/section_request.rs"
TESTS = [
    ["-p", "axiolid-reference", "--test", "section_request"],
    ["-p", "axiolid-reference", "--test", "section_conformance"],
    ["-p", "axiolid-reference", "--lib", "section"],
    ["-p", "axiolid-dispatch", "--all-features", "--test", "mesh_section_contract"],
]

MUTANTS = [
    (
        'side test: any off-plane neighbour bounds the cut', Q,
        '    opposite(facing, third)',
        '    facing == third', TESTS,
    ),
    (
        'bounding region edge reversed', Q,
        'cut_edges.insert((*from, *to));',
        'cut_edges.insert((*to, *from));', TESTS,
    ),
    (
        'sheet request not refused by the provider', Q,
        '    if sheet && request.cut_regions.is_some() {',
        '    if false && sheet && request.cut_regions.is_some() {', TESTS,
    ),
    (
        'pinch not refused', K,
        '        if out > 1 {',
        '        if false && out > 1 {', TESTS,
    ),
    (
        'region run vertex always dropped', K,
        'NodeKey::Vertex(vertex) => !inputs.region_ring_vertices.contains(&vertex),',
        'NodeKey::Vertex(_) => true,', TESTS,
    ),
    (
        'region run vertex never dropped', K,
        'NodeKey::Vertex(vertex) => !inputs.region_ring_vertices.contains(&vertex),',
        'NodeKey::Vertex(_) => false,', TESTS,
    ),
    (
        'junction of a chain and a region edge dropped', K,
        '            _ => false,\n        };\n        if !droppable {',
        '            _ => true,\n        };\n        if !droppable {', TESTS,
    ),
    (
        'collinear chain nodes kept in the ring', K,
        '(Link::Cut, Link::Cut) => is_droppable(inputs.chains, previous, node, next)?,',
        '(Link::Cut, Link::Cut) => false,', TESTS,
    ),
    (
        'ring orientation by float area', K,
        '        match exact_area_sign(&points)? {',
        '        match if super::signed_area(&points) < 0.0 {\n'
        '            Sign::Negative\n'
        '        } else {\n'
        '            Sign::Positive\n'
        '        } {', TESTS,
    ),
    (
        'innermost container taken as the outermost', N,
        'let deepest = containers[index].iter().map(|&c| depth(c)).max();',
        'let deepest = containers[index].iter().map(|&c| depth(c)).min();', TESTS,
    ),
    (
        'outer ring inside an outer ring accepted', N,
        '(false, Some(container)) if rings[container].hole => {}',
        '(false, Some(_)) => {}', TESTS,
    ),
    (
        'mixed containment accepted', N,
        '                (false, true) => {}\n                _ => return Err(refused()),',
        '                (false, true) => {}\n                _ => {}', TESTS,
    ),
    (
        'winding walk counts upward edges against the point', N,
        'if b.y > point.y && side == Sign::Positive {',
        'if b.y > point.y && side == Sign::Negative {', TESTS,
    ),
    (
        'proper crossing accepted', N,
        'return Ok(around_first[0] != around_first[1] && around_second[0] != around_second[1]);',
        'return Ok(false);', TESTS,
    ),
    (
        'collinear overlap direction flipped', N,
        'let same_way = (along(b) > along(a)) == (along(d) > along(c));',
        'let same_way = (along(b) > along(a)) != (along(d) > along(c));', TESTS,
    ),
    (
        'pair budget unchecked', N,
        '    if *pairs > limits.max_pairs {',
        '    if false && *pairs > limits.max_pairs {', TESTS,
    ),
    (
        'vertex budget unchecked', N,
        '    if vertices > limits.max_region_vertices {',
        '    if false && vertices > limits.max_region_vertices {', TESTS,
    ),
    (
        'region budget unchecked', N,
        '    if outer_count > limits.max_regions {',
        '    if false && outer_count > limits.max_regions {', TESTS,
    ),
    (
        'registry accepts an answer without the requested cut regions', D,
        '        (Some(_), None) => Err(violation(\n'
        '            backend,\n'
        '            "provider did not trace the requested cut regions",\n'
        '        )),',
        '        (Some(_), None) => Ok(()),', TESTS,
    ),
    (
        'registry accepts cut regions that were not requested', D,
        '        (None, Some(_)) => Err(violation(\n'
        '            backend,\n'
        '            "provider returned cut regions that were not requested",\n'
        '        )),',
        '        (None, Some(_)) => Ok(()),', TESTS,
    ),
    (
        'registry accepts cut regions beside a Both region', D,
        '.any(|region| region.side == InPlaneSide::Both)',
        '.any(|_| false)', TESTS,
    ),
    (
        'registry lets a sheet request reach the provider', D,
        '    if request.cut_regions.is_some() && request.role == SectionRole::Sheet {',
        '    if false && request.cut_regions.is_some() && request.role == SectionRole::Sheet {',
        TESTS,
    ),
]

probe_lib.main(MUTANTS, timeout=1800, reject_uncompiled=True)
