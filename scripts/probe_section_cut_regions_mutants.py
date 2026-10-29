"""Mutation probe for the cut regions of a section request (#372, ADR 0092).

Each mutant breaks one decision that turns the cut of a solid into closed
regions: which region edges bound the cut and in which direction, how the
rings are linked and simplified, the pinch refusal, the exact area sign, the
nesting (innermost container, containment, crossings, collinear overlaps), the
budgets, and the registry checks of the answer. A mutant must turn a test red:
the hand-derived cases of `tests/section_request.rs`, the region conformance
suite, the unit tests of the ring tracing and nesting, or the registry tests
of `mesh_section_contract.rs`.

Not mutated, because no input reaches the difference: the mixed-containment
and tied-container refusals (a ring with vertices strictly inside and outside
another needs a crossing or contact, which the sweep refuses first), the
same-ring clause of the contact rule (the local test also sees both passages
of a repeated vertex), and the order of the two local tests (a crossing is seen
from both rings).
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
        'registry does not fall back past a provider that ignores the request',
        "crates/execution/dispatch/src/section.rs",
        'Ok(section) if request.cut_regions.is_some() && section.cut_regions.is_none() => {',
        'Ok(section) if false && request.cut_regions.is_some() && section.cut_regions.is_none() => {',
        TESTS,
    ),
    (
        'registry cut region budget unchecked', D,
        '            if cut_regions.len() > limits.max_regions {',
        '            if false && cut_regions.len() > limits.max_regions {', TESTS,
    ),
    (
        'registry cut pair budget unchecked', D,
        '            if section.evidence.cut_pairs > limits.max_pairs {',
        '            if false && section.evidence.cut_pairs > limits.max_pairs {', TESTS,
    ),
    (
        'conformance accepts rings that do not cover the chains',
        "crates/contracts/operations/mesh-section/src/region_conformance.rs",
        'fn cut_rings_cover_the_chains(section: &PlaneSection, cut: &[CutRegion]) -> bool {',
        'fn cut_rings_cover_the_chains(section: &PlaneSection, cut: &[CutRegion]) -> bool {\n'
        '    if true {\n        return true;\n    }', TESTS,
    ),
    (
        'ring edges of one ring never compared', N,
        '                if gap == 1 || gap == length - 1 {',
        '                if true {', TESTS,
    ),
    (
        'contact between rings never crosses', N,
        '    Ok(inside && outside)',
        '    Ok(false)', TESTS,
    ),
    (
        'spike accepted', N,
        '        refuse_spikes(&ring.points)?;',
        '        let _ = ring;', TESTS,
    ),
    (
        'hole inside a hole accepted', N,
        '(true, Some(container)) if !rings[container].hole => {',
        '(true, Some(container)) => {', TESTS,
    ),
    (
        'unequal degrees not an internal error', K,
        '        if out != inn {',
        '        if false && out != inn {', TESTS,
    ),
    (
        'ring without area taken for a hole', K,
        '            _ => {\n'
        '                return Err(GeomError::Degenerate(\n'
        '                    "cut region ring has no signed area".into(),\n'
        '                ))\n'
        '            }',
        '            _ => holes.push(canonical_rotation(points)),', TESTS,
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
