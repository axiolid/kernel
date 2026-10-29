"""Mutation probe for the cut regions of a section request (#372, ADR 0092).

Each mutant breaks one decision that turns the cut of a solid into closed
regions: which region edges bound the cut and in which direction, how the
rings are linked and simplified, the pinch refusal, the exact points the rings
are decided on (#377: the frame sign, the weight normalisation, the exact
fallback, the box of a point), the outer-or-hole turn at the smallest vertex,
the nesting (innermost container, containment, crossings, collinear overlaps),
the budgets, and the registry checks of the answer. A mutant must turn a test red:
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
P = "crates/algorithms/reference/src/section/plane_points.rs"
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
        'ring orientation by float area', N,
        '        ring.hole = is_hole(exact, &ring.nodes)?;',
        '        ring.hole = super::signed_area(&ring.points) < 0.0;', TESTS,
    ),
    (
        'ring role taken from the first vertex, not the smallest', N,
        '        if exact.lexicographic(nodes[index], nodes[smallest])? == Ordering::Less {',
        '        if false && exact.lexicographic(nodes[index], nodes[smallest])? == Ordering::Less {',
        TESTS,
    ),
    (
        'orientation ignores the frame sign', P,
        '        Ok(sign_times(sign, self.frame_sign))',
        '        Ok(sign)', TESTS,
    ),
    (
        'crossing weight not made positive', P,
        '        if point[2].sign() == Some(Sign::Negative) {',
        '        if false && point[2].sign() == Some(Sign::Negative) {', TESTS,
    ),
    (
        'orientation fallback guesses a sign', P,
        '            None => determinant(&a.exact, &b.exact, &c.exact)\n'
        '                .sign()\n'
        '                .ok_or_else(undecided)?,',
        '            None => Sign::Positive,', TESTS,
    ),
    (
        'comparison fallback guesses a sign', P,
        '            None => difference(&a.exact, &b.exact, axis)\n'
        '                .sign()\n'
        '                .ok_or_else(undecided)?,',
        '            None => Sign::Positive,', TESTS,
    ),
    (
        'vertex off the plane accepted', P,
        '        if self.plane_value(point)?.sign() != Some(Sign::Zero) {',
        '        if false && self.plane_value(point)?.sign() != Some(Sign::Zero) {', TESTS,
    ),
    (
        'edge on one side of the plane crossed', P,
        '        if !opposite {',
        '        if false && !opposite {', TESTS,
    ),
    (
        'box of a crossing not divided by its weight', P,
        '            enclosure[0].quotient(enclosure[2]),',
        '            enclosure[0],', TESTS,
    ),
    (
        'fast path orients with a float determinant', P,
        '            return Ok(sign_times(orient_sign(fa, fb, fc)?, self.frame_sign));',
        '            let det = (fb.x - fa.x) * (fc.y - fa.y) - (fb.y - fa.y) * (fc.x - fa.x);\n'
        '            let sign = if det > 0.0 {\n'
        '                Sign::Positive\n'
        '            } else if det < 0.0 {\n'
        '                Sign::Negative\n'
        '            } else {\n'
        '                Sign::Zero\n'
        '            };\n'
        '            return Ok(sign_times(sign, self.frame_sign));', TESTS,
    ),
    (
        'fast path orientation ignores the frame sign', P,
        '            return Ok(sign_times(orient_sign(fa, fb, fc)?, self.frame_sign));',
        '            return Ok(orient_sign(fa, fb, fc)?);', TESTS,
    ),
    (
        'fast path comparison reversed', P,
        '            return coordinate(fa, axis)\n'
        '                .partial_cmp(&coordinate(fb, axis))',
        '            return coordinate(fb, axis)\n'
        '                .partial_cmp(&coordinate(fa, axis))', TESTS,
    ),
    (
        'rounded crossing stored as a binary64 point', P,
        '    (exact_candidate.mul(weight) == *numerator).then_some(candidate)',
        '    Some(candidate)', TESTS,
    ),
    (
        'edge fixed in a coordinate taken from its other end', P,
        '            exact_quotient(&point[0], &point[2], fixed(a.x, b.x)),',
        '            exact_quotient(&point[0], &point[2], fixed(a.x, b.x).map(|x| x.next_up())),',
        TESTS,
    ),
    (
        'vertex node not checked against its certified side', P,
        '                if certified.sign != Sign::Zero {',
        '                if false && certified.sign != Sign::Zero {', TESTS,
    ),
    (
        'lexicographic order by the second coordinate only', P,
        '        match self.compare(a, b, 0)? {\n'
        '            Ordering::Equal => self.compare(a, b, 1),\n'
        '            order => Ok(order),\n'
        '        }',
        '        self.compare(a, b, 1)', TESTS,
    ),
    (
        'box pruning inverted', N,
        '        self.core[0] >= other.bounds[0]',
        '        self.core[0] <= other.bounds[0]', TESTS,
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
        'if exact.compare(b, point, 1)? == Ordering::Greater && side == Sign::Positive {',
        'if exact.compare(b, point, 1)? == Ordering::Greater && side == Sign::Negative {',
        TESTS,
    ),
    (
        'proper crossing accepted', N,
        'return Ok(around_first[0] != around_first[1] && around_second[0] != around_second[1]);',
        'return Ok(false);', TESTS,
    ),
    (
        'collinear overlap direction flipped', N,
        '        == (exact.compare(d, c, axis)? == Ordering::Greater);',
        '        != (exact.compare(d, c, axis)? == Ordering::Greater);', TESTS,
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
        '        refuse_spikes(exact, &ring.nodes)?;',
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
        'ring without area not refused as degenerate', K,
        '        if !turns {',
        '        if false && !turns {', TESTS,
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
