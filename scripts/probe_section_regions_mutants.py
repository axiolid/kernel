"""Mutation probe for in-plane section regions (#359, ADR 0091).

Each mutant breaks one decision of `ScalarSection::section_request`: the side
of an in-plane face, the dropped-axis orientation, what an edge in the plane
contributes, the overlap refusals, the fan walk that links boundary edges, the
outer/hole decision, ring orientation and canonical form, or the chain
normalisation. The sheet mutants break what makes an open mesh
different: what a chain end lies on, the registry's topology check for a
sheet, the fold test and the orientation of a sheet's closed chain. A mutant
must turn a test red: the hand-derived cases of `tests/section_request.rs`,
the region conformance suite, or the unit tests of the region tracing.
"""
import probe_lib

R = "crates/algorithms/reference/src/section/regions.rs"
Q = "crates/algorithms/reference/src/section/request.rs"
D = "crates/execution/dispatch/src/section_request.rs"
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
        '                if region_edge_borders_cut(*facing, *third) {',
        '                ' + SEGMENT + '\n'
        '                if region_edge_borders_cut(*facing, *third) {',
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
        'std::cmp::Ordering::Greater => {\n'
        '                points.reverse();\n'
        '                ends.reverse();\n'
        '            }',
        'std::cmp::Ordering::Greater => {}', TESTS,
    ),
    (
        'open chain ends not swapped with the points', C,
        '                points.reverse();\n                ends.reverse();',
        '                points.reverse();', TESTS,
    ),
    (
        'saddle vertex accepted', C,
        '(_, 1) | (_, 2) => {}',
        '(_, 1) | (_, 2) | (NodeKey::Vertex(_), _) => {}', TESTS,
    ),
    (
        'protected chain vertex simplified', C,
        'if ctx.region_vertices.contains(&vertex) {',
        'if false && ctx.region_vertices.contains(&vertex) {', TESTS,
    ),
    (
        'open chain endpoints simplified', C,
        'if !ring && (index == 0 || index == last) {',
        'if false && !ring && (index == 0 || index == last) {', TESTS,
    ),
    (
        'edge node dropped without coplanarity', C,
        'NodeKey::Edge(edge) => Ok(ctx.coplanar_edges.contains(&edge)),',
        'NodeKey::Edge(_) => Ok(true),', TESTS,
    ),
    (
        'closed chain direction ignored', C,
        'if forward == Some(false) {', 'if false && forward == Some(false) {', TESTS,
    ),
    (
        'boundary edge end misclassified', C,
        '            NodeKey::Edge(edge) if self.boundary_edges.contains(&edge) => {\n'
        '                Ok(ChainEnd::MeshBoundary)',
        '            NodeKey::Edge(edge) if self.boundary_edges.contains(&edge) => {\n'
        '                Ok(ChainEnd::InPlaneRegion)',
        TESTS,
    ),
    (
        'boundary vertex end missing', C,
        'NodeKey::Vertex(vertex) if self.boundary_vertices.contains(&vertex) => {',
        'NodeKey::Vertex(vertex) if false && self.boundary_vertices.contains(&vertex) => {',
        TESTS,
    ),
    (
        'Sheet accepts inconsistent winding', D,
        'if health.non_manifold_edges > 0 || health.inconsistent_winding_edges > 0 {',
        'if health.non_manifold_edges > 0 {', TESTS,
    ),
    (
        'fold pair skipped', R,
        'if shared >= 2 && first.facing == second.facing {', 'if shared >= 2 {', TESTS,
    ),
    (
        'sheet boundary edge in the plane not counted as touching', Q,
        '[Incidence::Off(..)] if sheet => touching_edges += 1,',
        '[Incidence::Off(..)] if sheet => touching_edges += 0,',
        TESTS,
    ),
    (
        'sheet closed chain orientation ignored', C,
        '        Sign::Negative => points.reverse(),',
        '        Sign::Negative => {}', TESTS,
    ),
    (
        'sheet chain orientation by float area', C,
        '    match exact_area_sign(points)? {',
        '    match if signed_area(points) < 0.0 {\n'
        '        Sign::Negative\n'
        '    } else {\n'
        '        Sign::Positive\n'
        '    } {',
        TESTS,
    ),
    (
        'pinched extreme vertex has no fallback order', R,
        'for order in EXTREME_ORDERS {',
        'for order in EXTREME_ORDERS.into_iter().take(1) {', TESTS,
    ),
]

probe_lib.main(MUTANTS, timeout=1800, reject_uncompiled=True)
