"""Mutation probe for ruled surfaces (#402, ADR 0085 format 1.5).

Each mutant breaks one decision `SurfaceRelation::Ruled` rests on: in
`axiolid-model` the builder's validation (an empty or non-finite interval, a
non-finite apex, two apexes and one curve run twice admitted, the curve
sides widened beyond an atomic `Curve3`, the second side left out of the
references) and the wire format (the variant read in a 1.4 payload, a
payload never labelled 1.5, the feature test missing the variant); in
`axiolid-mesh-compile` the sheet (the curve term without the interval's
square, the twist term dropped, wrongly scaled or not divided by the rows,
the apex shortcut applied to every surface, the sides' tangent breaks not
forced as columns, a reversed interval read forward, the conic's turn and
the spline's domain checks dropped, sides that coincide admitted, the winding
reversed, the apex vertex not shared, the closure reported as a solid) and
the columns not split by the twist estimate, a triangle of no area skipped
silently, a collapsed triangle returned, a column whose ruling vanishes
given its rows' vertices, a cell whose end rulings both vanish left whole,
sides called coincident without the midpoint ruling or the gap bound, each
of the four cell caps dropped)
and the refusals that name what a surface relation is not (a B-rep face on
one, the exact compiler's family name, a curve-bounded basis, a sweep's
reference surface). Each must turn a test red.

Every mutant runs in its own process group with a timeout that kills the
whole group, and the file is restored whatever happens.
"""
import probe_lib

TIMEOUT = 1800
V = "crates/representations/modeling/graph/src/validation.rs"
S = "crates/representations/modeling/graph/src/surface_relation.rs"
G = "crates/representations/modeling/graph/src/wire.rs"
R = "crates/execution/compile/src/ruled.rs"
M = "crates/execution/compile/src/compiler.rs"
F = "crates/execution/compile/src/brep/view.rs"
E = "crates/execution/compile/src/exact.rs"
B = "crates/execution/compile/src/bounded.rs"

RULED = ["-p", "axiolid-model", "--all-features", "--test", "ruled_surface"]
WIRE = ["-p", "axiolid-model", "--all-features", "--test", "wire"]
SHEET = ["-p", "axiolid-mesh-compile", "--all-features", "--test", "ruled_surface"]

MUTANTS = [
    # The builder's validation.
    ("an empty interval is admitted", V,
     "            if interval.start == interval.end {",
     "            if false {", [RULED]),
    ("a non-finite interval is admitted", V,
     "            if !interval.start.is_finite() || !interval.end.is_finite() {",
     "            if false {", [RULED]),
    ("a non-finite apex is admitted", V,
     "            if point.is_finite() {",
     "            if true {", [RULED]),
    ("two apexes are admitted", V,
     "        (RuledSide::Apex(_), RuledSide::Apex(_)) => Err(GraphError::InvalidRuledSurface {",
     "        (RuledSide::Apex(_), RuledSide::Apex(_)) if false => Err(GraphError::InvalidRuledSurface {",
     [RULED]),
    ("one curve run twice over one interval is admitted", V,
     "        ) if ia == ib && (a == b || nodes[a.index()] == nodes[b.index()]) => {",
     "        ) if false => {", [RULED]),
    ("equal curve values under two nodes are told apart", V,
     "(a == b || nodes[a.index()] == nodes[b.index()])",
     "(a == b)", [RULED]),
    ("a side is widened beyond an atomic curve3", V,
     "            Self::AtomicCurve3 => return matches!(node, GeometryNode::Curve3(_)),",
     "            Self::AtomicCurve3 => return Self::Curve3.accepts(reference, node, nodes),",
     [RULED]),
    ("the second side is left out of the references", S,
     "                for side in [first, second] {",
     "                for side in [first] {", [RULED]),
    # The wire format (ADR 0085, format 1.5).
    ("the variant is read in a 1.4 payload", G,
     "    if version < RULED_SINCE && is_ruled(node) {",
     "    if false && is_ruled(node) {", [WIRE]),
    ("a payload is never labelled 1.5", G,
     "        version = version.max(RULED_SINCE);",
     "        version = version.max(ALONG_Z_SINCE);", [WIRE]),
    ("the feature test misses the variant", G,
     "        GeometryNode::SurfaceRelation(SurfaceRelation::Ruled { .. })",
     "        GeometryNode::SurfaceRelation(SurfaceRelation::Offset { .. })", [WIRE]),
    # The sheet.
    ("the curve term drops the interval's square", R,
     "                Ok((end - start) * (end - start) * bounds.second)",
     "                Ok(bounds.second)", [SHEET]),
    ("the twist term is dropped", R,
     "            0.25 * twist_at(&rims, i) / rows as Scalar",
     "            0.0", [SHEET]),
    ("the twist term is an eighth", R,
     "            0.25 * twist_at(&rims, i) / rows as Scalar",
     "            0.125 * twist_at(&rims, i) / rows as Scalar", [SHEET]),
    ("the twist term is not divided by the rows", R,
     "            0.25 * twist_at(&rims, i) / rows as Scalar",
     "            0.25 * twist_at(&rims, i)", [SHEET]),
    ("the curve term is a quarter", R,
     "            curve_term = curve_term.max(0.125 * width * width * second);",
     "            curve_term = curve_term.max(0.25 * width * width * second);", [SHEET]),
    ("the apex shortcut is applied to every surface", R,
     "    let apex = sides.iter().any(|side| matches!(side, Side::Apex(_)));",
     "    let apex = true;", [SHEET]),
    ("the sides' tangent breaks are not forced", R,
     "        for t in axiolid_reference::bound::continuity_breaks3(curve, 1) {",
     "        for t in Vec::<Scalar>::new() {", [SHEET]),
    ("a reversed interval is read forward", R,
     "            let (start, end) = (interval.start, interval.end);",
     "            let (start, end) = (\n"
     "                interval.start.min(interval.end),\n"
     "                interval.start.max(interval.end),\n"
     "            );", [SHEET]),
    ("a conic over more than a turn is let through", R,
     "                    if (end - start).abs() > TAU {",
     "                    if false {", [SHEET]),
    ("a spline interval outside its domain is let through", R,
     "                    if !(inside(start) && inside(end)) {",
     "                    if false {", [SHEET]),
    ("sides that coincide are admitted", R,
     "    if equal == 0 {",
     "    if true {", [SHEET]),
    ("the winding is reversed", R,
     "                indices.extend([p, q, r]);",
     "                indices.extend([p, r, q]);", [SHEET]),
    ("the apex vertex is not shared", R,
     "                Some(s) if matches!(sides[s], Side::Apex(_)) => {",
     "                Some(s) if false && matches!(sides[s], Side::Apex(_)) => {", [SHEET]),
    ("the closure is reported as a solid", M,
     "            }) => crate::ruled::ruled(graph, first, second, options).map(|(mesh, bound)| {\n"
     "                Built::with_closure(mesh, axiolid_mesh_compile_contract::MeshClosure::Surface)",
     "            }) => crate::ruled::ruled(graph, first, second, options).map(|(mesh, bound)| {\n"
     "                Built::with_closure(mesh, axiolid_mesh_compile_contract::MeshClosure::Solid)",
     [SHEET]),
    # What a surface relation is not.
    ("a B-rep face on a surface relation is a malformed graph again", F,
     "        Some(GeometryNode::SurfaceRelation(_)) => Err(GeomError::UnsupportedInput {",
     "        Some(GeometryNode::SurfaceRelation(_)) if false => Err(GeomError::UnsupportedInput {",
     [SHEET]),
    ("the exact compiler names a surface relation as a surface", E,
     '        GeometryNode::SurfaceRelation(_) => "surface relation",',
     '        GeometryNode::SurfaceRelation(_) => "surface",', [SHEET]),
    ("a curve-bounded basis that is a relation is a malformed graph again", B,
     "        Some(GeometryNode::SurfaceRelation(_)) => {",
     "        Some(GeometryNode::SurfaceRelation(_)) if false => {", [SHEET]),
    ("a sweep's reference that is a relation is a malformed graph again", M,
     "            Some(GeometryNode::SurfaceRelation(_)) => Err(GeomError::UnsupportedInput {",
     "            Some(GeometryNode::SurfaceRelation(_)) if false => Err(GeomError::UnsupportedInput {",
     [SHEET]),
    # The grid.
    ("the columns are not split by the twist estimate", R,
     "        if split > 1 {",
     "        if false {", [SHEET]),
    ("a triangle of no area is skipped", R,
     "                if p == q || q == r || r == p {",
     "                if p == q\n"
     "                    || q == r\n"
     "                    || r == p\n"
     "                    || positions[p as usize] == positions[q as usize]\n"
     "                    || positions[q as usize] == positions[r as usize]\n"
     "                    || positions[r as usize] == positions[p as usize]\n"
     "                {", [SHEET]),
    ("a cell whose end rulings both vanish is left whole", R,
     "        if !(vanish[i] && vanish[i + 1]) {",
     "        if true {", [SHEET]),
    ("sides are called coincident without the midpoint ruling", R,
     "        if k == 0.0 || (mid_vanishes && 0.125 * width * width * k <= 0.5 * chord) {",
     "        if k == 0.0 || 0.125 * width * width * k <= 0.5 * chord {", [SHEET]),
    ("sides are called coincident without the gap bound", R,
     "        if k == 0.0 || (mid_vanishes && 0.125 * width * width * k <= 0.5 * chord) {",
     "        if k == 0.0 || mid_vanishes {", [SHEET]),
    ("the halving of vanishing cells has no cap", R,
     "        if columns.len() > MAX_CELLS_PER_DIRECTION {",
     "        if false {", [SHEET]),
    ("a column whose ruling vanishes keeps its rows' vertices", R,
     "        if a == b {",
     "        if false {", [SHEET]),
    ("a collapsed triangle is returned", R,
     "                if (qp - pp).cross(rp - pp) == Vec3::ZERO {",
     "                if false {", [SHEET]),
    ("the column cap is dropped", R,
     "        if columns.len() + n > MAX_CELLS_PER_DIRECTION + 1 {",
     "        if false {", [SHEET]),
    ("the split cap is dropped", R,
     "    if (columns.len() - 1).saturating_mul(split) > MAX_CELLS_PER_DIRECTION {",
     "    if false {", [SHEET]),
    ("the total cap is dropped", R,
     "    if cells.saturating_mul(rows) > MAX_CELLS {",
     "    if false {", [SHEET]),
]

probe_lib.main(MUTANTS, timeout=TIMEOUT, reject_uncompiled=True)
