"""Mutation probe for a boolean's curved operand meshed by its section angle (#380).

Each mutant weakens one decision of `section_angle.rs` or of its wiring in
`compiler.rs`; the `section_angle` unit tests or `tests/boolean_deviation.rs`
must fail on it. A build error also counts as killed.

Not probed, because no fixture separates them:

- The second-order term of a face's ball radius (`0.5 * (duu hx^2 + ...)`):
  the balls of the fixtures meet or miss by far more than that term, so no
  factor a test can read changes. The ball is a Taylor bound
  (`section_angle::ball`), not a value a test pins.
- The `1e-9` slacks on a cosine and on a feasibility test (`SLACK`): they only
  keep a tangent contact read with rounding from falling out, which the
  depth test then drops anyway; a placed tangent hole sits on both sides of
  them.
- The frame each operand is read in by `boolean_in_shared_frame` and
  `clip_in_subject_frame` (`subject_budget`, `tool_budget`): read at the
  boolean's own budget they find no mesh at that key, return an error, and
  `build_boolean` falls back to the world-frame cut with the operands
  meshed as before, which is never worse by design (ADR 0096). The meshes
  and the bounds are the same to rounding.
- The visit limit of `inventory` and the pcurve quarter points of
  `parameter_box`: a limit and a sampling, not a decision.
"""
import probe_lib

S = "crates/execution/compile/src/section_angle.rs"
C = "crates/execution/compile/src/compiler.rs"
UNIT = ["-p", "axiolid-mesh-compile", "--lib", "section_angle"]
BUDGET = ["-p", "axiolid-mesh-compile", "--lib", "compiler::tests"]
INTEGRATION = ["-p", "axiolid-mesh-compile", "--test", "boolean_deviation"]
NESTED = ["-p", "axiolid-mesh-compile", "--test", "curved_wall_openings"]
FAR = ["-p", "axiolid-mesh-compile", "--test", "far_placed_beam"]
TESTS = [UNIT, BUDGET, INTEGRATION, NESTED, FAR]

MUTANTS = [
    ("the operands keep the boolean's chord", C,
     "            chord: self.chord * (factor / self.shrink),",
     "            chord: self.chord * 1.0,"),
    ("nested factors compound to their product", C,
     "        if factor >= self.shrink {\n"
     "            return self;\n"
     "        }\n"
     "        Self {\n"
     "            chord: self.chord * (factor / self.shrink),\n"
     "            shrink: factor,\n"
     "            ..self\n"
     "        }",
     "        Self {\n"
     "            chord: self.chord * factor,\n"
     "            shrink: self.shrink * factor,\n"
     "            ..self\n"
     "        }"),
    ("the factors of the two operands swapped", C,
     [("                        EvalKey::new(*left, budget.shrunk(factors.left)),\n"
       "                        EvalKey::new(*right, budget.shrunk(factors.right)),",
       "                        EvalKey::new(*left, budget.shrunk(factors.right)),\n"
       "                        EvalKey::new(*right, budget.shrunk(factors.left)),"),
      ("            _ => (factors.left, factors.right),",
       "            _ => (factors.right, factors.left),")],
     TESTS),
    ("a clip's subject read with the region's factor", C,
     [("                    let factor = if clip.subject == *left {",
       "                    let factor = if clip.subject != *left {"),
      ("            Some(clip) if clip.subject != left => (factors.right, factors.left),",
       "            Some(clip) if clip.subject == left => (factors.right, factors.left),")],
     TESTS),
    ("the floor dropped", S,
     "    ((s * GRID + NOISE).floor() / GRID).max(FLOOR)",
     "    ((s * GRID + NOISE).floor() / GRID)"),
    ("quantization rounds up", S,
     "    ((s * GRID + NOISE).floor() / GRID).max(FLOOR)",
     "    ((s * GRID + NOISE).ceil() / GRID).max(FLOOR)"),
    ("rounding noise not credited before the floor", S,
     "    ((s * GRID + NOISE).floor() / GRID).max(FLOOR)",
     "    ((s * GRID).floor() / GRID).max(FLOOR)"),
    ("a graze within the tolerance read as a section", S,
     "            if depth(s, curvature_radius(surface, *lo, *hi)) <= reach {",
     "            if depth(s, curvature_radius(surface, *lo, *hi)) <= 0.0 * reach {"),
    ("the torus's minor angle ignored", S,
     "    let (whole_u, whole_v) = (u.1 - u.0 >= TAU - 1e-9, v.1 - v.0 >= TAU - 1e-9);",
     "    let (whole_u, whole_v) = (u.1 - u.0 >= TAU - 1e-9, true);"),
    ("the ends of the torus's minor angle skipped", S,
     "        if !whole_v {",
     "        if false && !whole_v {"),
    ("the torus feasibility roots skipped", S,
     "    if discriminant >= 0.0 {",
     "    if false && discriminant >= 0.0 {"),
    ("a torus face's angle range ignored", S,
     "    let (whole_u, whole_v) = (u.1 - u.0 >= TAU - 1e-9, v.1 - v.0 >= TAU - 1e-9);",
     "    let (whole_u, whole_v) = (true, v.1 - v.0 >= TAU - 1e-9);"),
    ("a torus section read past its feasibility", S,
     "        if q.abs() > reach * (1.0 + SLACK) {",
     "        if q.abs() > reach * 1e9 {"),
    ("the torus's axis term dropped from the roots", S,
     "    let c = rel.delta * rel.delta - small * small * rel.alpha * rel.alpha;",
     "    let c = rel.delta * rel.delta;"),
    ("the cylinder ignores the axial range", S,
     "    Some(sine(rel.m * lo.abs().max(hi.abs())))",
     "    Some(sine(rel.m))"),
    ("the cone's axis term dropped (square to the axis)", S,
     "        return (v0 <= height && height <= v1).then(|| sine(sin * rel.alpha));",
     "        return (v0 <= height && height <= v1).then(|| sine(0.0));"),
    ("the cone's axis term dropped (oblique)", S,
     "    let cosine = |c: Scalar| (cos * rel.m * c - sin * rel.alpha).abs();",
     "    let cosine = |c: Scalar| (cos * rel.m * c).abs();"),
    ("the cone read past its apex", S,
     "        v0 = v0.max(-r0 / slope);",
     "        v0 = v0.max(Scalar::NEG_INFINITY);"),
    ("the sphere's sine without the square root", S,
     "    (cos.abs() <= 1.0 + SLACK).then(|| sine(cos))",
     "    (cos.abs() <= 1.0 + SLACK).then(|| 1.0 - cos.abs())"),
    ("balls never meet", S,
     "            (Some(a), Some(b)) => a.centre.distance(b.centre) <= a.radius + b.radius + reach,",
     "            (Some(a), Some(b)) => a.centre.distance(b.centre) < 0.0 * (a.radius + b.radius + reach),"),
    ("balls always meet", S,
     "            (Some(a), Some(b)) => a.centre.distance(b.centre) <= a.radius + b.radius + reach,",
     "            (Some(a), Some(b)) => a.centre.distance(b.centre) <= f64::INFINITY * (a.radius + b.radius + reach),"),
    ("a half-space leaf contributes no plane", S,
     "                GeometryNode::HalfSpace(half_space) => {\n"
     "                    faces.extend(half_space_plane(half_space, &transform));",
     "                GeometryNode::HalfSpace(half_space) => {\n"
     "                    let _ = half_space;"),
    ("a half-space's normal moved by the matrix, not its inverse transpose", S,
     "    let normal = (matrix.inverse().transpose() * half_space.boundary.normal).normalize_or_zero();",
     "    let normal = (matrix * half_space.boundary.normal).normalize_or_zero();"),
    ("nested placements composed the wrong way round", S,
     "                    transform *= instance.transform;",
     "                    transform = instance.transform * transform;"),
    ("the planar pre-filter never taken", S,
     "    if certainly_planar(graph, *left) && certainly_planar(graph, *right) {",
     "    if false && certainly_planar(graph, *left) && certainly_planar(graph, *right) {"),
    ("a rounded rectangle read as planar", S,
     "        Profile::Rectangle(r) => r.outer_radius.is_none() && r.inner_radius.is_none(),",
     "        Profile::Rectangle(r) => r.inner_radius.is_none(),"),
    ("an arc read as straight", S,
     "                    .all(|s| matches!(s.curve, Curve2::Line(_) | Curve2::Polyline(_)))",
     "                    .all(|s| matches!(s.curve, Curve2::Line(_) | Curve2::Polyline(_) | Curve2::Circle(_)))"),
]

probe_lib.main(MUTANTS, tests=TESTS, timeout=1800, reject_uncompiled=True)
