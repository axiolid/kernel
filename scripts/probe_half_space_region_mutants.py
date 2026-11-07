"""Mutation probe for half-space regions as clip tools, and empty operands
(#326).

Each mutant breaks one decision of a region clip -- how the stand-in tool is
built leaf by leaf, what each leaf decides alone, the swap of an
intersection's operands, the dispatch -- or of the empty operands: the
identity table both compilers read, the empty exact body the booleans
absorb, the mesh compiler's identity cut and the pending skin it passes
through, the deviation of an empty result. Each must turn a test red.

Equivalent mutant, deliberately not listed: sizing every mesh slab by the
subject instead of by the stand-in so far. The results are the same volumes;
at a site of (6e5, 5.6e6, 200) it cuts 40 triangles instead of 28 and
certifies 9.1e-5 instead of 3.2e-5, both within the budget. (E2 is the
exact compiler's version: there the far faces coincide, and the report
stops reading exact.)
"""
import pathlib

import probe_lib

ROOT = pathlib.Path(__file__).resolve().parents[1]
CLIP = "crates/execution/compile/src/exact/clip.rs"
BOOLEAN = "crates/execution/compile/src/exact/boolean.rs"
EXACT = "crates/execution/compile/src/exact.rs"
COMPILER = "crates/execution/compile/src/compiler.rs"
EMPTY = "crates/execution/compile/src/empty.rs"
DEVIATION = "crates/execution/compile/src/deviation/boolean.rs"
REGIONS = "crates/execution/compile/src/regions.rs"

TESTS = [
    ["-p", "axiolid-mesh-compile", "--test", "half_space_region"],
    ["-p", "axiolid-mesh-compile", "--test", "empty_operand"],
    ["-p", "axiolid-mesh-compile", "--lib"],
]
TURNED = [["-p", "axiolid-mesh-compile", "--test", "turned_clip"]]
DEEP = [["-p", "axiolid-mesh-compile", "--test", "compiler", "a_deep_region"]]

REGIONS_LIB = [["-p", "axiolid-mesh-compile", "--lib", "regions"]]

# The memo's walk of a nest of intersections written recursively: a nest
# 100000 deep aborts the process on the test's stack, and no node is
# remembered, so a chain is quadratic.
WALK = (
    "                })) => {\n"
    "                    fn walk(\n"
    "                        graph: &GeometryGraph,\n"
    "                        id: NodeId,\n"
    "                        leaf: &dyn Fn(NodeId) -> bool,\n"
    "                    ) -> bool {\n"
    "                        match graph.get(id) {\n"
    "                            Some(GeometryNode::SolidOperation(SolidOperation::Boolean {\n"
    "                                left,\n"
    "                                right,\n"
    "                                operator: BooleanOperator::Intersection,\n"
    "                            })) => walk(graph, *left, leaf) && walk(graph, *right, leaf),\n"
    "                            Some(GeometryNode::Instance(i)) => walk(graph, i.source, leaf),\n"
    "                            Some(_) => leaf(id),\n"
    "                            None => false,\n"
    "                        }\n"
    "                    }\n"
    "                    let _ = (left, right);\n"
    "                    Ok(walk(graph, top, &leaf))\n"
    "                }\n"
)

MUTANTS = [
    # The exact stand-in.
    ("E1 only the first crossing leaf is used", CLIP,
     "            match self.leaf_cover(leaf, around)? {",
     "            match match self.leaf_cover(leaf, around)? {\n"
     "                Cover::Cuts(..) if tool.is_some() => Cover::Holds,\n"
     "                other => other,\n"
     "            } {"),
    ("E2 a further leaf is sized over the envelope", CLIP,
     "                Some((stand_in, _)) => edge_box(stand_in)?,",
     "                Some((_, _)) => envelope,"),
    ("E3 a leaf that misses decides nothing", CLIP,
     "                Cover::Misses => return Ok(Cover::Misses),\n"
     "                Cover::Holds => {}",
     "                Cover::Misses => {}\n"
     "                Cover::Holds => {}"),
    ("E4 a leaf that holds the box misses it", CLIP,
     "                Cover::Misses => return Ok(Cover::Misses),\n"
     "                Cover::Holds => {}",
     "                Cover::Misses => return Ok(Cover::Misses),\n"
     "                Cover::Holds => return Ok(Cover::Misses),"),
    ("E5 stand-ins with nothing in common hold", CLIP,
     "                                Err(BooleanError::EmptyResult) => {\n"
     "                                    match at_zero_tolerance(\n"
     "                                        &stand_in,\n"
     "                                        &side,\n"
     "                                        BooleanOperator::Intersection,\n"
     "                                    ) {\n"
     "                                        Ok(done) => done,\n"
     "                                        Err(BooleanError::EmptyResult) => return Ok(Cover::Misses),",
     "                                Err(BooleanError::EmptyResult) => {\n"
     "                                    match at_zero_tolerance(\n"
     "                                        &stand_in,\n"
     "                                        &side,\n"
     "                                        BooleanOperator::Intersection,\n"
     "                                    ) {\n"
     "                                        Ok(done) => done,\n"
     "                                        Err(BooleanError::EmptyResult) => (stand_in, BooleanReport::default()),"),
    ("E6 the exact swap is not a clip", CLIP,
     "            (true, false) if operator == BooleanOperator::Intersection => (right, left),",
     "            (true, false) if false && operator == BooleanOperator::Intersection => (right, left),"),
    ("E7 the exact clip is not dispatched", BOOLEAN,
     "        if let Some((subject, region)) = self.clip_of(left, right, operator)? {",
     "        if let Some((subject, region)) = self.clip_of(left, right, operator)?.filter(|_| false) {"),
    ("E10 a clip that leaves nothing is a refusal", CLIP,
     "return if emptiness_is_exact(&subject, &solid, operator) {\n"
     "                    Ok(Body::Empty(CLIP_EMPTIED))",
     "return if false {\n"
     "                    Ok(Body::Empty(CLIP_EMPTIED))"),
    ("E11 a boolean that leaves nothing is a refusal", BOOLEAN,
     "return if emptiness_is_exact(&subject, &tool, operator) {\n"
     "                        Ok(Body::Empty(emptied(operator)))",
     "return if false {\n"
     "                        Ok(Body::Empty(emptied(operator)))"),
    ("E12 an instance of nothing is refused", EXACT,
     "                // Nothing placed is still nothing.\n"
     "                let Body::Solid(source) = source else {\n"
     "                    return Ok(source);",
     "                // Nothing placed is still nothing.\n"
     "                let Body::Solid(source) = source else {\n"
     "                    return Err(unsupported(\"exact instance of nothing\"));"),
    # The identity table, both compilers.
    ("E8 a difference by nothing is nothing", EMPTY,
     "        BooleanOperator::Difference => Some(if subject_empty {\n"
     "            Identity::Empty\n"
     "        } else {\n"
     "            Identity::Subject\n"
     "        }),",
     "        BooleanOperator::Difference => Some(Identity::Empty),"),
    ("E9 a union with nothing is nothing", EMPTY,
     "            Some(if subject_empty && tool_empty {\n"
     "                Identity::Empty\n",
     "            Some(if true {\n"
     "                Identity::Empty\n"),
    ("M6 an empty subject is not read", COMPILER,
     "            subject.mesh.indices.is_empty(),\n"
     "            built_tool.is_some_and(",
     "            false,\n"
     "            built_tool.is_some_and("),
    # The deviation of an empty result.
    ("D1 an empty result against an empty mesh is unbounded", DEVIATION,
     "    with_readings(EMPTY_RESULT, DeviationBound::Certified(0.0), report)",
     "    with_readings(EMPTY_RESULT, DeviationBound::Unbounded(\"mutant\"), report)"),
    ("D2 an empty mesh is never read", DEVIATION,
     "    if mesh.indices.is_empty() {\n"
     "        return with_readings(",
     "    if true {\n"
     "        return with_readings("),
    # The mesh stand-in.
    ("M1 the mesh tool takes the first leaf only", COMPILER,
     "                self.region_tool_mesh(graph, &clip.leaves, &subject.mesh, options, budget.placed)?",
     "                self.region_tool_mesh(graph, &clip.leaves[..1], &subject.mesh, options, budget.placed)?"),
    ("M2 an empty mesh stand-in is a refusal", COMPILER,
     "            (Some(_), None, None) => crate::empty::identity(operator, false, true),",
     "            (Some(_), None, None) => {\n"
     "                return Err(GeomError::InvalidInput(\"mutant\".into()))\n"
     "            }"),
    ("M3 a region is scheduled as a dependency", COMPILER,
     "                    vec![same(clip.subject)]",
     "                    vec![same(clip.subject), same(clip.tool)]"),
    ("M4 a region of several leaves is left out of the subject's frame", COMPILER,
     "        // Each leaf of a region (#326) is placed in the subject's frame; one\n",
     "        if clip.leaves.len() > 1 {\n"
     "            return Ok(None);\n"
     "        }\n"
     "        // Each leaf of a region (#326) is placed in the subject's frame; one\n",
     TURNED + TESTS[:1]),
    ("M5 a pending skin is dropped at an identity", COMPILER,
     "            let skin = built.pending_skin.take();",
     "            let skin = {\n"
     "                built.pending_skin.take();\n"
     "                None\n"
     "            };",
     TURNED),
    ("F1 the memo reads a nest of intersections recursively", REGIONS,
     "                })) => match (known.get(left), known.get(right)) {\n"
     "                    (Some(false), _) | (_, Some(false)) => Ok(false),\n"
     "                    (Some(true), Some(true)) => Ok(true),\n"
     "                    (None, _) => Err(*left),\n"
     "                    (_, None) => Err(*right),\n"
     "                },\n",
     WALK, DEEP + REGIONS_LIB),
    ("F5e the exact region is read recursively", CLIP,
     "                        work.push((*right, placement));\n"
     "                        work.push((*left, placement));\n"
     "                        break;",
     "                        leaves.extend(self.region(*left)?.leaves);\n"
     "                        leaves.extend(self.region(*right)?.leaves);\n"
     "                        break;", DEEP),
    ("F5m the mesh region is read recursively", COMPILER,
     "                    work.push((*right, outer));\n"
     "                    work.push((*left, outer));\n"
     "                    break;",
     "                    leaves.extend(region(graph, *left)?);\n"
     "                    leaves.extend(region(graph, *right)?);\n"
     "                    break;", DEEP),
    ("F2a a boolean emptied within the tolerance is an identity", BOOLEAN,
     "return if emptiness_is_exact(&subject, &tool, operator) {",
     "return if true {"),
    ("F2b stand-ins meeting at the margin within the tolerance are refused", CLIP,
     "                                        Ok(done) => done,\n"
     "                                        Err(BooleanError::EmptyResult) => return Ok(Cover::Misses),",
     "                                        Ok(_) => return Err(unsupported(STAND_IN_WITHIN_TOLERANCE)),\n"
     "                                        Err(BooleanError::EmptyResult) => return Ok(Cover::Misses),"),
    ("F2d the prism along z meeting the plane's side at the margin is refused", CLIP,
     "                        Ok(done) => done,\n"
     "                        Err(BooleanError::EmptyResult) => return Ok(Cover::Misses),",
     "                        Ok(_) => return Err(unsupported(STAND_IN_WITHIN_TOLERANCE)),\n"
     "                        Err(BooleanError::EmptyResult) => return Ok(Cover::Misses),"),
    ("F4 a region is forgotten as soon as it is decided", REGIONS,
     "        let mut known = self.known.borrow_mut();\n",
     "        let mut known = self.known.borrow_mut();\n"
     "        known.clear();\n",
     REGIONS_LIB),
    ("F2c a clip emptied within the tolerance is an identity", CLIP,
     "return if emptiness_is_exact(&subject, &solid, operator) {",
     "return if true {"),
    ("F3a an empty result ignores its readings", DEVIATION,
     "    with_readings(EMPTY_RESULT, DeviationBound::Certified(0.0), report)",
     "    Deviation::one(\n"
     "        DeviationPath::Boolean,\n"
     "        EMPTY_RESULT,\n"
     "        DeviationBound::Certified(0.0),\n"
     "    )"),
    ("F3b an empty result of flattened operands is certified", DEVIATION,
     "    if flattened > 0.0 {\n"
     "        return Deviation::one(",
     "    if false && flattened > 0.0 {\n"
     "        return Deviation::one("),
    ("M7 the mesh swap is not a clip", COMPILER,
     "        (true, false) if operator == Intersection => (right, left),",
     "        (true, false) if false && operator == Intersection => (right, left),"),
]

probe_lib.main(MUTANTS, tests=TESTS, timeout=1800, reject_uncompiled=True)
