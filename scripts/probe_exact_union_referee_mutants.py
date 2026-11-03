"""Mutation probe for the exact referee of unions (#391, ADR 0098).

Each mutant reintroduces a way the referee could fail -- never asked, asked
always or for contact it need not decide, its refusal turned into the
caller's, its refusal not named, its details swapped or its coincident
reading dropped, its channels kept, its result measured again, the exact
B-rep read by the graph's domains or against its intervals' direction, and
`near_coplanar` blind to the exactness test, to one of the two planes or to
the distance -- and must turn a test red. `exact_union_referee.rs` checks
the routing and the reports, `union_contact.rs` the flipped refusal, the
`referee` unit test the channels, `brep/exact_tests.rs` the tessellator,
and the `axiolid-inspect` tests `near_coplanar`. `--only NAME` runs only the
mutants whose names contain NAME.

Not probed: the referee's closed-solid and unmerged-contact guards. They are
defence in depth: the exact compiler reads contact at the same tolerance the
guard measures, so no exact union that tessellates open, or holds two bodies
within the tolerance, exists to reach them. The pinch guard is probed: an
edge contact is exactly what it refuses.
"""
import probe_lib

COMPILER = "crates/execution/compile/src/compiler.rs"
REFEREE = "crates/execution/compile/src/referee.rs"
VIEW = "crates/execution/compile/src/brep/view.rs"
COPLANAR = "crates/algorithms/query/inspect/src/coplanar.rs"

E = ["-p", "axiolid-mesh-compile", "--test", "exact_union_referee"]
U = ["-p", "axiolid-mesh-compile", "--test", "union_contact"]
L = ["-p", "axiolid-mesh-compile", "--lib"]
I = ["-p", "axiolid-inspect", "--test", "coplanar"]

MUTANTS = [
    ("referee never asked", COMPILER,
     "        if operator != axiolid_core::BooleanOperator::Union {\n            return cut",
     "        if true {\n            return cut",
     [E, U]),
    ("unsettled contact never read", COMPILER,
     "            && axiolid_inspect::near_coplanar(cut_subject, cut_tool, eps).is_some();",
     "            && false;",
     [E]),
    ("every union is unsettled", COMPILER,
     "            && axiolid_inspect::near_coplanar(cut_subject, cut_tool, eps).is_some();",
     "            && (axiolid_inspect::near_coplanar(cut_subject, cut_tool, eps).is_some() || true);",
     [E]),
    ("unsettled read on the given operands, not the pair cut", COMPILER,
     "        let snapped_pair = moved > 0.0;",
     "        let snapped_pair = false;",
     [E, U]),
    ("a shared frame ignores the sources' composed turns", COMPILER,
     "            if !crate::frame::shared_linear(a, b) {",
     "            if false && !crate::frame::shared_linear(a, b) {",
     [E]),
    ("contact refusals never reach the referee", COMPILER,
     "            Err(CutRefusal::Contact(error)) => match referee(&mut reason)? {",
     "            Err(CutRefusal::Contact(error)) => match { let _ = &referee; None::<Built> } {",
     [E, U]),
    ("a refused referee turns a valid mesh result into an error", COMPILER,
     "                None => (cut.built, Decided::Mesh),",
     "                None => return Err(GeomError::Degenerate(reason.unwrap_or_default())),",
     [E]),
    ("a contact refusal does not name the referee's reason", COMPILER,
     '                        "{text}; the exact referee refused: {}",\n                        reason.unwrap_or_default()',
     '                        "{text}{}",\n                        reason.map(|_| "").unwrap_or_default()',
     [E, U]),
    ("a contact refusal loses its own message", COMPILER,
     '                        "{text}; the exact referee refused: {}",\n                        reason.unwrap_or_default()',
     '                        "the exact referee refused: {}",\n                        reason.unwrap_or_default()',
     [E]),
    ("a refereed union is measured again", COMPILER,
     "Measure::Deviation(_) | Measure::Skip if decided == Decided::Referee => built,",
     "Measure::Deviation(_) | Measure::Skip if decided == Decided::Mesh && false => built,",
     [E]),
    ("referee drops the coincident reading", REFEREE,
     "        .find(|decision| decision.kind == ToleranceDecisionKind::CoincidentSupports)",
     "        .find(|_| false)",
     [E]),
    ("referee swaps its two details", REFEREE,
     "    let detail = if report.is_exact() {",
     "    let detail = if !report.is_exact() {",
     [E]),
    ("referee reports no bound of its own", REFEREE,
     "Deviation::one(DeviationPath::Boolean, detail, tessellation.worst_bound())",
     "Deviation::one(DeviationPath::Boolean, detail, DeviationBound::Certified(0.0))",
     [E, L]),
    ("referee keeps its operands' channels", REFEREE,
     "AttributeFate::Dropped(DropReason::ProviderLimitation)))",
     "AttributeFate::Preserved))",
     [L]),
    ("referee accepts a union that touches itself", REFEREE,
     "    if crate::pinch::find(&mesh).is_some() {",
     "    if crate::pinch::find(&mesh).is_some() && false {",
     [E]),
    ("exact view reads the graph's domains", VIEW,
     "                        let span = exact\n                            .pcurve_interval(loop_id, position)\n                            .ok_or_else(|| missing(\"pcurve interval\"))?;",
     "                        let span = axiolid_reference::curve::domain2(curve);",
     [L]),
    ("exact view ignores the interval's direction", VIEW,
     "                            .ok_or_else(|| missing(\"pcurve interval\"))?;",
     "                            .ok_or_else(|| missing(\"pcurve interval\"))\n                            .map(|s| Interval::new(s.start.min(s.end), s.start.max(s.end)))?;",
     [L]),
    ("exact view ranks pcurves by use, not handle", VIEW,
     "                            key: ranks[&id],\n                        })\n                    })\n                    .transpose()?;\n                Ok(EdgeUse {\n                    edge: use_.edge,\n                    orientation: use_.orientation,\n                    pcurve,\n                })\n            })\n            .collect::<GeomResult<Vec<_>>>()?;\n        view.add_loop(Loop { edges });\n    }\n    for face in brep.faces() {\n        let surface = face",
     "                            key: position as u32,\n                        })\n                    })\n                    .transpose()?;\n                Ok(EdgeUse {\n                    edge: use_.edge,\n                    orientation: use_.orientation,\n                    pcurve,\n                })\n            })\n            .collect::<GeomResult<Vec<_>>>()?;\n        view.add_loop(Loop { edges });\n    }\n    for face in brep.faces() {\n        let surface = face",
     [L]),
    ("near_coplanar reports exact contact", COPLANAR,
     "    (!exactly_on_plane).then_some(offset)",
     "    let _ = exactly_on_plane;\n    Some(offset)",
     [I]),
    ("near_coplanar reads one plane only", COPLANAR,
     "    let offset = to_a.max(to_b);",
     "    let offset = to_a;",
     [I]),
    ("near_coplanar reads the other plane only", COPLANAR,
     "    let offset = to_a.max(to_b);",
     "    let offset = to_b;",
     [I]),
    ("near_coplanar ignores the distance between triangles", COPLANAR,
     "    if offset > linear || triangle_distance(a, b) > linear {",
     "    if offset > linear {",
     [I]),
    ("near_coplanar ignores the offset", COPLANAR,
     "    if offset > linear || triangle_distance(a, b) > linear {",
     "    if triangle_distance(a, b) > linear {",
     [I]),
    ("near_coplanar reports a later pair first", COPLANAR,
     "        seconds.sort_unstable();",
     "        seconds.sort_unstable_by(|a, b| b.cmp(a));",
     [I]),
    ("near_coplanar reads a degenerate triangle", COPLANAR,
     "    if !length.is_finite() || length <= 0.0 {",
     "    if !length.is_finite() {",
     [I]),
]

probe_lib.main(MUTANTS, timeout=1800, reject_uncompiled=True)
