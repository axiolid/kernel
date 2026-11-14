"""Mutation probe for the clip frame and the skin guard (#400, ADR 0102).

Each mutant reintroduces a way the cut could fail -- a clip left out of its
subject's frame, its offset dropped or its sign flipped, the guard gated to
unions, the guard reading the tolerance instead of the rounding reach, and
`thin_sheet` blind to the overlap of the pair, to the offset or to the linear
distance -- and must turn a test red. `turned_clip.rs` checks the frame, the
offset, the authored skin and the sweep; `frame::relative` has unit tests;
`axiolid-inspect`'s `sheet.rs` checks the query. `--only NAME` runs only the
mutants whose names contain NAME.

The referee decides the guard's refusal, accepted only within the reach: a
referee that reads above the reach, a skin that never reaches it, and one
that is refereed at every cut of a chain, not once at its root (ADR 0102),
are mutants too.
"""
import probe_lib

COMPILER = "crates/execution/compile/src/compiler.rs"
REFEREE = "crates/execution/compile/src/referee.rs"
FRAME = "crates/execution/compile/src/frame.rs"
SHEET = "crates/algorithms/query/inspect/src/sheet.rs"

T = ["-p", "axiolid-mesh-compile", "--test", "turned_clip"]
L = ["-p", "axiolid-mesh-compile", "--lib"]
F = ["-p", "axiolid-mesh-compile", "--test", "flush_through_cut"]
S = ["-p", "axiolid-inspect", "--test", "sheet"]

GUARD_HEAD = (
    "            axiolid_core::BooleanOperator::Difference | axiolid_core::BooleanOperator::Intersection\n"
    "        ) {\n"
    "            axiolid_inspect::thin_sheet("
)

MUTANTS = [
    ("clips left out of the subject's frame", COMPILER,
     "        let Some(clip) = clip_of(graph, regions, (left, right), operator)? else {\n            return Ok(None);",
     "        let Some(clip) = clip_of(graph, regions, (left, right), operator)?.filter(|_| false) else {\n            return Ok(None);",
     [T]),
    ("relative placement drops the offset", FRAME,
     "        let offset = relative_offset(subject, tool)?;\n        return Some(Transform3::from_translation(offset));",
     "        let _ = relative_offset(subject, tool)?;\n        return Some(Transform3::from_translation(Vec3::ZERO));",
     [T, L]),
    ("relative placement flips the offset's sign", FRAME,
     "        return Some(Transform3::from_translation(offset));",
     "        return Some(Transform3::from_translation(-offset));",
     [T, L]),
    ("the skin guard reads unions only", COMPILER,
     GUARD_HEAD,
     GUARD_HEAD.replace(
         "axiolid_core::BooleanOperator::Difference | axiolid_core::BooleanOperator::Intersection",
         "axiolid_core::BooleanOperator::Union"),
     [T, F]),
    ("the skin guard reads the tolerance, not the rounding reach", COMPILER,
     "axiolid_inspect::thin_sheet(&outcome.mesh, reach)",
     "axiolid_inspect::thin_sheet(&outcome.mesh, crate::snap::contact_reach(options.tolerance()))",
     [T]),
    ("thin_sheet ignores how far the pair overlaps", SHEET,
     "    (area2 > linear * perimeter).then_some(offset)",
     "    let _ = (area2, perimeter);\n    Some(offset)",
     [S]),
    ("thin_sheet reads a sheet only of zero offset", SHEET,
     "    if offset > linear {",
     "    if offset > 0.0 {",
     [S]),
    ("thin_sheet ignores the linear distance", SHEET,
     "    if !linear.is_finite() || linear < 0.0 {",
     "    let linear = linear.max(1e-5);\n    if !linear.is_finite() || linear < 0.0 {",
     [S]),
    ("the referee is accepted above the rounding reach", REFEREE,
     "        if !within_reach(report.linear(), reach) {",
     "        if false && !within_reach(report.linear(), reach) {",
     [T, F]),
    ("a skin never reaches the referee", COMPILER,
     "                skin: Some(skin),\n                ..\n            }) => {",
     "                skin: Some(skin),\n                ..\n            }) if false => {",
     [T, F]),
    ("every skin is refereed where it is cut, not at the chain's root", COMPILER,
     "    if id == root || graph.roots().contains(&id) {\n        return false;\n    }",
     "    if true || id == root || graph.roots().contains(&id) {\n        return false;\n    }",
     [L + ["--", "many_flush_openings_do_not_compile_the_subtree_each"]]),
]

probe_lib.main(MUTANTS, timeout=1800, reject_uncompiled=True)
