"""Mutation probe for union contact at the linear tolerance and for booleans
cut in their operands' shared frame (#356, #357, ADR 0096).

Each mutant reintroduces a way the contact handling failed or could fail --
the guard off or reading no tolerance, the union snap back at the rounding
or at every operator or at `eps / sqrt 3`, the shared frame off, loose
about the linear parts, dropping the offset, placing by the tool's
transform, unaware of chains, or `surface_gap` blind to an edge through a
face -- and must turn a test red. `union_contact.rs` and
`rotated_face_contact.rs` check the merged and refused unions, `thin_skin.rs`
the skins kept by differences, the library tests the helpers, and the
`axiolid-inspect` tests the clearance. `--only NAME` runs only the mutants
whose names contain NAME.

A guard that only a looser rule would reach is probed together with that
rule: the stretch-by-direction check (`is_similarity`) and the world-
tolerance guard on the placed result back each other up, so removing one
alone is equivalent (a unit test pins the first, and the other is killed
with the first removed); likewise the world-frame fallback of
`build_boolean` is only reached by a refused shared-frame cut, which the
similarity check keeps out of reach, so it is probed with that check
removed.

Not probed: the `placement_retry` exclusion in `frame_of` (an instance whose
source is re-meshed for its placement's rounding is left to the world
frame). Dropping it changes only how finely the shared frame meshes such a
source, which the bound of a measured boolean would show, but that bound is
unbounded for the boolean of curved operands that retry needs, so no test
can see it.
"""
import probe_lib
from probe_lib import Mutant, ROOT

COMPILER = "crates/execution/compile/src/compiler.rs"
CONTACT = "crates/execution/compile/src/contact.rs"
SNAP = "crates/execution/compile/src/snap.rs"
FRAME = "crates/execution/compile/src/frame.rs"
CLEARANCE = "crates/algorithms/query/inspect/src/clearance.rs"

U = ["-p", "axiolid-mesh-compile", "--test", "union_contact"]
R = ["-p", "axiolid-mesh-compile", "--test", "rotated_face_contact"]
S = ["-p", "axiolid-mesh-compile", "--test", "thin_skin"]
L = ["-p", "axiolid-mesh-compile", "--lib"]
I = ["-p", "axiolid-inspect"]
T = ["-p", "axiolid-mesh-compile", "--test", "surface_models"]
PR = ["-p", "axiolid-mesh-compile", "--test", "placement_rounding"]

def edits(*changes):
    return [(ROOT / path, old, new) for path, old, new in changes]


SIMILAR = (FRAME, "        && is_similarity(t)\n", "        && true\n")
FALLBACK = (
    COMPILER,
    "        if let Ok(Some(built)) =\n            self.boolean_in_shared_frame(graph, (left, right), operator, options, budget, cache)\n        {",
    "        if let Some(built) =\n            self.boolean_in_shared_frame(graph, (left, right), operator, options, budget, cache)?\n        {",
)
GUARD = (
    COMPILER,
    "            if let Some(found) = crate::contact::unmerged(&placed.mesh, eps) {",
    "            if let Some(found) = crate::contact::unmerged(&placed.mesh, eps).filter(|_| false) {",
)

MUTANTS = [
    ("guard off", COMPILER,
     "        if operator == axiolid_core::BooleanOperator::Union {\n            let eps = options.tolerance().linear();\n            if let Some(found) = crate::contact::unmerged(&outcome",
     "        if false && operator == axiolid_core::BooleanOperator::Union {\n            let eps = options.tolerance().linear();\n            if let Some(found) = crate::contact::unmerged(&outcome",
     [U, R, L]),
    ("guard reads no tolerance", COMPILER,
     "crate::contact::unmerged(&outcome.mesh, eps)",
     "crate::contact::unmerged(&outcome.mesh, 0.0)",
     [U]),
    ("guard skips the first pair of bodies", CONTACT,
     "for j in i + 1..parts.len()",
     "for j in i + 2..parts.len()",
     [L, U]),
    ("union reach back at the rounding", COMPILER,
     "            crate::snap::contact_reach(options.tolerance())\n",
     "            crate::snap::reach(&subject.mesh, tool, options.tolerance(), placed)\n",
     [U, L]),
    ("contact reach for every operator", COMPILER,
     "let reach = if operator == axiolid_core::BooleanOperator::Union {",
     "let reach = if true {",
     [S, U]),
    ("contact reach divided by sqrt 3", SNAP,
     "    if eps > 0.0 && eps.is_finite() {\n        eps\n",
     "    if eps > 0.0 && eps.is_finite() {\n        eps / 3.0_f64.sqrt()\n",
     [U]),
    ("shared frame off", COMPILER,
     "        if clipping && clip_tool(graph, right).is_some() {\n            return Ok(None);",
     "        if true || clipping && clip_tool(graph, right).is_some() {\n            return Ok(None);",
     [R]),
    ("linear parts compared loosely", FRAME,
     "    similarity(a) && a.matrix3.to_cols_array() == b.matrix3.to_cols_array()\n",
     "    similarity(a) && b.matrix3.determinant().is_finite()\n",
     [R, L]),
    ("tool offset dropped", FRAME,
     "    if same {\n        return Some(Vec3::ZERO);",
     "    if true {\n        return Some(Vec3::ZERO);",
     [R, L]),
    ("result placed by the tool's transform", COMPILER,
     "channels::transform(&local, subject_frame, rounding)",
     "channels::transform(&local, tool_frame, rounding)",
     [R]),
    ("chains not framed", COMPILER,
     "            .map(|local| (local.transform, &local.built)))",
     "            .map(|local| (local.transform, &local.built))\n            .filter(|_| false))",
     [R]),
    ("surface_gap blind to an edge through a face", CLEARANCE,
     "    if edge_crosses(a, b) || edge_crosses(b, a) {",
     "    if false && (edge_crosses(a, b) || edge_crosses(b, a)) {",
     [I]),
    ("stretch by direction read as a shared frame", FRAME,
     "        && is_similarity(t)\n", "        && true\n", [L]),
    Mutant("stretched frame let through without the world guard",
           edits(SIMILAR, GUARD),
           [["cargo", "test", "-q", *U]]),
    Mutant("refused shared-frame cut returned, stretch let through",
           edits(SIMILAR, FALLBACK),
           [["cargo", "test", "-q", *U]]),
    ("shared frame cut at the world tolerance", COMPILER,
     ".with_tolerance(local_budget.tolerance)",
     ".with_tolerance(options.tolerance())", [U]),
    ("shared frame takes surface operands", COMPILER,
     '        refuse_surface_operand(subject, "subject")?;\n        refuse_surface_operand(tool, "tool")?;\n',
     "        let _ = (subject, tool);\n", [T]),
    ("union snap moves inward", SNAP,
     "            if !whole && outward[k][i] * (face.at - p[k]) < 0.0 || !whole && outward[k][i] == 0.0 {",
     "            if !whole && outward[k][i] == 0.0 {", [U]),
    ("union snap moves a body away from the contact", SNAP,
     "        if whole {\n            key(v)\n        } else {\n            patch[k][i] as u64\n        }",
     "        key(v)", [U]),
]

probe_lib.main(MUTANTS, timeout=1800, reject_uncompiled=True)
