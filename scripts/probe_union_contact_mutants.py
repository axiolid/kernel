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

Equivalent mutant, deliberately not listed: removing the world-frame
fallback of `build_boolean` (a refused shared-frame cut returning its error
instead of cutting the operands in the world frame). No test configuration
is refused in the shared frame but succeeds in the world frame, because both
frames cut the same operands with the same guard.
"""
import probe_lib

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

MUTANTS = [
    ("guard off", COMPILER,
     "        if operator == axiolid_core::BooleanOperator::Union {\n            let eps",
     "        if false && operator == axiolid_core::BooleanOperator::Union {\n            let eps",
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
     "        && columns == b.matrix3.to_cols_array()\n",
     "        && b.matrix3.determinant().is_finite()\n",
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
]

probe_lib.main(MUTANTS, timeout=1800, reject_uncompiled=True)
