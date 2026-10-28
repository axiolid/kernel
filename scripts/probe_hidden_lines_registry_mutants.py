"""Mutation probe for the hidden-line registry (#138, ADR 0090).

Each mutant removes one check the registry owes the contract, or lets a
refusal lose the mesh it is about, and must turn a registry test red.
"""
import probe_lib

R = "crates/execution/dispatch/src/hidden_lines.rs"
TESTS = [
    "-p", "axiolid-dispatch", "--no-default-features", "--features", "mesh-hidden-lines",
    "--test", "mesh_hidden_lines_contract",
]

MUTANTS = [
    ('duplicate ids accepted', R, '        if !ids.insert(entry.id) {', '        if !ids.insert(entry.id) && false {', TESTS),
    ('Solid boundary accepted', R, '        if entry.role == MeshRole::Solid && health.boundary_edges > 0 {', '        if false && entry.role == MeshRole::Solid && health.boundary_edges > 0 {', TESTS),
    ('view change accepted', R, '    if outcome.view != request.view {', '    if false && outcome.view != request.view {', TESTS),
    ('segment order unchecked', R, '        if previous.is_some_and(|last| last > key) {', '        if false && previous.is_some_and(|last| last > key) {', TESTS),
    ('hidden emitted when not requested', R, '            Visibility::Hidden if request.edges.hidden => hidden += 1,', '            Visibility::Hidden => hidden += 1,', TESTS),
    ('same-visibility neighbours accepted', R, '        if request.edges.hidden\n            && previous == Some(key)', '        if false\n            && previous == Some(key)', TESTS),
    ('refusal loses its mesh id', R, '        let refuse = |error: GeomError| HiddenLineError::mesh(entry.id, error);', '        let refuse = |error: GeomError| HiddenLineError::request(error);', TESTS),
    ('refusal naming a foreign mesh accepted', R, '                if !request.meshes.iter().any(|entry| entry.id == *id) {', '                if false && !request.meshes.iter().any(|entry| entry.id == *id) {', TESTS),
    ('edge kind change accepted', R, '        if previous == Some(key) && previous_kind != Some(segment.kind) {', '        if false && previous == Some(key) && previous_kind != Some(segment.kind) {', TESTS),
    ('partial triangle accepted', R, '        if entry.mesh.indices.is_empty() || !entry.mesh.indices.len().is_multiple_of(3) {', '        if entry.mesh.indices.is_empty() {', TESTS),
    ('fallback on any error', R, '                Err(error) => return Err(error),', '                Err(error) => {\n                    last_retryable = Some(error);\n                }', TESTS),
]

probe_lib.main(MUTANTS, timeout=1800, reject_uncompiled=True)
