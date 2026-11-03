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

LOOP = '        for entry in request.meshes {\n            options.check_cancelled()?;\n            if let Some(error) = mesh_verdict'

MUTANTS = [
    ('duplicate ids accepted', R, '        if !ids.insert(entry.id) {', '        if !ids.insert(entry.id) && false {', TESTS),
    ('Solid boundary accepted', R, '    if entry.role == MeshRole::Solid && health.boundary_edges > 0 {', '    if false && entry.role == MeshRole::Solid && health.boundary_edges > 0 {', TESTS),
    ('view change accepted', R, '    if outcome.view != request.view {', '    if false && outcome.view != request.view {', TESTS),
    ('segment order unchecked', R, '        if previous.is_some_and(|last| last > key) {', '        if false && previous.is_some_and(|last| last > key) {', TESTS),
    ('hidden emitted when not requested', R, '            Visibility::Hidden if request.edges.hidden => hidden += 1,', '            Visibility::Hidden => hidden += 1,', TESTS),
    ('same-visibility neighbours accepted', R, '        if request.edges.hidden\n            && previous == Some(key)', '        if false\n            && previous == Some(key)', TESTS),
    ('refusal loses its mesh id', R, '            return Err(HiddenLineError::mesh(entry.id, error));', '            return Err(HiddenLineError::request(error));', TESTS),
    ('refusal naming a foreign mesh accepted', R, '                if !request.meshes.iter().any(|entry| entry.id == *id) {', '                if false && !request.meshes.iter().any(|entry| entry.id == *id) {', TESTS),
    ('edge kind change accepted', R, '        if previous == Some(key) && previous_kind != Some(segment.kind) {', '        if false && previous == Some(key) && previous_kind != Some(segment.kind) {', TESTS),
    ('partial triangle accepted', R, '        if entry.mesh.indices.is_empty() || !entry.mesh.indices.len().is_multiple_of(3) {', '        if entry.mesh.indices.is_empty() {', TESTS),
    ('fallback on any error', R, '                Err(error) => return Err(error),', '                Err(error) => {\n                    last_retryable = Some(error);\n                }', TESTS),
    ('audit stops at the first refusal', R, '                refused.push(RefusedMesh::new(entry.id, error));', '                refused.push(RefusedMesh::new(entry.id, error));\n                break;', TESTS),
    ('audit skips the first mesh', R, LOOP, LOOP.replace('request.meshes {', 'request.meshes.iter().skip(1) {'), TESTS),
    ('audit skips the last mesh', R, LOOP, LOOP.replace('request.meshes {', 'request.meshes.iter().rev().skip(1).rev() {'), TESTS),
    ('audit names the wrong mesh', R, 'refused.push(RefusedMesh::new(entry.id, error))', 'refused.push(RefusedMesh::new(request.meshes[0].id, error))', TESTS),
    ('audit lists every mesh', R, 'if let Some(error) = mesh_verdict(entry, options)? {\n                refused', 'if let Some(error) = mesh_verdict(entry, options)?.or(Some(GeomError::Cancelled)) {\n                refused', TESTS),
    ('audit skips the shape check', R, '        validate_request_shape(request, options)?;\n        if let Some(budget)', '        if let Some(budget)', TESTS),
    ('audit ignores the memory budget', R, '            if !fits {', '            if false && !fits {', TESTS),
    ('audit reversed', R, '        Ok(refused)', '        refused.reverse();\n        Ok(refused)', TESTS),
]

probe_lib.main(MUTANTS, timeout=1800, reject_uncompiled=True)
