"""Mutation probe for the hidden-line registry (#138, ADR 0090).

Each mutant removes one check the registry owes the contract, or lets a
refusal lose the mesh it is about, and must turn a registry test red.
"""
import probe_lib

R = "crates/execution/dispatch/src/hidden_lines.rs"
V = "crates/execution/dispatch/src/view_checks.rs"
TESTS = [
    "-p", "axiolid-dispatch", "--no-default-features", "--features", "mesh-hidden-lines",
    "--test", "mesh_hidden_lines_contract",
]

LOOP = '        for entry in request.meshes {\n            options.check_cancelled()?;\n            if let Some(error) = mesh_verdict'

MUTANTS = [
    ('duplicate ids accepted', V, '        if !ids.insert(entry.id) {', '        if !ids.insert(entry.id) && false {', TESTS),
    ('Solid boundary accepted', V, '    if entry.role == MeshRole::Solid && health.boundary_edges > 0 {', '    if false && entry.role == MeshRole::Solid && health.boundary_edges > 0 {', TESTS),
    ('view change accepted', R, '    if outcome.view != request.view {', '    if false && outcome.view != request.view {', TESTS),
    ('segment order unchecked', R, '        if previous.is_some_and(|last| last > key) {', '        if false && previous.is_some_and(|last| last > key) {', TESTS),
    ('hidden emitted when not requested', R, '            Visibility::Hidden if request.edges.hidden => hidden += 1,', '            Visibility::Hidden => hidden += 1,', TESTS),
    ('same-visibility neighbours accepted', R, '        if request.edges.hidden\n            && previous == Some(key)', '        if false\n            && previous == Some(key)', TESTS),
    ('refusal loses its mesh id', V, '        if let Some(error) = mesh_verdict(entry, groups, options)? {\n            return Err(HiddenLineError::mesh(entry.id, error));', '        if let Some(error) = mesh_verdict(entry, groups, options)? {\n            return Err(HiddenLineError::request(error));', TESTS),
    ('refusal naming a foreign mesh accepted', R, '                if !request.meshes.iter().any(|entry| entry.id == *id) {', '                if false && !request.meshes.iter().any(|entry| entry.id == *id) {', TESTS),
    ('edge kind change accepted', R, '        if previous == Some(key) && previous_kind != Some(segment.kind) {', '        if false && previous == Some(key) && previous_kind != Some(segment.kind) {', TESTS),
    ('an unselected group boundary accepted', R, '        EdgeKind::GroupBoundary => request.edges.group_boundary,', '        EdgeKind::GroupBoundary => true,', TESTS),
    ('a selected group boundary refused', R, '        EdgeKind::GroupBoundary => request.edges.group_boundary,', '        EdgeKind::GroupBoundary => false,', TESTS),
    ('hidden_lines does not read group lengths', R, '            request.edges.group_boundary.then_some(LABEL),\n            options,', '            None,\n            options,', TESTS),
    ('hidden_lines reads group lengths when group boundaries are not selected', R, '            request.edges.group_boundary.then_some(LABEL),\n            options,', '            Some(LABEL),\n            options,', TESTS),
    ('the audit does not read group lengths', R, '        let groups = request.edges.group_boundary.then_some(LABEL);', '        let groups: Option<&str> = None;', TESTS),
    ('the verdict ignores a short group slice', V, '    if let Some(fault) = groups.and_then(|label| group_length_fault(entry, label)) {', '    if let Some(fault) = groups.and_then(|label| group_length_fault(entry, label)).filter(|_| false) {', TESTS),
    ('partial triangle accepted', V, '        if entry.mesh.indices.is_empty() || !entry.mesh.indices.len().is_multiple_of(3) {', '        if entry.mesh.indices.is_empty() {', TESTS),
    ('fallback on any error', R, '                Err(error) => return Err(error),', '                Err(error) => {\n                    last_retryable = Some(error);\n                }', TESTS),
    ('audit stops at the first refusal', R, '                refused.push(RefusedMesh::new(entry.id, error));', '                refused.push(RefusedMesh::new(entry.id, error));\n                break;', TESTS),
    ('audit skips the first mesh', R, LOOP, LOOP.replace('request.meshes {', 'request.meshes.iter().skip(1) {'), TESTS),
    ('audit skips the last mesh', R, LOOP, LOOP.replace('request.meshes {', 'request.meshes.iter().rev().skip(1).rev() {'), TESTS),
    ('audit names the wrong mesh', R, 'refused.push(RefusedMesh::new(entry.id, error))', 'refused.push(RefusedMesh::new(request.meshes[0].id, error))', TESTS),
    ('audit lists every mesh', R, 'if let Some(error) = mesh_verdict(entry, groups, options)? {\n                refused', 'if let Some(error) = mesh_verdict(entry, groups, options)?.or(Some(GeomError::Cancelled)) {\n                refused', TESTS),
    ('audit skips the shape check', R, '        validate_request_shape(request, options)?;\n        if let Some(budget)', '        if let Some(budget)', TESTS),
    ('audit ignores the memory budget', R, '            if !fits {', '            if false && !fits {', TESTS),
    ('audit reversed', R, '        Ok(refused)', '        refused.reverse();\n        Ok(refused)', TESTS),
]

probe_lib.main(MUTANTS, timeout=1800, reject_uncompiled=True)
