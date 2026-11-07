"""Mutation probe for the C ABI of hidden lines (#473).

Each mutant relaxes one boundary check, crosses one mapping or drops one
piece of recorded context in `axiolid-capi`'s `drawing.rs` (or the entry
points it shares) and must turn `tests/drawing.rs` red.

Not a mutant: the `Unsupported` mapping of a refusal for a capability the
provider lacks. The portable provider never refuses that way, so no test can
provoke it; the budget mapping beside it is probed.

The short-capacity mutant also truncates the copy: dropping only the check
would let the copy write past the test's seven-slot buffer.
"""
import probe_lib

DRAWING = "crates/facade/axiolid-capi/src/drawing.rs"
LIB = "crates/facade/axiolid-capi/src/lib.rs"

TEST = ["cargo", "test", "-q", "-p", "axiolid-capi", "--test", "drawing"]

MUTANTS = [
    ('finiteness is not checked', DRAWING,
     '    if values.iter().all(|value| value.is_finite()) {',
     '    if true {', [TEST]),
    ('has_near is ignored', DRAWING,
     '    if has_near {\n        result = result.with_near(view.near);',
     '    if true {\n        result = result.with_near(view.near);', [TEST]),
    ('has_far is ignored', DRAWING,
     '    if has_far {\n        result = result.with_far(view.far);',
     '    if true {\n        result = result.with_far(view.far);', [TEST]),
    ('flags are not checked', DRAWING,
     '        1 => Ok(true),\n        _ => Err(Refusal::invalid(message)),',
     '        _ => Ok(true),', [TEST]),
    ('view reserved byte accepted', DRAWING,
     '    if view.reserved != [0; 2] {', '    if false {', [TEST]),
    ('edge reserved byte accepted', DRAWING,
     '    if edges.reserved != [0; 3] {', '    if false {', [TEST]),
    ('unknown projection accepted', DRAWING,
     '    if view.projection != AXIOLID_PROJECTION_ORTHOGRAPHIC {', '    if false {', [TEST]),
    ('unknown mesh role is a sheet', DRAWING,
     '        _ => {\n            return Err(Refusal::mesh(\n                AxiolidStatus::InvalidArgument,\n                entry.id,\n                "unknown mesh role",\n            ))\n        }',
     '        _ => MeshRole::Sheet,', [TEST]),
    ('group count with a single mode accepted', DRAWING,
     '            if entry.group_count != 0 {', '            if false {', [TEST]),
    ('planar mode draws single groups', DRAWING,
     '                    AXIOLID_GROUP_MODE_PLANAR => input.with_planar_groups(),',
     '                    AXIOLID_GROUP_MODE_PLANAR => input,', [TEST]),
    ('mesh count not limited by the context', DRAWING,
     '        if mesh_count > context.config.max_meshes as usize {', '        if false {', [TEST]),
    ('region ignored', DRAWING,
     '            let region = region.map(read_region).transpose()?;',
     '            let region = None;', [TEST]),
    ('copy writes despite a short capacity', DRAWING,
     [('            if capacity < outcome.segments.len() {', '            if false {'),
      ('            for (index, segment) in outcome.segments.iter().enumerate() {',
       '            for (index, segment) in outcome.segments.iter().enumerate().take(capacity) {')],
     [TEST]),
    ('refusal drops the mesh', DRAWING,
     '        message: error.to_string(),\n        mesh,\n', '        message: error.to_string(),\n        mesh: None,\n', [TEST]),
    ('budget refusal is a plain failure', DRAWING,
     '        Some(GeomError::BudgetExceeded { .. }) => AxiolidStatus::LimitExceeded,',
     '        Some(GeomError::BudgetExceeded { .. }) => AxiolidStatus::OperationFailed,', [TEST]),
    ('silhouette reported as a crease', DRAWING,
     '        EdgeKind::Silhouette => AxiolidEdgeKind::Silhouette,',
     '        EdgeKind::Silhouette => AxiolidEdgeKind::Crease,', [TEST]),
    ('hidden reported as visible', DRAWING,
     '        Visibility::Hidden => AxiolidVisibility::Hidden,',
     '        Visibility::Hidden => AxiolidVisibility::Visible,', [TEST]),
    ('last error mesh reports no mesh', DRAWING,
     '            out_has_mesh.write(u8::from(error.mesh.is_some()));',
     '            out_has_mesh.write(0);', [TEST]),
    ('capability not listed', LIB,
     '                    | axiolid::contracts::Operation::MeshHiddenLines\n', '', [TEST]),
    ('operation id unmapped', LIB,
     '        Operation::MeshHiddenLines => AxiolidOperation::MeshHiddenLines,',
     '        Operation::MeshHiddenLines => AxiolidOperation::Unknown,', [TEST]),
]

probe_lib.main(MUTANTS, timeout=1800, command=(), reject_uncompiled=True)
