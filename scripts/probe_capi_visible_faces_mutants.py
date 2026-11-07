"""Mutation probe for the C ABI of visible faces (#386).

Each mutant offsets a table, swaps two rings, evidence fields or limits,
crosses a mapping, relaxes a capacity check or drops recorded context in
`axiolid-capi`'s `faces.rs` (or the entry points it shares) and must turn
`tests/faces.rs` red.

Not a mutant: the `Unsupported` mapping of a refusal for a capability the
provider lacks. The portable provider never refuses that way, so no test can
provoke it; the budget mapping beside it is probed.

The snap fields are probed through a thin-occluder scene and a seeded
isometric scene (`rerouted_edges`, `collapsed_regions` and `collapsed_holes`
are told apart by their different values there).
"""
import probe_lib

FACES = "crates/facade/axiolid-capi/src/faces.rs"
DRAWING = "crates/facade/axiolid-capi/src/drawing.rs"
LIB = "crates/facade/axiolid-capi/src/lib.rs"
OPERATIONS = "crates/facade/axiolid-capi/src/operations.rs"

TEST = ["cargo", "test", "-q", "-p", "axiolid-capi", "--test", "faces"]


def swap(name, first, second):
    """Evidence mutant: two fields read each other's source value."""
    a = f"                {first}: evidence.{first},\n                {second}: evidence.{second},"
    b = f"                {first}: evidence.{second},\n                {second}: evidence.{first},"
    return (name, FACES, a, b, [TEST])


def zero(name, field, value="0"):
    """Evidence mutant: one field is hard-wired to zero."""
    source = f"                {field}: evidence.{field}" + (".get()," if field == "snap_displacement" else ",")
    return (name, FACES, source, f"                {field}: {value},", [TEST])


def cross(name, field, other):
    """Evidence mutant: one field reads another's source value."""
    return (name, FACES, f"                {field}: evidence.{field},", f"                {field}: evidence.{other},", [TEST])


MUTANTS = [
    # Offsets and ring order.
    ('first ring offset is one too high', FACES,
     '                    first_ring: ring_at,', '                    first_ring: ring_at + 1,', [TEST]),
    ('first point offset is one too high', FACES,
     '                        first_point: point_at,', '                        first_point: point_at + 1,', [TEST]),
    ('first triangle offset is one too high', FACES,
     '                    first_triangle: triangle_at,', '                    first_triangle: triangle_at + 1,', [TEST]),
    ('outer and hole rings swapped', FACES,
     'std::iter::once(&region.outer).chain(&region.holes)',
     'region.holes.iter().chain(std::iter::once(&region.outer))', [TEST]),
    ('ring count forgets the outer ring', FACES,
     '                    ring_count: 1 + region.holes.len(),', '                    ring_count: region.holes.len(),', [TEST]),
    ('triangle count is the ring count', FACES,
     '                    triangle_count: region.triangles.len(),', '                    triangle_count: region.holes.len(),', [TEST]),
    ('point x and y swapped', FACES,
     '                            points_xy.add(2 * point_at).write(point.x);\n                            points_xy.add(2 * point_at + 1).write(point.y);',
     '                            points_xy.add(2 * point_at).write(point.y);\n                            points_xy.add(2 * point_at + 1).write(point.x);', [TEST]),
    ('point y shifted', FACES,
     '                            points_xy.add(2 * point_at + 1).write(point.y);',
     '                            points_xy.add(2 * point_at + 1).write(point.y + 1.0);', [TEST]),
    # Counts.
    ('counts forget the outer rings', FACES,
     '        counts.rings += 1 + region.holes.len();', '        counts.rings += region.holes.len();', [TEST]),
    ('counts forget the hole points', FACES,
     '        counts.points += region.outer.len() + region.holes.iter().map(Vec::len).sum::<usize>();',
     '        counts.points += region.outer.len();', [TEST]),
    ('counts forget the triangles', FACES,
     '        counts.triangles += region.triangles.len();', '        counts.triangles += 0;', [TEST]),
    ('region count is the ring count', FACES,
     '        regions: outcome.regions.len(),', '        regions: outcome.regions.iter().map(|r| 1 + r.holes.len()).sum(),', [TEST]),
    # Capacities and null pointers.
    ('copy ignores a short region capacity', FACES, 'region_cap < need.regions', 'false', [TEST]),
    ('copy ignores a short ring capacity', FACES, 'ring_cap < need.rings', 'false', [TEST]),
    ('copy ignores a short point capacity', FACES, 'point_cap < need.points', 'false', [TEST]),
    ('copy ignores a short triangle capacity', FACES, 'triangle_cap < need.triangles', 'false', [TEST]),
    ('copy refuses a null table with zero capacity', FACES,
     '(region_cap != 0 && regions.is_null())', 'regions.is_null()', [TEST]),
    ('copy accepts a null ring table', FACES,
     '|| (ring_cap != 0 && rings.is_null())', '|| false', [TEST]),
    ('copy accepts a null point table', FACES,
     '|| (point_cap != 0 && points_xy.is_null())', '|| false', [TEST]),
    ('copy accepts a null triangle table', FACES,
     '|| (triangle_cap != 0 && triangles.is_null())', '|| false', [TEST]),
    # Region mappings.
    ('group reported as a cap', FACES,
     '        FaceSource::Group(group) => (AxiolidFaceSource::Group, group.0),',
     '        FaceSource::Group(group) => (AxiolidFaceSource::Cap, group.0),', [TEST]),
    ('cap reported as a group', FACES,
     '        FaceSource::Cap => (AxiolidFaceSource::Cap, 0),',
     '        FaceSource::Cap => (AxiolidFaceSource::Group, 0),', [TEST]),
    ('group id dropped', FACES,
     '        FaceSource::Group(group) => (AxiolidFaceSource::Group, group.0),',
     '        FaceSource::Group(_) => (AxiolidFaceSource::Group, 0),', [TEST]),
    ('cap carries a group id', FACES,
     '        FaceSource::Cap => (AxiolidFaceSource::Cap, 0),',
     '        FaceSource::Cap => (AxiolidFaceSource::Cap, 1),', [TEST]),
    ('mesh id dropped', FACES,
     '                    mesh_id: region.mesh.0,', '                    mesh_id: 0,', [TEST]),
    ('result kind is a mesh', OPERATIONS,
     '            StoredResult::VisibleFaces(_) => AxiolidGeometryKind::VisibleFaces,',
     '            StoredResult::VisibleFaces(_) => AxiolidGeometryKind::TriangleMesh,', [TEST]),
    # Evidence.
    swap('evidence meshes and triangles swapped', 'source_meshes', 'source_triangles'),
    swap('evidence edge-on and removed swapped', 'edge_on_triangles', 'removed_triangles'),
    swap('evidence boundary and intersection swapped', 'boundary_segments', 'intersection_segments'),
    swap('evidence crossings and pieces swapped', 'crossings', 'pieces'),
    swap('evidence tied sides and pairs swapped', 'tied_sides', 'candidate_pairs'),
    swap('evidence pairs and escalations swapped', 'candidate_pairs', 'exact_escalations'),
    swap('evidence regions and cap regions swapped', 'regions', 'cap_regions'),
    cross('evidence meshes read the crossings', 'source_meshes', 'crossings'),
    cross('evidence intersections read the cap regions', 'intersection_segments', 'cap_regions'),
    swap('evidence rerouted and collapsed regions swapped', 'rerouted_edges', 'collapsed_regions'),
    swap('evidence collapsed regions and holes swapped', 'collapsed_regions', 'collapsed_holes'),
    zero('evidence rerouted edges hard-wired to zero', 'rerouted_edges'),
    zero('evidence collapsed regions hard-wired to zero', 'collapsed_regions'),
    zero('evidence collapsed holes hard-wired to zero', 'collapsed_holes'),
    zero('evidence snap displacement hard-wired to zero', 'snap_displacement', '0.0'),
    # Limits.
    ('limits: regions and vertices swapped', FACES,
     '            limits.max_regions,\n            limits.max_region_vertices,',
     '            limits.max_region_vertices,\n            limits.max_regions,', [TEST]),
    ('limits: positions and triangles swapped', FACES,
     '            limits.max_source_vertices,\n            limits.max_source_triangles,',
     '            limits.max_source_triangles,\n            limits.max_source_vertices,', [TEST]),
    ('limits: meshes and pairs swapped', FACES,
     '            limits.max_meshes,\n            limits.max_source_vertices,',
     '            limits.max_candidate_pairs,\n            limits.max_source_vertices,', [TEST]),
    ('limits: pairs ignored', FACES,
     '            limits.max_candidate_pairs,\n            limits.max_regions,',
     '            usize::MAX,\n            limits.max_regions,', [TEST]),
    # Boundary.
    ('mesh count not limited by the context', FACES,
     '        if mesh_count > context.config.max_meshes as usize {', '        if false {', [TEST]),
    ('region ignored', FACES,
     '            let region = region.map(read_region).transpose()?;', '            let region = None;', [TEST]),
    ('planar mode draws single groups', DRAWING,
     '            AXIOLID_GROUP_MODE_PLANAR => input.with_planar_groups(),',
     '            AXIOLID_GROUP_MODE_PLANAR => input,', [TEST]),
    ('group count not checked against the mesh', DRAWING,
     '            if entry.group_count != mesh.indices.len() / 3 {', '            if false {', [TEST]),
    ('view reserved byte accepted', DRAWING,
     '    if view.reserved != [0; 2] {', '    if false {', [TEST]),
    ('unknown projection accepted', DRAWING,
     '    if view.projection != AXIOLID_PROJECTION_ORTHOGRAPHIC {', '    if false {', [TEST]),
    ('unknown mesh role is a sheet', DRAWING,
     '        _ => {\n            return Err(Refusal::mesh(\n                AxiolidStatus::InvalidArgument,\n                entry.id,\n                "unknown mesh role",\n            ))\n        }',
     '        _ => MeshRole::Sheet,', [TEST]),
    # Errors.
    ('boundary refusal recorded as hidden lines', DRAWING,
     '        status: refusal.status,\n        operation,', '        status: refusal.status,\n        operation: AxiolidOperation::MeshHiddenLines,', [TEST]),
    ('contract refusal recorded as hidden lines', DRAWING,
     '        status,\n        operation,', '        status,\n        operation: AxiolidOperation::MeshHiddenLines,', [TEST]),
    ('refusal drops the mesh', DRAWING,
     '        message: error.to_string(),\n        mesh,\n', '        message: error.to_string(),\n        mesh: None,\n', [TEST]),
    ('budget refusal is a plain failure', DRAWING,
     '        Some(GeomError::BudgetExceeded { .. }) => AxiolidStatus::LimitExceeded,',
     '        Some(GeomError::BudgetExceeded { .. }) => AxiolidStatus::OperationFailed,', [TEST]),
    # Capability.
    ('capability not listed', LIB,
     '                    | axiolid::contracts::Operation::MeshVisibleFaces\n', '', [TEST]),
    ('operation id unmapped', LIB,
     '        Operation::MeshVisibleFaces => AxiolidOperation::MeshVisibleFaces,',
     '        Operation::MeshVisibleFaces => AxiolidOperation::Unknown,', [TEST]),
]

probe_lib.main(MUTANTS, timeout=1800, command=(), reject_uncompiled=True)
