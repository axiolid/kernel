use super::*;

/// Import an indexed triangle mesh. Input buffers are borrowed for this call only.
/// # Safety
/// Every non-null pointer must be aligned and valid for the documented read or write extent.
#[no_mangle]
pub unsafe extern "C" fn axiolid_v0_4_mesh_import(
    context: AxiolidContextHandle,
    positions_xyz: *const f64,
    vertex_count: usize,
    triangle_indices: *const u32,
    triangle_count: usize,
    out_mesh: *mut AxiolidMeshHandle,
) -> AxiolidStatus {
    boundary(|| {
        if out_mesh.is_null()
            || (vertex_count != 0 && positions_xyz.is_null())
            || (triangle_count != 0 && triangle_indices.is_null())
        {
            return AxiolidStatus::NullPointer;
        }
        let Some(position_len) = vertex_count.checked_mul(3) else {
            return AxiolidStatus::LimitExceeded;
        };
        let Some(index_len) = triangle_count.checked_mul(3) else {
            return AxiolidStatus::LimitExceeded;
        };
        let Some(context) = context_entry(context.0) else {
            return AxiolidStatus::InvalidHandle;
        };
        let mut context = lock_unpoisoned(&context);
        if vertex_count > context.config.max_vertices_per_mesh as usize
            || triangle_count > context.config.max_triangles_per_mesh as usize
            || context.meshes.live_count() >= context.config.max_meshes as usize
        {
            return AxiolidStatus::LimitExceeded;
        }
        let positions: &[f64] = if position_len == 0 {
            &[]
        } else {
            // SAFETY: non-zero pointer/length combinations and overflow were rejected; the
            // caller promises a readable array of exactly the declared length for this call.
            unsafe { std::slice::from_raw_parts(positions_xyz, position_len) }
        };
        let indices: &[u32] = if index_len == 0 {
            &[]
        } else {
            // SAFETY: same invariant as positions, for three indices per triangle.
            unsafe { std::slice::from_raw_parts(triangle_indices, index_len) }
        };
        if positions.iter().any(|value| !value.is_finite()) {
            return record_error(
                context,
                AxiolidStatus::InvalidArgument,
                AxiolidOperation::Healing,
                AxiolidTolerance {
                    linear: 0.0,
                    angular: 0.0,
                },
                "axiolid-capi",
                "mesh positions must all be finite",
            );
        }
        if indices
            .iter()
            .any(|index| usize::try_from(*index).map_or(true, |index| index >= vertex_count))
        {
            return record_error(
                context,
                AxiolidStatus::InvalidArgument,
                AxiolidOperation::Healing,
                AxiolidTolerance {
                    linear: 0.0,
                    angular: 0.0,
                },
                "axiolid-capi",
                "triangle index is outside the imported vertex array",
            );
        }
        let positions = positions
            .as_chunks::<3>()
            .0
            .iter()
            .map(|value| Point3::new(value[0], value[1], value[2]))
            .collect();
        let handle = context
            .meshes
            .insert(TriMesh::new(positions, indices.to_vec()));
        // SAFETY: null was rejected; caller provides one writable handle.
        unsafe { out_mesh.write(AxiolidMeshHandle(handle)) };
        AxiolidStatus::Ok
    })
}

/// Query vertex and triangle counts before allocating export buffers.
/// # Safety
/// Every non-null pointer must be aligned and valid for the documented read or write extent.
#[no_mangle]
pub unsafe extern "C" fn axiolid_v0_4_mesh_counts(
    context: AxiolidContextHandle,
    mesh: AxiolidMeshHandle,
    out_vertex_count: *mut usize,
    out_triangle_count: *mut usize,
) -> AxiolidStatus {
    boundary(|| {
        if out_vertex_count.is_null() || out_triangle_count.is_null() {
            return AxiolidStatus::NullPointer;
        }
        let Some(context) = context_entry(context.0) else {
            return AxiolidStatus::InvalidHandle;
        };
        let context = lock_unpoisoned(&context);
        let Some(mesh) = context.meshes.get(mesh.0) else {
            return AxiolidStatus::InvalidHandle;
        };
        // SAFETY: null was rejected; caller provides one writable count each.
        unsafe {
            out_vertex_count.write(mesh.positions.len());
            out_triangle_count.write(mesh.indices.len() / 3);
        }
        AxiolidStatus::Ok
    })
}

/// Copy an indexed mesh into exactly-sized, caller-owned output buffers.
/// # Safety
/// Every non-null pointer must be aligned and valid for the documented write extent. The
/// position and index output ranges must not overlap.
#[no_mangle]
pub unsafe extern "C" fn axiolid_v0_4_mesh_copy(
    context: AxiolidContextHandle,
    mesh: AxiolidMeshHandle,
    out_positions_xyz: *mut f64,
    vertex_capacity: usize,
    out_triangle_indices: *mut u32,
    triangle_capacity: usize,
) -> AxiolidStatus {
    boundary(|| {
        if (vertex_capacity != 0 && out_positions_xyz.is_null())
            || (triangle_capacity != 0 && out_triangle_indices.is_null())
        {
            return AxiolidStatus::NullPointer;
        }
        let Some(context) = context_entry(context.0) else {
            return AxiolidStatus::InvalidHandle;
        };
        let context = lock_unpoisoned(&context);
        let Some(mesh) = context.meshes.get(mesh.0) else {
            return AxiolidStatus::InvalidHandle;
        };
        let triangle_count = mesh.indices.len() / 3;
        if vertex_capacity < mesh.positions.len() || triangle_capacity < triangle_count {
            return AxiolidStatus::BufferTooSmall;
        }
        let Some(position_len) = mesh.positions.len().checked_mul(3) else {
            return AxiolidStatus::LimitExceeded;
        };
        let positions: &mut [f64] = if position_len == 0 {
            &mut []
        } else {
            // SAFETY: capacities were checked against exact output lengths; caller
            // promises a writable, non-overlapping buffer for the duration of the call.
            unsafe { std::slice::from_raw_parts_mut(out_positions_xyz, position_len) }
        };
        let indices: &mut [u32] = if mesh.indices.is_empty() {
            &mut []
        } else {
            // SAFETY: same invariant as positions, with three indices per triangle.
            unsafe { std::slice::from_raw_parts_mut(out_triangle_indices, mesh.indices.len()) }
        };
        for (target, point) in positions
            .as_chunks_mut::<3>()
            .0
            .iter_mut()
            .zip(&mesh.positions)
        {
            target.copy_from_slice(&[point.x, point.y, point.z]);
        }
        indices.copy_from_slice(&mesh.indices);
        AxiolidStatus::Ok
    })
}

/// Destroy a mesh handle. Repeated destruction returns `InvalidHandle`.
#[no_mangle]
pub extern "C" fn axiolid_v0_4_mesh_destroy(
    context: AxiolidContextHandle,
    mesh: AxiolidMeshHandle,
) -> AxiolidStatus {
    boundary(|| {
        let Some(context) = context_entry(context.0) else {
            return AxiolidStatus::InvalidHandle;
        };
        let mut context = lock_unpoisoned(&context);
        match context.meshes.remove(mesh.0) {
            Some(_) => AxiolidStatus::Ok,
            None => AxiolidStatus::InvalidHandle,
        }
    })
}

/// Caller-supplied tolerance in model units and radians.
#[repr(C)]
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct AxiolidTolerance {
    pub linear: f64,
    pub angular: f64,
}

/// Column-major affine transform.
#[repr(C)]
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct AxiolidTransform {
    pub columns: [f64; 16],
}

/// Axis-aligned mesh bounds.
#[repr(C)]
#[derive(Debug, Clone, Copy, PartialEq, Default)]
pub struct AxiolidBounds {
    pub min: [f64; 3],
    pub max: [f64; 3],
}

/// Stable mesh-audit summary.
#[repr(C)]
#[derive(Debug, Clone, Copy, PartialEq, Eq, Default)]
pub struct AxiolidMeshAudit {
    pub vertex_count: usize,
    pub triangle_count: usize,
    pub closed_two_manifold: u8,
    /// Must be zero; reserves layout space without exposing padding bytes.
    pub reserved: [u8; 7],
}

/// Surface and signed-volume measurements.
#[repr(C)]
#[derive(Debug, Clone, Copy, PartialEq, Default)]
pub struct AxiolidMeasurements {
    pub surface_area: f64,
    pub signed_volume: f64,
    pub surface_centroid: [f64; 3],
    pub volume_centroid: [f64; 3],
}

/// Audit mesh topology under an explicit tolerance.
/// # Safety
/// Every non-null pointer must be aligned and valid for the documented read or write extent.
#[no_mangle]
pub unsafe extern "C" fn axiolid_v0_4_mesh_audit(
    context: AxiolidContextHandle,
    mesh: AxiolidMeshHandle,
    tolerance_value: AxiolidTolerance,
    out_audit: *mut AxiolidMeshAudit,
) -> AxiolidStatus {
    boundary(|| {
        if out_audit.is_null() {
            return AxiolidStatus::NullPointer;
        }
        let Ok(tolerance) = tolerance(tolerance_value) else {
            return AxiolidStatus::InvalidArgument;
        };
        let Some(context) = context_entry(context.0) else {
            return AxiolidStatus::InvalidHandle;
        };
        let context = lock_unpoisoned(&context);
        let Some(mesh) = context.meshes.get(mesh.0) else {
            return AxiolidStatus::InvalidHandle;
        };
        let mesh = mesh.clone();
        let health = match context.application.validate_mesh(&mesh, tolerance) {
            Ok(health) => health,
            Err(error) => return record_application_error(context, error),
        };
        let result = AxiolidMeshAudit {
            vertex_count: mesh.positions.len(),
            triangle_count: mesh.indices.len() / 3,
            closed_two_manifold: u8::from(health.is_closed_two_manifold()),
            reserved: [0; 7],
        };
        // SAFETY: null was rejected; caller provides one writable audit record.
        unsafe { out_audit.write(result) };
        AxiolidStatus::Ok
    })
}

/// Compute axis-aligned bounds for a non-empty mesh.
/// # Safety
/// Every non-null pointer must be aligned and valid for the documented read or write extent.
#[no_mangle]
pub unsafe extern "C" fn axiolid_v0_4_mesh_bounds(
    context: AxiolidContextHandle,
    mesh: AxiolidMeshHandle,
    out_bounds: *mut AxiolidBounds,
) -> AxiolidStatus {
    boundary(|| {
        if out_bounds.is_null() {
            return AxiolidStatus::NullPointer;
        }
        let Some(context) = context_entry(context.0) else {
            return AxiolidStatus::InvalidHandle;
        };
        let context = lock_unpoisoned(&context);
        let Some(mesh) = context.meshes.get(mesh.0) else {
            return AxiolidStatus::InvalidHandle;
        };
        let Some(first) = mesh.positions.first() else {
            return AxiolidStatus::InvalidArgument;
        };
        let mut min = [first.x, first.y, first.z];
        let mut max = min;
        for point in &mesh.positions[1..] {
            min[0] = min[0].min(point.x);
            min[1] = min[1].min(point.y);
            min[2] = min[2].min(point.z);
            max[0] = max[0].max(point.x);
            max[1] = max[1].max(point.y);
            max[2] = max[2].max(point.z);
        }
        // SAFETY: null was rejected; caller provides one writable bounds record.
        unsafe { out_bounds.write(AxiolidBounds { min, max }) };
        AxiolidStatus::Ok
    })
}

/// Measure surface area and signed volume under an explicit tolerance.
/// # Safety
/// Every non-null pointer must be aligned and valid for the documented read or write extent.
#[no_mangle]
pub unsafe extern "C" fn axiolid_v0_4_mesh_measure(
    context: AxiolidContextHandle,
    mesh: AxiolidMeshHandle,
    tolerance_value: AxiolidTolerance,
    out_measurements: *mut AxiolidMeasurements,
) -> AxiolidStatus {
    boundary(|| {
        if out_measurements.is_null() {
            return AxiolidStatus::NullPointer;
        }
        let Ok(tolerance) = tolerance(tolerance_value) else {
            return AxiolidStatus::InvalidArgument;
        };
        let Some(context) = context_entry(context.0) else {
            return AxiolidStatus::InvalidHandle;
        };
        let context = lock_unpoisoned(&context);
        let Some(mesh) = context.meshes.get(mesh.0) else {
            return AxiolidStatus::InvalidHandle;
        };
        let mesh = mesh.clone();
        let measurements = match context.application.measure_mesh(&mesh, tolerance) {
            Ok(measurements) => measurements,
            Err(error) => return record_application_error(context, error),
        };
        let result = AxiolidMeasurements {
            surface_area: measurements.surface.area,
            signed_volume: measurements.volume.signed_volume,
            surface_centroid: [
                measurements.surface.centroid.x,
                measurements.surface.centroid.y,
                measurements.surface.centroid.z,
            ],
            volume_centroid: [
                measurements.volume.centroid.x,
                measurements.volume.centroid.y,
                measurements.volume.centroid.z,
            ],
        };
        // SAFETY: null was rejected; caller provides one writable measurement record.
        unsafe { out_measurements.write(result) };
        AxiolidStatus::Ok
    })
}

/// Apply a finite affine transform to mesh positions in place.
/// # Safety
/// Every non-null pointer must be aligned and valid for the documented read or write extent.
#[no_mangle]
pub unsafe extern "C" fn axiolid_v0_4_mesh_transform(
    context: AxiolidContextHandle,
    mesh: AxiolidMeshHandle,
    transform: *const AxiolidTransform,
) -> AxiolidStatus {
    boundary(|| {
        if transform.is_null() {
            return AxiolidStatus::NullPointer;
        }
        // SAFETY: null was rejected; caller provides one readable POD transform.
        let transform = unsafe { transform.read() };
        let m = transform.columns;
        if m.iter().any(|value| !value.is_finite())
            || m[3] != 0.0
            || m[7] != 0.0
            || m[11] != 0.0
            || m[15] != 1.0
        {
            return AxiolidStatus::InvalidArgument;
        }
        let Some(context) = context_entry(context.0) else {
            return AxiolidStatus::InvalidHandle;
        };
        let mut context = lock_unpoisoned(&context);
        let Some(source) = context.meshes.get(mesh.0) else {
            return AxiolidStatus::InvalidHandle;
        };
        let mut transformed = Vec::with_capacity(source.positions.len());
        for point in &source.positions {
            let [x, y, z] = [point.x, point.y, point.z];
            let point = Point3::new(
                m[0].mul_add(x, m[4].mul_add(y, m[8].mul_add(z, m[12]))),
                m[1].mul_add(x, m[5].mul_add(y, m[9].mul_add(z, m[13]))),
                m[2].mul_add(x, m[6].mul_add(y, m[10].mul_add(z, m[14]))),
            );
            if !point.is_finite() {
                return record_error(
                    context,
                    AxiolidStatus::InvalidArgument,
                    AxiolidOperation::Healing,
                    AxiolidTolerance {
                        linear: 0.0,
                        angular: 0.0,
                    },
                    "axiolid-capi",
                    "transform produced non-finite coordinates",
                );
            }
            transformed.push(point);
        }
        let Some(mesh) = context.meshes.get_mut(mesh.0) else {
            return AxiolidStatus::InvalidHandle;
        };
        mesh.positions = transformed;
        AxiolidStatus::Ok
    })
}

/// Globally unique opaque operation-result token owned by a context.
#[repr(transparent)]
#[derive(Debug, Clone, Copy, PartialEq, Eq, Default)]
pub struct AxiolidResultHandle(pub u64);

impl AxiolidResultHandle {
    pub const INVALID: Self = Self(0);
}

/// Representation stored by an operation result.
#[repr(i32)]
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum AxiolidGeometryKind {
    TriangleMesh = 1,
    ExactBrep = 2,
    Unknown = 255,
}
