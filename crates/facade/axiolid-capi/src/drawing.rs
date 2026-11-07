//! Hidden-line drawings through the C ABI (#473, ADR 0105).
//!
//! One call validates every plain-data argument at the boundary (finite
//! floats, known enum values, 0/1 flags, zero reserved bytes) and refuses with
//! `InvalidArgument` before the contract sees it. Frame orthonormality, group
//! lengths and the depth range stay with the contract. A refusal that names an
//! input mesh keeps the mesh id in the context's error record, so a caller can
//! drop that mesh and retry; [`axiolid_v0_4_context_last_error_mesh`] reads it.

use axiolid::application::{ApplicationError, ApplicationErrorSource};
use axiolid::contracts::{ExecutionOptions, GeomError};
use axiolid::core::{Aabb2, Frame3, Point2, Vec3};
use axiolid::mesh_hidden_lines::{
    EdgeKind, EdgeSelection, FaceGroup, HiddenLineLimits, HiddenLineMesh, HiddenLineRequest,
    MeshId, MeshRole, View, Visibility,
};

use super::*;

/// Right-handed view frame: origin and three axes in model space.
#[repr(C)]
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct AxiolidViewFrame {
    pub origin: [f64; 3],
    pub x: [f64; 3],
    pub y: [f64; 3],
    pub z: [f64; 3],
}

/// Stable projection selector. Only orthographic is defined.
pub type AxiolidProjection = i32;
pub const AXIOLID_PROJECTION_ORTHOGRAPHIC: AxiolidProjection = 1;

/// A view: frame, depth range and projection. `near` and `far` are read only
/// when their `has_*` flag is 1.
#[repr(C)]
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct AxiolidView {
    pub frame: AxiolidViewFrame,
    pub near: f64,
    pub far: f64,
    pub projection: AxiolidProjection,
    pub has_near: u8,
    pub has_far: u8,
    /// Must be zero.
    pub reserved: [u8; 2],
}

/// Closed image-plane crop box.
#[repr(C)]
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct AxiolidImageBox {
    pub min: [f64; 2],
    pub max: [f64; 2],
}

/// Occlusion role of an input mesh.
pub type AxiolidMeshRole = i32;
/// Closed, consistently wound two-manifold.
pub const AXIOLID_MESH_ROLE_SOLID: AxiolidMeshRole = 1;
/// Edge-manifold with consistent winding, open or closed.
pub const AXIOLID_MESH_ROLE_SHEET: AxiolidMeshRole = 2;

/// How the triangles of an input mesh are grouped.
pub type AxiolidGroupMode = i32;
/// Every triangle is in group 0; `groups` is not read.
pub const AXIOLID_GROUP_MODE_SINGLE: AxiolidGroupMode = 1;
/// One id per triangle in `groups`; `group_count` is the triangle count.
pub const AXIOLID_GROUP_MODE_PER_TRIANGLE: AxiolidGroupMode = 2;
/// Planar patches derived exactly from the geometry; `groups` is not read.
pub const AXIOLID_GROUP_MODE_PLANAR: AxiolidGroupMode = 3;

/// One input mesh of a drawing. `groups` is borrowed for the call.
#[repr(C)]
#[derive(Debug, Clone, Copy)]
pub struct AxiolidViewMesh {
    pub mesh: AxiolidMeshHandle,
    /// Caller-chosen id, unique within the call, echoed in every segment.
    pub id: u64,
    pub groups: *const u32,
    pub group_count: usize,
    pub role: AxiolidMeshRole,
    pub group_mode: AxiolidGroupMode,
}

/// Which edges to draw. Flags are 0 or 1; `crease_cosine` is read only when
/// `has_crease` is 1.
#[repr(C)]
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct AxiolidEdgeSelection {
    pub crease_cosine: f64,
    pub has_crease: u8,
    pub boundary: u8,
    pub silhouette: u8,
    pub hidden: u8,
    pub group_boundary: u8,
    /// Must be zero.
    pub reserved: [u8; 3],
}

/// Hard bounds for one drawing.
#[repr(C)]
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct AxiolidHiddenLineLimits {
    pub max_meshes: usize,
    pub max_source_vertices: usize,
    pub max_source_triangles: usize,
    pub max_candidate_pairs: usize,
    pub max_output_segments: usize,
}

/// Why an edge is drawn.
#[non_exhaustive]
#[repr(i32)]
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum AxiolidEdgeKind {
    Boundary = 1,
    Silhouette = 2,
    Crease = 3,
    GroupBoundary = 4,
    Unknown = 255,
}

/// Whether a piece of an edge can be seen.
#[non_exhaustive]
#[repr(i32)]
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum AxiolidVisibility {
    Visible = 1,
    Hidden = 2,
    Unknown = 255,
}

/// One piece of a source edge in image coordinates.
#[repr(C)]
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct AxiolidHiddenLineSegment {
    pub mesh_id: u64,
    /// Position indices of the source edge, ascending.
    pub vertices: [u32; 2],
    pub start: [f64; 2],
    pub end: [f64; 2],
    pub kind: AxiolidEdgeKind,
    pub visibility: AxiolidVisibility,
}

/// Auditable counts of one drawing.
#[repr(C)]
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct AxiolidHiddenLineEvidence {
    pub source_meshes: usize,
    pub source_triangles: usize,
    pub selected_edges: usize,
    pub point_edges: usize,
    pub clipped_edges: usize,
    pub candidate_pairs: usize,
    pub exact_escalations: usize,
    pub collapsed_segments: usize,
    pub visible_segments: usize,
    pub hidden_segments: usize,
}

fn edge_kind(value: EdgeKind) -> AxiolidEdgeKind {
    match value {
        EdgeKind::Boundary => AxiolidEdgeKind::Boundary,
        EdgeKind::Silhouette => AxiolidEdgeKind::Silhouette,
        EdgeKind::Crease => AxiolidEdgeKind::Crease,
        EdgeKind::GroupBoundary => AxiolidEdgeKind::GroupBoundary,
        _ => AxiolidEdgeKind::Unknown,
    }
}

fn visibility(value: Visibility) -> AxiolidVisibility {
    match value {
        Visibility::Visible => AxiolidVisibility::Visible,
        Visibility::Hidden => AxiolidVisibility::Hidden,
        _ => AxiolidVisibility::Unknown,
    }
}

/// A boundary refusal: the status, which mesh it names and why.
struct Refusal {
    status: AxiolidStatus,
    mesh: Option<u64>,
    message: &'static str,
}

impl Refusal {
    const fn invalid(message: &'static str) -> Self {
        Self {
            status: AxiolidStatus::InvalidArgument,
            mesh: None,
            message,
        }
    }

    const fn mesh(status: AxiolidStatus, mesh: u64, message: &'static str) -> Self {
        Self {
            status,
            mesh: Some(mesh),
            message,
        }
    }
}

fn flag(value: u8, message: &'static str) -> Result<bool, Refusal> {
    match value {
        0 => Ok(false),
        1 => Ok(true),
        _ => Err(Refusal::invalid(message)),
    }
}

fn finite(values: &[f64], message: &'static str) -> Result<(), Refusal> {
    if values.iter().all(|value| value.is_finite()) {
        Ok(())
    } else {
        Err(Refusal::invalid(message))
    }
}

fn vec3(value: [f64; 3]) -> Vec3 {
    Vec3::new(value[0], value[1], value[2])
}

fn read_view(view: &AxiolidView) -> Result<View, Refusal> {
    if view.projection != AXIOLID_PROJECTION_ORTHOGRAPHIC {
        return Err(Refusal::invalid("unknown view projection"));
    }
    if view.reserved != [0; 2] {
        return Err(Refusal::invalid("view reserved bytes must be zero"));
    }
    let frame = &view.frame;
    for axis in [&frame.origin, &frame.x, &frame.y, &frame.z] {
        finite(axis, "view frame must be finite")?;
    }
    finite(&[view.near, view.far], "view depth range must be finite")?;
    let has_near = flag(view.has_near, "has_near must be 0 or 1")?;
    let has_far = flag(view.has_far, "has_far must be 0 or 1")?;
    let mut result = View::new(Frame3 {
        origin: vec3(frame.origin),
        x: vec3(frame.x),
        y: vec3(frame.y),
        z: vec3(frame.z),
    });
    if has_near {
        result = result.with_near(view.near);
    }
    if has_far {
        result = result.with_far(view.far);
    }
    Ok(result)
}

fn read_edges(edges: &AxiolidEdgeSelection) -> Result<EdgeSelection, Refusal> {
    if edges.reserved != [0; 3] {
        return Err(Refusal::invalid(
            "edge selection reserved bytes must be zero",
        ));
    }
    finite(&[edges.crease_cosine], "crease cosine must be finite")?;
    let has_crease = flag(edges.has_crease, "has_crease must be 0 or 1")?;
    Ok(
        EdgeSelection::new(has_crease.then_some(edges.crease_cosine))
            .with_boundary(flag(edges.boundary, "boundary must be 0 or 1")?)
            .with_silhouette(flag(edges.silhouette, "silhouette must be 0 or 1")?)
            .with_hidden(flag(edges.hidden, "hidden must be 0 or 1")?)
            .with_group_boundary(flag(edges.group_boundary, "group_boundary must be 0 or 1")?),
    )
}

fn read_region(region: &AxiolidImageBox) -> Result<Aabb2, Refusal> {
    finite(&[region.min, region.max].concat(), "region must be finite")?;
    Ok(Aabb2 {
        min: Point2::new(region.min[0], region.min[1]),
        max: Point2::new(region.max[0], region.max[1]),
    })
}

/// One decoded input mesh before it borrows the context's mesh.
struct DecodedMesh {
    id: u64,
    handle: u64,
    role: MeshRole,
    mode: AxiolidGroupMode,
    groups: Vec<FaceGroup>,
}

fn read_mesh(entry: &AxiolidViewMesh, context: &Context) -> Result<DecodedMesh, Refusal> {
    let role = match entry.role {
        AXIOLID_MESH_ROLE_SOLID => MeshRole::Solid,
        AXIOLID_MESH_ROLE_SHEET => MeshRole::Sheet,
        _ => {
            return Err(Refusal::mesh(
                AxiolidStatus::InvalidArgument,
                entry.id,
                "unknown mesh role",
            ))
        }
    };
    let Some(mesh) = context.meshes.get(entry.mesh.0) else {
        return Err(Refusal::mesh(
            AxiolidStatus::InvalidHandle,
            entry.id,
            "unknown mesh handle",
        ));
    };
    let groups = match entry.group_mode {
        AXIOLID_GROUP_MODE_SINGLE | AXIOLID_GROUP_MODE_PLANAR => {
            if entry.group_count != 0 {
                return Err(Refusal::mesh(
                    AxiolidStatus::InvalidArgument,
                    entry.id,
                    "group_count must be zero unless the group mode is per-triangle",
                ));
            }
            Vec::new()
        }
        AXIOLID_GROUP_MODE_PER_TRIANGLE => {
            // The count is checked against the mesh before any caller memory is read.
            if entry.group_count != mesh.indices.len() / 3 {
                return Err(Refusal::mesh(
                    AxiolidStatus::InvalidArgument,
                    entry.id,
                    "group_count must equal the mesh's triangle count",
                ));
            }
            if entry.group_count == 0 {
                Vec::new()
            } else if entry.groups.is_null() {
                return Err(Refusal::mesh(
                    AxiolidStatus::NullPointer,
                    entry.id,
                    "per-triangle groups are null",
                ));
            } else {
                // SAFETY: non-null and non-zero length were checked; the caller promises a
                // readable array of exactly `group_count` ids for this call.
                let ids = unsafe { std::slice::from_raw_parts(entry.groups, entry.group_count) };
                ids.iter().copied().map(FaceGroup).collect()
            }
        }
        _ => {
            return Err(Refusal::mesh(
                AxiolidStatus::InvalidArgument,
                entry.id,
                "unknown group mode",
            ))
        }
    };
    Ok(DecodedMesh {
        id: entry.id,
        handle: entry.mesh.0,
        role,
        mode: entry.group_mode,
        groups,
    })
}

fn record_refusal(
    mut context: MutexGuard<'_, Context>,
    refusal: Refusal,
    tolerance: AxiolidTolerance,
) -> AxiolidStatus {
    context.last_error = Some(ErrorRecord {
        status: refusal.status,
        operation: AxiolidOperation::MeshHiddenLines,
        tolerance,
        provider: "axiolid-capi".into(),
        message: refusal.message.into(),
        mesh: refusal.mesh,
    });
    refusal.status
}

/// Map a refused drawing to a status and the mesh it names.
///
/// A capability the provider lacks, or an input family it does not support, is
/// `Unsupported` and an exhausted budget
/// `LimitExceeded`; anything else, including a provider that failed its
/// conformance suite, is `OperationFailed`.
fn record_drawing_error(
    mut context: MutexGuard<'_, Context>,
    error: ApplicationError,
) -> AxiolidStatus {
    let (geometry, mesh) = match &*error.source {
        ApplicationErrorSource::HiddenLines(refusal) => {
            (Some(&refusal.error), refusal.mesh.map(|MeshId(id)| id))
        }
        ApplicationErrorSource::Geometry(geometry) => (Some(geometry), None),
        _ => (None, None),
    };
    let status = match geometry {
        Some(
            GeomError::Unsupported { .. }
            | GeomError::UnsupportedInput { .. }
            | GeomError::Unavailable { .. },
        ) => AxiolidStatus::Unsupported,
        Some(GeomError::BudgetExceeded { .. }) => AxiolidStatus::LimitExceeded,
        _ => AxiolidStatus::OperationFailed,
    };
    context.last_error = Some(ErrorRecord {
        status,
        operation: AxiolidOperation::MeshHiddenLines,
        tolerance: AxiolidTolerance {
            linear: error.context.tolerance.linear(),
            angular: error.context.tolerance.angular(),
        },
        provider: error.context.provider.as_str().into(),
        message: error.to_string(),
        mesh,
    });
    status
}

/// Draw the edges of borrowed meshes with hidden-line removal and return an
/// owned result handle (read it with the `hidden_lines_*` functions).
///
/// `meshes` and every `groups` array are borrowed for this call. A null
/// `region` is no crop. The context's `max_meshes` bounds `mesh_count`.
/// # Safety
/// Every non-null pointer must be aligned and valid for the documented read or write extent.
#[no_mangle]
pub unsafe extern "C" fn axiolid_v0_4_hidden_lines(
    context: AxiolidContextHandle,
    meshes: *const AxiolidViewMesh,
    mesh_count: usize,
    view: *const AxiolidView,
    edges: *const AxiolidEdgeSelection,
    region: *const AxiolidImageBox,
    limits: *const AxiolidHiddenLineLimits,
    tolerance_value: AxiolidTolerance,
    out_result: *mut AxiolidResultHandle,
) -> AxiolidStatus {
    boundary(|| {
        if out_result.is_null()
            || view.is_null()
            || edges.is_null()
            || limits.is_null()
            || (mesh_count != 0 && meshes.is_null())
        {
            return AxiolidStatus::NullPointer;
        }
        let Some(context) = context_entry(context.0) else {
            return AxiolidStatus::InvalidHandle;
        };
        let mut context = lock_unpoisoned(&context);
        if mesh_count > context.config.max_meshes as usize {
            return AxiolidStatus::LimitExceeded;
        }
        // SAFETY: null was rejected; the caller promises readable aligned POD storage.
        let (view, edges, limits) = unsafe { (view.read(), edges.read(), limits.read()) };
        // SAFETY: a null region is no crop; otherwise the caller promises readable POD storage.
        let region = unsafe { region.as_ref() };
        let entries: &[AxiolidViewMesh] = if mesh_count == 0 {
            &[]
        } else {
            let Some(bytes) = mesh_count.checked_mul(std::mem::size_of::<AxiolidViewMesh>()) else {
                return AxiolidStatus::LimitExceeded;
            };
            if bytes > isize::MAX as usize {
                return AxiolidStatus::LimitExceeded;
            }
            // SAFETY: non-zero pointer/length was validated; the caller promises a readable
            // array of `mesh_count` entries for this call.
            unsafe { std::slice::from_raw_parts(meshes, mesh_count) }
        };
        let Ok(tolerance) = tolerance(tolerance_value) else {
            return record_refusal(
                context,
                Refusal::invalid("tolerance must be finite and non-negative"),
                tolerance_value,
            );
        };
        let decoded = (|| {
            let view = read_view(&view)?;
            let edges = read_edges(&edges)?;
            let region = region.map(read_region).transpose()?;
            let meshes = entries
                .iter()
                .map(|entry| read_mesh(entry, &context))
                .collect::<Result<Vec<_>, _>>()?;
            Ok((view, edges, region, meshes))
        })();
        let (view, edges, region, decoded) = match decoded {
            Ok(decoded) => decoded,
            Err(refusal) => return record_refusal(context, refusal, tolerance_value),
        };
        let request_limits = HiddenLineLimits::new(
            limits.max_meshes,
            limits.max_source_vertices,
            limits.max_source_triangles,
            limits.max_candidate_pairs,
            limits.max_output_segments,
        );
        let options = ExecutionOptions::new(tolerance);
        // The request borrows the context's meshes; it ends with this block so the
        // guard can be moved into the error record and the result table afterwards.
        let outcome = {
            let state: &Context = &context;
            let mut inputs = Vec::with_capacity(decoded.len());
            for entry in &decoded {
                let Some(mesh) = state.meshes.get(entry.handle) else {
                    return record_refusal(
                        context,
                        Refusal::mesh(
                            AxiolidStatus::InvalidHandle,
                            entry.id,
                            "unknown mesh handle",
                        ),
                        tolerance_value,
                    );
                };
                let input = HiddenLineMesh::new(MeshId(entry.id), entry.role, mesh);
                inputs.push(match entry.mode {
                    AXIOLID_GROUP_MODE_PER_TRIANGLE => input.with_groups(&entry.groups),
                    AXIOLID_GROUP_MODE_PLANAR => input.with_planar_groups(),
                    _ => input,
                });
            }
            let mut request = HiddenLineRequest::new(&inputs, view, edges, request_limits);
            if let Some(region) = region {
                request = request.with_region(region);
            }
            state.application.hidden_lines(&request, &options)
        };
        let outcome = match outcome {
            Ok(outcome) => outcome,
            Err(error) => return record_drawing_error(context, error),
        };
        let Ok(handle) = insert_result(&mut context, StoredResult::HiddenLines(outcome)) else {
            return AxiolidStatus::LimitExceeded;
        };
        // SAFETY: null was rejected; caller provides one writable result handle.
        unsafe { out_result.write(AxiolidResultHandle(handle)) };
        AxiolidStatus::Ok
    })
}

/// Run `read` on the hidden-line outcome behind `result`.
fn with_outcome(
    context: AxiolidContextHandle,
    result: AxiolidResultHandle,
    read: impl FnOnce(&axiolid::mesh_hidden_lines::HiddenLineOutcome) -> AxiolidStatus,
) -> AxiolidStatus {
    let Some(context) = context_entry(context.0) else {
        return AxiolidStatus::InvalidHandle;
    };
    let context = lock_unpoisoned(&context);
    match context.results.get(result.0) {
        Some(StoredResult::HiddenLines(outcome)) => read(outcome),
        Some(_) => AxiolidStatus::WrongResultKind,
        None => AxiolidStatus::InvalidHandle,
    }
}

/// Query the number of segments of a hidden-line result.
/// # Safety
/// Every non-null pointer must be aligned and valid for the documented read or write extent.
#[no_mangle]
pub unsafe extern "C" fn axiolid_v0_4_hidden_lines_count(
    context: AxiolidContextHandle,
    result: AxiolidResultHandle,
    out_count: *mut usize,
) -> AxiolidStatus {
    boundary(|| {
        if out_count.is_null() {
            return AxiolidStatus::NullPointer;
        }
        with_outcome(context, result, |outcome| {
            // SAFETY: null was rejected; caller provides one writable count.
            unsafe { out_count.write(outcome.segments.len()) };
            AxiolidStatus::Ok
        })
    })
}

/// Copy the segments of a hidden-line result, in the outcome's order, into
/// caller-owned storage.
///
/// A `capacity` below the count returns `BufferTooSmall` and writes nothing.
/// # Safety
/// Every non-null pointer must be aligned and valid for the documented read or write extent.
#[no_mangle]
pub unsafe extern "C" fn axiolid_v0_4_hidden_lines_copy(
    context: AxiolidContextHandle,
    result: AxiolidResultHandle,
    out: *mut AxiolidHiddenLineSegment,
    capacity: usize,
) -> AxiolidStatus {
    boundary(|| {
        if capacity != 0 && out.is_null() {
            return AxiolidStatus::NullPointer;
        }
        with_outcome(context, result, |outcome| {
            if capacity < outcome.segments.len() {
                return AxiolidStatus::BufferTooSmall;
            }
            for (index, segment) in outcome.segments.iter().enumerate() {
                let value = AxiolidHiddenLineSegment {
                    mesh_id: segment.source.mesh.0,
                    vertices: segment.source.vertices,
                    start: [segment.start.x, segment.start.y],
                    end: [segment.end.x, segment.end.y],
                    kind: edge_kind(segment.kind),
                    visibility: visibility(segment.visibility),
                };
                // SAFETY: index < len <= capacity; the caller promises `capacity` writable slots.
                unsafe { out.add(index).write(value) };
            }
            AxiolidStatus::Ok
        })
    })
}

/// Copy the evidence of a hidden-line result.
/// # Safety
/// Every non-null pointer must be aligned and valid for the documented read or write extent.
#[no_mangle]
pub unsafe extern "C" fn axiolid_v0_4_hidden_lines_evidence(
    context: AxiolidContextHandle,
    result: AxiolidResultHandle,
    out: *mut AxiolidHiddenLineEvidence,
) -> AxiolidStatus {
    boundary(|| {
        if out.is_null() {
            return AxiolidStatus::NullPointer;
        }
        with_outcome(context, result, |outcome| {
            let evidence = &outcome.evidence;
            let value = AxiolidHiddenLineEvidence {
                source_meshes: evidence.source_meshes,
                source_triangles: evidence.source_triangles,
                selected_edges: evidence.selected_edges,
                point_edges: evidence.point_edges,
                clipped_edges: evidence.clipped_edges,
                candidate_pairs: evidence.candidate_pairs,
                exact_escalations: evidence.exact_escalations,
                collapsed_segments: evidence.collapsed_segments,
                visible_segments: evidence.visible_segments,
                hidden_segments: evidence.hidden_segments,
            };
            // SAFETY: null was rejected; caller provides one writable evidence record.
            unsafe { out.write(value) };
            AxiolidStatus::Ok
        })
    })
}

/// The input mesh the context's most recent error names.
///
/// `out_has_mesh` is 1 and `out_mesh_id` the caller's id when the error names
/// a mesh, otherwise 0 and 0. Returns `NoError` when no error was recorded.
/// # Safety
/// Every non-null pointer must be aligned and valid for the documented read or write extent.
#[no_mangle]
pub unsafe extern "C" fn axiolid_v0_4_context_last_error_mesh(
    context: AxiolidContextHandle,
    out_has_mesh: *mut u8,
    out_mesh_id: *mut u64,
) -> AxiolidStatus {
    boundary(|| {
        if out_has_mesh.is_null() || out_mesh_id.is_null() {
            return AxiolidStatus::NullPointer;
        }
        let Some(context) = context_entry(context.0) else {
            return AxiolidStatus::InvalidHandle;
        };
        let context = lock_unpoisoned(&context);
        let Some(error) = &context.last_error else {
            return AxiolidStatus::NoError;
        };
        // SAFETY: null was rejected; caller provides one writable value each.
        unsafe {
            out_has_mesh.write(u8::from(error.mesh.is_some()));
            out_mesh_id.write(error.mesh.unwrap_or(0));
        }
        AxiolidStatus::Ok
    })
}
