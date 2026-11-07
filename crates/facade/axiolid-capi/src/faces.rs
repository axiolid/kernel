//! Visible face regions through the C ABI (#386, ADR 0094, ADR 0105).
//!
//! The call takes the meshes, view, crop and tolerance of
//! [`axiolid_v0_4_hidden_lines`] and validates them the same way. The result
//! is a flat, offset-indexed table: regions index into rings, rings into
//! points, and each region owns a run of triangle indices. A region's first
//! ring is its outer ring (counter-clockwise); the rest are holes.
//! `visible_faces_copy` writes nothing unless every capacity covers its
//! count. A refusal that names a mesh keeps its id in the same error record
//! as hidden lines.

use axiolid::contracts::ExecutionOptions;
use axiolid::mesh_hidden_lines::{
    FaceSource, VisibleFaceLimits, VisibleFaceOutcome, VisibleFaceRequest,
};

use super::*;

/// Hard bounds for one visible-face call.
#[repr(C)]
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct AxiolidVisibleFaceLimits {
    pub max_meshes: usize,
    pub max_source_vertices: usize,
    pub max_source_triangles: usize,
    pub max_candidate_pairs: usize,
    pub max_regions: usize,
    pub max_region_vertices: usize,
}

/// What a visible region is a region of.
#[non_exhaustive]
#[repr(i32)]
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum AxiolidFaceSource {
    /// The triangles of one face group; `AxiolidVisibleRegion::group` names it.
    Group = 1,
    /// The cap of a solid cut by the near plane; `group` is 0.
    Cap = 2,
    Unknown = 255,
}

/// One visible region: a run of rings and a run of triangle indices.
#[repr(C)]
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct AxiolidVisibleRegion {
    /// The caller's mesh id.
    pub mesh_id: u64,
    /// Index of the region's outer ring in the ring table; its holes follow.
    pub first_ring: usize,
    /// Rings of the region: the outer ring plus its holes.
    pub ring_count: usize,
    /// Index of the region's first triangle in the triangle table.
    pub first_triangle: usize,
    /// Frontmost triangles of the region; zero for a cap.
    pub triangle_count: usize,
    /// The face group when `source` is `Group`, otherwise 0.
    pub group: u32,
    pub source: AxiolidFaceSource,
}

/// One ring: a run of points in the point table.
#[repr(C)]
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct AxiolidRegionRing {
    pub first_point: usize,
    pub point_count: usize,
}

/// Table sizes of a visible-face result.
#[repr(C)]
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct AxiolidVisibleFaceCounts {
    pub regions: usize,
    pub rings: usize,
    pub points: usize,
    pub triangles: usize,
}

/// Auditable counts of one visible-face call.
#[repr(C)]
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct AxiolidVisibleFaceEvidence {
    pub source_meshes: usize,
    pub source_triangles: usize,
    pub edge_on_triangles: usize,
    pub removed_triangles: usize,
    pub boundary_segments: usize,
    pub intersection_segments: usize,
    pub crossings: usize,
    pub pieces: usize,
    pub tied_sides: usize,
    pub candidate_pairs: usize,
    pub exact_escalations: usize,
    pub regions: usize,
    pub cap_regions: usize,
    pub rerouted_edges: usize,
    pub collapsed_regions: usize,
    pub collapsed_holes: usize,
    /// Bound, in the maximum norm of the image, on how far snap rounding
    /// moved region boundaries.
    pub snap_displacement: f64,
}

fn face_source(value: FaceSource) -> (AxiolidFaceSource, u32) {
    match value {
        FaceSource::Group(group) => (AxiolidFaceSource::Group, group.0),
        FaceSource::Cap => (AxiolidFaceSource::Cap, 0),
        _ => (AxiolidFaceSource::Unknown, 0),
    }
}

fn counts(outcome: &VisibleFaceOutcome) -> AxiolidVisibleFaceCounts {
    let mut counts = AxiolidVisibleFaceCounts {
        regions: outcome.regions.len(),
        rings: 0,
        points: 0,
        triangles: 0,
    };
    for region in &outcome.regions {
        counts.rings += 1 + region.holes.len();
        counts.points += region.outer.len() + region.holes.iter().map(Vec::len).sum::<usize>();
        counts.triangles += region.triangles.len();
    }
    counts
}

/// Find the regions where each face group is frontmost and return an owned
/// result handle (read it with the `visible_faces_*` functions).
///
/// `meshes` and every `groups` array are borrowed for this call. A null
/// `region` is no crop. The context's `max_meshes` bounds `mesh_count`.
/// # Safety
/// Every non-null pointer must be aligned and valid for the documented read or write extent.
#[no_mangle]
pub unsafe extern "C" fn axiolid_v0_4_visible_faces(
    context: AxiolidContextHandle,
    meshes: *const AxiolidViewMesh,
    mesh_count: usize,
    view: *const AxiolidView,
    region: *const AxiolidImageBox,
    limits: *const AxiolidVisibleFaceLimits,
    tolerance_value: AxiolidTolerance,
    out_result: *mut AxiolidResultHandle,
) -> AxiolidStatus {
    const OPERATION: AxiolidOperation = AxiolidOperation::MeshVisibleFaces;
    boundary(|| {
        if out_result.is_null()
            || view.is_null()
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
        let (view, limits) = unsafe { (view.read(), limits.read()) };
        // SAFETY: a null region is no crop; otherwise the caller promises readable POD storage.
        let region = unsafe { region.as_ref() };
        // SAFETY: a non-null pointer was validated for a non-zero count; the caller promises a
        // readable array of `mesh_count` entries for this call.
        let entries = match unsafe { mesh_entries(meshes, mesh_count) } {
            Ok(entries) => entries,
            Err(status) => return status,
        };
        let Ok(tolerance) = tolerance(tolerance_value) else {
            return record_refusal(
                context,
                OPERATION,
                Refusal::invalid("tolerance must be finite and non-negative"),
                tolerance_value,
            );
        };
        let decoded = (|| {
            let view = read_view(&view)?;
            let region = region.map(read_region).transpose()?;
            let meshes = entries
                .iter()
                .map(|entry| read_mesh(entry, &context))
                .collect::<Result<Vec<_>, _>>()?;
            Ok((view, region, meshes))
        })();
        let (view, region, decoded) = match decoded {
            Ok(decoded) => decoded,
            Err(refusal) => return record_refusal(context, OPERATION, refusal, tolerance_value),
        };
        let request_limits = VisibleFaceLimits::new(
            limits.max_meshes,
            limits.max_source_vertices,
            limits.max_source_triangles,
            limits.max_candidate_pairs,
            limits.max_regions,
            limits.max_region_vertices,
        );
        let options = ExecutionOptions::new(tolerance);
        // The request borrows the context's meshes; it ends with this block so the
        // guard can be moved into the error record and the result table afterwards.
        let outcome = {
            let state: &Context = &context;
            let inputs = match mesh_inputs(state, &decoded) {
                Ok(inputs) => inputs,
                Err(refusal) => {
                    return record_refusal(context, OPERATION, refusal, tolerance_value)
                }
            };
            let mut request = VisibleFaceRequest::new(&inputs, view, request_limits);
            if let Some(region) = region {
                request = request.with_region(region);
            }
            state.application.visible_faces(&request, &options)
        };
        let outcome = match outcome {
            Ok(outcome) => outcome,
            Err(error) => return record_drawing_error(context, OPERATION, error),
        };
        let Ok(handle) = insert_result(&mut context, StoredResult::VisibleFaces(outcome)) else {
            return AxiolidStatus::LimitExceeded;
        };
        // SAFETY: null was rejected; caller provides one writable result handle.
        unsafe { out_result.write(AxiolidResultHandle(handle)) };
        AxiolidStatus::Ok
    })
}

/// Run `read` on the visible-face outcome behind `result`.
fn with_outcome(
    context: AxiolidContextHandle,
    result: AxiolidResultHandle,
    read: impl FnOnce(&VisibleFaceOutcome) -> AxiolidStatus,
) -> AxiolidStatus {
    let Some(context) = context_entry(context.0) else {
        return AxiolidStatus::InvalidHandle;
    };
    let context = lock_unpoisoned(&context);
    match context.results.get(result.0) {
        Some(StoredResult::VisibleFaces(outcome)) => read(outcome),
        Some(_) => AxiolidStatus::WrongResultKind,
        None => AxiolidStatus::InvalidHandle,
    }
}

/// Query the table sizes of a visible-face result.
/// # Safety
/// Every non-null pointer must be aligned and valid for the documented read or write extent.
#[no_mangle]
pub unsafe extern "C" fn axiolid_v0_4_visible_faces_counts(
    context: AxiolidContextHandle,
    result: AxiolidResultHandle,
    out_counts: *mut AxiolidVisibleFaceCounts,
) -> AxiolidStatus {
    boundary(|| {
        if out_counts.is_null() {
            return AxiolidStatus::NullPointer;
        }
        with_outcome(context, result, |outcome| {
            // SAFETY: null was rejected; caller provides one writable record.
            unsafe { out_counts.write(counts(outcome)) };
            AxiolidStatus::Ok
        })
    })
}

/// Copy the regions, rings, points and triangle indices of a visible-face
/// result into caller-owned storage.
///
/// `points_xy` holds two `f64` per point (x, then y) and `point_cap` counts
/// points. Regions keep the outcome's order; a region's rings are its outer
/// ring then its holes, in order; triangle indices are the region's, ascending.
/// A capacity below its count returns `BufferTooSmall` and writes nothing at
/// all.
/// # Safety
/// Every non-null pointer must be aligned and valid for the documented read or write extent.
#[no_mangle]
pub unsafe extern "C" fn axiolid_v0_4_visible_faces_copy(
    context: AxiolidContextHandle,
    result: AxiolidResultHandle,
    regions: *mut AxiolidVisibleRegion,
    region_cap: usize,
    rings: *mut AxiolidRegionRing,
    ring_cap: usize,
    points_xy: *mut f64,
    point_cap: usize,
    triangles: *mut u32,
    triangle_cap: usize,
) -> AxiolidStatus {
    boundary(|| {
        if (region_cap != 0 && regions.is_null())
            || (ring_cap != 0 && rings.is_null())
            || (point_cap != 0 && points_xy.is_null())
            || (triangle_cap != 0 && triangles.is_null())
        {
            return AxiolidStatus::NullPointer;
        }
        with_outcome(context, result, |outcome| {
            let need = counts(outcome);
            if region_cap < need.regions
                || ring_cap < need.rings
                || point_cap < need.points
                || triangle_cap < need.triangles
            {
                return AxiolidStatus::BufferTooSmall;
            }
            let (mut ring_at, mut point_at, mut triangle_at) = (0_usize, 0_usize, 0_usize);
            for (index, region) in outcome.regions.iter().enumerate() {
                let (source, group) = face_source(region.source);
                let value = AxiolidVisibleRegion {
                    mesh_id: region.mesh.0,
                    first_ring: ring_at,
                    ring_count: 1 + region.holes.len(),
                    first_triangle: triangle_at,
                    triangle_count: region.triangles.len(),
                    group,
                    source,
                };
                // SAFETY: index < regions <= region_cap; the caller promises that many slots.
                unsafe { regions.add(index).write(value) };
                for ring in std::iter::once(&region.outer).chain(&region.holes) {
                    let value = AxiolidRegionRing {
                        first_point: point_at,
                        point_count: ring.len(),
                    };
                    // SAFETY: ring_at < rings <= ring_cap.
                    unsafe { rings.add(ring_at).write(value) };
                    ring_at += 1;
                    for point in ring {
                        // SAFETY: point_at < points <= point_cap; two f64 per point.
                        unsafe {
                            points_xy.add(2 * point_at).write(point.x);
                            points_xy.add(2 * point_at + 1).write(point.y);
                        }
                        point_at += 1;
                    }
                }
                for &triangle in &region.triangles {
                    // SAFETY: triangle_at < triangles <= triangle_cap.
                    unsafe { triangles.add(triangle_at).write(triangle) };
                    triangle_at += 1;
                }
            }
            AxiolidStatus::Ok
        })
    })
}

/// Copy the evidence of a visible-face result.
/// # Safety
/// Every non-null pointer must be aligned and valid for the documented read or write extent.
#[no_mangle]
pub unsafe extern "C" fn axiolid_v0_4_visible_faces_evidence(
    context: AxiolidContextHandle,
    result: AxiolidResultHandle,
    out: *mut AxiolidVisibleFaceEvidence,
) -> AxiolidStatus {
    boundary(|| {
        if out.is_null() {
            return AxiolidStatus::NullPointer;
        }
        with_outcome(context, result, |outcome| {
            let evidence = &outcome.evidence;
            let value = AxiolidVisibleFaceEvidence {
                source_meshes: evidence.source_meshes,
                source_triangles: evidence.source_triangles,
                edge_on_triangles: evidence.edge_on_triangles,
                removed_triangles: evidence.removed_triangles,
                boundary_segments: evidence.boundary_segments,
                intersection_segments: evidence.intersection_segments,
                crossings: evidence.crossings,
                pieces: evidence.pieces,
                tied_sides: evidence.tied_sides,
                candidate_pairs: evidence.candidate_pairs,
                exact_escalations: evidence.exact_escalations,
                regions: evidence.regions,
                cap_regions: evidence.cap_regions,
                rerouted_edges: evidence.rerouted_edges,
                collapsed_regions: evidence.collapsed_regions,
                collapsed_holes: evidence.collapsed_holes,
                snap_displacement: evidence.snap_displacement.get(),
            };
            // SAFETY: null was rejected; caller provides one writable evidence record.
            unsafe { out.write(value) };
            AxiolidStatus::Ok
        })
    })
}
