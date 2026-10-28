//! The scalar exact hidden-line provider.

use std::collections::BTreeSet;

use axiolid_contracts::{
    Backend, BackendDescriptor, BackendId, CancellationGranularity, ExecutionOptions,
    ExecutionTarget, GeomError, Operation, ScratchRequirement,
};
use axiolid_core::{Frame3, SpaceFrame, Tolerance};
use axiolid_mesh::try_audit_mesh;
use axiolid_mesh_hidden_lines_contract::{
    HiddenLineError, HiddenLineEvidence, HiddenLineMesh, HiddenLineOutcome, HiddenLineRequest,
    HiddenLineSegment, MeshHiddenLines, MeshRole, SourceEdge, Visibility,
};

use crate::edges::{collect, DrawEdge, MeshEdges, Tri};
use crate::occlusion::{resolve, Resolved, Scene, Work};
use crate::output::endpoint;
use crate::view::Vert;

/// Portable exact orthographic hidden-line provider.
///
/// Visibility is decided exactly in the view's affine frame: nothing is
/// rounded before a decision, only the endpoints of the output pieces are
/// (each the exact point rounded once). It tests every edge against every
/// triangle whose bounding box can matter, so the cost grows with edges times
/// triangles; `max_candidate_pairs` bounds it.
///
/// Refused by name: a near clipping plane with any [`MeshRole::Solid`] in the
/// request (the cap of a cut solid occludes, and is not built yet).
#[derive(Debug, Default, Clone, Copy)]
pub struct ScalarHiddenLines;

impl ScalarHiddenLines {
    /// Stable identity of this provider.
    pub const ID: BackendId = BackendId::new("scalar-hidden-lines");

    /// Construct the provider.
    #[must_use]
    pub const fn new() -> Self {
        Self
    }
}

impl Backend for ScalarHiddenLines {
    fn descriptor(&self) -> BackendDescriptor {
        BackendDescriptor::new(Self::ID, ExecutionTarget::PortableCpu)
    }
}

impl MeshHiddenLines for ScalarHiddenLines {
    fn scratch_requirement(&self) -> ScratchRequirement {
        // Exact view-space vertices, triangles, edge records and the spans
        // of one edge.
        ScratchRequirement::PerElement {
            bytes_per_element: 1024,
        }
    }

    fn cancellation_granularity(&self) -> CancellationGranularity {
        CancellationGranularity::Incremental
    }

    fn hidden_lines(
        &self,
        request: &HiddenLineRequest<'_>,
        options: &ExecutionOptions,
    ) -> Result<HiddenLineOutcome, HiddenLineError> {
        options.check_cancelled()?;
        let (total_vertices, total_triangles) = validate_shape(request, options)?;
        for entry in request.meshes {
            options.check_cancelled()?;
            audit(entry, options)?;
        }
        if request.view.near.is_some() {
            if let Some(solid) = request
                .meshes
                .iter()
                .find(|entry| entry.role == MeshRole::Solid)
            {
                return Err(HiddenLineError::mesh(
                    solid.id,
                    GeomError::UnsupportedInput {
                        backend: Self::ID,
                        operation: Operation::MeshHiddenLines,
                        input: "near clip through a capped solid",
                    },
                ));
            }
        }

        let memory = || HiddenLineError::from(GeomError::BudgetExceeded { resource: "memory" });
        let mut verts: Vec<Vert> = Vec::new();
        verts
            .try_reserve_exact(total_vertices)
            .map_err(|_| memory())?;
        let mut tris: Vec<Tri> = Vec::new();
        let mut edges: Vec<DrawEdge> = Vec::new();
        let mut escalations = 0;
        let frame = &request.view.frame;
        for (index, entry) in request.meshes.iter().enumerate() {
            let base = u32::try_from(verts.len()).map_err(|_| {
                HiddenLineError::from(GeomError::BudgetExceeded {
                    resource: "hidden-line source vertices",
                })
            })?;
            verts.extend(entry.mesh.positions.iter().map(|p| Vert::new(*p, frame)));
            let input = MeshEdges {
                index,
                id: entry.id,
                role: entry.role,
                mesh: entry.mesh,
                base,
            };
            collect(
                &input,
                &verts,
                &request.edges,
                options,
                &mut escalations,
                &mut tris,
                &mut edges,
            )?;
        }
        debug_assert_eq!(tris.len(), total_triangles);

        let scene = Scene {
            verts: &verts,
            tris: &tris,
            near: request.view.near,
            far: request.view.far,
            region: request.region,
        };
        let mut evidence = HiddenLineEvidence::input_mesh(request.meshes.len(), total_triangles);
        evidence.selected_edges = edges.len();
        let mut work = Work::default();
        let mut segments: Vec<HiddenLineSegment> = Vec::new();
        for edge in &edges {
            options.check_cancelled()?;
            let id = request.meshes[edge.mesh].id;
            let resolved = resolve(&scene, edge, request.limits.max_candidate_pairs, &mut work)
                .map_err(|error| match error {
                    GeomError::Degenerate(_) => HiddenLineError::mesh(id, error),
                    other => HiddenLineError::from(other),
                })?;
            let pieces = match resolved {
                Resolved::Point => {
                    evidence.point_edges += 1;
                    continue;
                }
                Resolved::Clipped => {
                    evidence.clipped_edges += 1;
                    continue;
                }
                Resolved::Pieces(pieces) => pieces,
            };
            let [start, end] = edge.ends;
            for piece in pieces {
                if piece.hidden && !request.edges.hidden {
                    continue;
                }
                let from = endpoint(&verts, start, end, piece.lo)?;
                let to = endpoint(&verts, start, end, piece.hi)?;
                if from == to {
                    evidence.collapsed_segments += 1;
                    continue;
                }
                if segments.len() >= request.limits.max_output_segments {
                    return Err(GeomError::BudgetExceeded {
                        resource: "hidden-line output segments",
                    }
                    .into());
                }
                segments.try_reserve(1).map_err(|_| memory())?;
                let visibility = if piece.hidden {
                    evidence.hidden_segments += 1;
                    Visibility::Hidden
                } else {
                    evidence.visible_segments += 1;
                    Visibility::Visible
                };
                segments.push(HiddenLineSegment::new(
                    SourceEdge {
                        mesh: id,
                        vertices: edge.local,
                    },
                    edge.kind,
                    visibility,
                    from,
                    to,
                ));
            }
        }
        evidence.candidate_pairs = work.candidate_pairs;
        evidence.exact_escalations = work.escalations + escalations;
        Ok(HiddenLineOutcome::new(request.view, segments, evidence))
    }
}

/// Request-level checks. Returns the total positions and triangles.
fn validate_shape(
    request: &HiddenLineRequest<'_>,
    options: &ExecutionOptions,
) -> Result<(usize, usize), HiddenLineError> {
    let invalid = |detail: &str| HiddenLineError::from(GeomError::InvalidInput(detail.into()));
    let budget = |resource| HiddenLineError::from(GeomError::BudgetExceeded { resource });
    let limits = request.limits;
    if request.meshes.is_empty() {
        return Err(invalid("hidden-line request has no meshes"));
    }
    if request.meshes.len() > limits.max_meshes {
        return Err(budget("hidden-line meshes"));
    }
    let mut vertices = 0usize;
    let mut triangles = 0usize;
    for entry in request.meshes {
        vertices = vertices
            .checked_add(entry.mesh.positions.len())
            .ok_or_else(|| budget("hidden-line source vertices"))?;
        triangles = triangles
            .checked_add(entry.mesh.triangle_count())
            .ok_or_else(|| budget("hidden-line source triangles"))?;
    }
    if vertices > limits.max_source_vertices {
        return Err(budget("hidden-line source vertices"));
    }
    if triangles > limits.max_source_triangles {
        return Err(budget("hidden-line source triangles"));
    }

    let mut ids = BTreeSet::new();
    for entry in request.meshes {
        if !ids.insert(entry.id) {
            return Err(HiddenLineError::mesh(
                entry.id,
                GeomError::InvalidInput("duplicate mesh id".into()),
            ));
        }
        if entry.mesh.indices.is_empty() || !entry.mesh.indices.len().is_multiple_of(3) {
            return Err(HiddenLineError::mesh(
                entry.id,
                GeomError::InvalidInput(
                    "hidden-line mesh has no triangles or a partial triangle".into(),
                ),
            ));
        }
    }

    let view = &request.view;
    validate_frame(&view.frame, options.tolerance())?;
    for (name, bound) in [("near", view.near), ("far", view.far)] {
        if bound.is_some_and(|value| !value.is_finite()) {
            return Err(HiddenLineError::from(GeomError::InvalidInput(format!(
                "hidden-line {name} depth is not finite"
            ))));
        }
    }
    if let (Some(near), Some(far)) = (view.near, view.far) {
        if far >= near {
            return Err(invalid(
                "hidden-line far depth must be smaller than the near depth",
            ));
        }
    }
    if request
        .edges
        .crease_cosine
        .is_some_and(|cosine| !(-1.0..=1.0).contains(&cosine))
    {
        return Err(invalid("hidden-line crease cosine must lie in [-1, 1]"));
    }
    if let Some(region) = request.region {
        let finite = region.min.is_finite() && region.max.is_finite();
        if !finite || region.min.x > region.max.x || region.min.y > region.max.y {
            return Err(invalid(
                "hidden-line region must be finite with min not above max",
            ));
        }
    }
    Ok((vertices, triangles))
}

/// A finite right-handed orthonormal frame, as the registry requires.
fn validate_frame(frame: &Frame3, tolerance: Tolerance) -> Result<(), HiddenLineError> {
    // A coarse caller tolerance must not admit zero or scaled axes.
    let limit = tolerance.angular().clamp(128.0 * f64::EPSILON, 1.0e-6);
    let tolerance = Tolerance::new(limit, limit).map_err(|_| {
        HiddenLineError::from(GeomError::InvalidInput(
            "hidden-line view frame tolerance is invalid".into(),
        ))
    })?;
    SpaceFrame::new(frame.origin, frame.x, frame.y, frame.z, tolerance)
        .map(|_| ())
        .map_err(|error| {
            HiddenLineError::from(GeomError::InvalidInput(format!(
                "hidden-line view frame is invalid: {error}"
            )))
        })
}

/// Per-mesh audit: indices, positions, degenerate triangles, topology.
fn audit(entry: &HiddenLineMesh<'_>, options: &ExecutionOptions) -> Result<(), HiddenLineError> {
    let refuse = |error: GeomError| HiddenLineError::mesh(entry.id, error);
    let health = try_audit_mesh(entry.mesh, options.tolerance())
        .map_err(|_| refuse(GeomError::BudgetExceeded { resource: "memory" }))?;
    if let Some((triangle, index)) = health.first_invalid_index {
        return Err(refuse(GeomError::InvalidInput(format!(
            "triangle {triangle} references position {index}, but the mesh has {}",
            health.positions
        ))));
    }
    if let Some(position) = health.first_non_finite_position {
        return Err(refuse(GeomError::InvalidInput(format!(
            "position {position} is not finite"
        ))));
    }
    if health.degenerate_triangles > 0 {
        return Err(refuse(GeomError::Degenerate(format!(
            "{} degenerate triangle(s)",
            health.degenerate_triangles
        ))));
    }
    if health.non_manifold_edges > 0 || health.inconsistent_winding_edges > 0 {
        return Err(refuse(GeomError::NotManifold(format!(
            "non_manifold={}, inconsistent_winding={}",
            health.non_manifold_edges, health.inconsistent_winding_edges
        ))));
    }
    if entry.role == MeshRole::Solid && health.boundary_edges > 0 {
        return Err(refuse(GeomError::NotManifold(format!(
            "a solid must be closed; boundary={}",
            health.boundary_edges
        ))));
    }
    Ok(())
}
