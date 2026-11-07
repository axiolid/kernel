//! The scalar exact hidden-line provider.

use axiolid_contracts::{
    Backend, BackendDescriptor, BackendId, CancellationGranularity, ExecutionOptions,
    ExecutionTarget, GeomError, ScratchRequirement,
};
use axiolid_mesh_hidden_lines_contract::{
    HiddenLineError, HiddenLineEvidence, HiddenLineOutcome, HiddenLineRequest, HiddenLineSegment,
    MeshHiddenLines, SourceEdge, ViewProjection, Visibility,
};

use crate::broad::{Broad, BroadPhase};
use crate::edges::{collect, DrawEdge, MeshEdges};
use crate::groups::resolve_groups;
use crate::occlusion::{resolve, Resolved, Scene, Work};
use crate::output::endpoint;
use crate::prepare::{prepare, Prepared, ResourceNames, SharedLimits};

/// Portable exact hidden-line provider for orthographic and perspective views.
///
/// Visibility is decided exactly in the view's affine frame: nothing is
/// rounded before a decision, only the endpoints of the output pieces are
/// (each the exact point rounded once). For every edge it examines the
/// triangles whose view-space bounding box can matter, found in an
/// `axiolid-spatial` BVH; in the worst case (everything overlaps in the
/// image) the cost still grows with edges times triangles, and
/// `max_candidate_pairs` bounds it.
///
/// A near clipping plane that cuts a [`MeshRole::Solid`](axiolid_mesh_hidden_lines_contract::MeshRole::Solid) leaves a cap that
/// occludes (decided by ray parity with a symbolic perturbation of the image
/// point); a solid wholly beyond the plane is removed: it neither occludes
/// nor is drawn. A failure to produce an endpoint (a coordinate beyond the
/// binary64 range) is `Degenerate` and names the mesh of the edge.
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

    fn supports_projection(&self, _projection: ViewProjection) -> bool {
        true
    }

    fn hidden_lines(
        &self,
        request: &HiddenLineRequest<'_>,
        options: &ExecutionOptions,
    ) -> Result<HiddenLineOutcome, HiddenLineError> {
        draw(BroadPhase::Tree, request, options)
    }
}

fn draw(
    phase: BroadPhase,
    request: &HiddenLineRequest<'_>,
    options: &ExecutionOptions,
) -> Result<HiddenLineOutcome, HiddenLineError> {
    options.check_cancelled()?;
    if request
        .edges
        .crease_cosine
        .is_some_and(|cosine| !(-1.0..=1.0).contains(&cosine))
    {
        return Err(GeomError::InvalidInput(
            "hidden-line crease cosine must lie in [-1, 1]".into(),
        )
        .into());
    }
    let limits = request.limits;
    let Prepared {
        verts,
        tris,
        caps,
        ranges,
        mut escalations,
        total_triangles,
        ..
    } = prepare(
        request.meshes,
        &request.view,
        request.region,
        SharedLimits {
            max_meshes: limits.max_meshes,
            max_source_vertices: limits.max_source_vertices,
            max_source_triangles: limits.max_source_triangles,
        },
        &ResourceNames::HIDDEN_LINES,
        options,
    )?;
    let memory = || HiddenLineError::from(GeomError::BudgetExceeded { resource: "memory" });
    let mut edges: Vec<DrawEdge> = Vec::new();
    for (index, (entry, range)) in request.meshes.iter().zip(&ranges).enumerate() {
        // Groups are read only when their boundaries are drawn.
        let groups = if request.edges.group_boundary {
            Some(resolve_groups(entry, options, &mut escalations)?)
        } else {
            None
        };
        let input = MeshEdges {
            index,
            id: entry.id,
            role: entry.role,
            mesh: entry.mesh,
            base: range.base,
            first_tri: range.triangles.start,
            groups: groups.as_deref(),
        };
        collect(
            &input,
            &verts,
            &tris,
            &request.edges,
            options,
            &mut escalations,
            &mut edges,
        )?;
    }

    let broad = Broad::build(&tris, phase);
    let scene = Scene {
        verts: &verts,
        tris: &tris,
        near: request.view.near,
        far: request.view.far,
        region: request.region,
        caps: &caps,
        broad: &broad,
        edge_on_occludes: true,
        perspective: request.view.projection == ViewProjection::Perspective,
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
        let source = SourceEdge {
            mesh: id,
            vertices: edge.local,
        };
        let first = segments.len();
        // A collapsed piece since the last pushed segment: its neighbours
        // of one visibility are one piece after rounding.
        let mut collapsed = false;
        for piece in pieces {
            if piece.hidden && !request.edges.hidden {
                collapsed = false;
                continue;
            }
            let round = |param| {
                endpoint(&verts, start, end, param)
                    .map_err(|error| HiddenLineError::mesh(id, error))
            };
            let (from, to) = (round(piece.lo)?, round(piece.hi)?);
            if from == to {
                evidence.collapsed_segments += 1;
                collapsed = true;
                continue;
            }
            let visibility = if piece.hidden {
                Visibility::Hidden
            } else {
                Visibility::Visible
            };
            if collapsed && segments.len() > first {
                let count = segments.len();
                let last = &mut segments[count - 1];
                if last.visibility == visibility {
                    last.end = to;
                    collapsed = false;
                    continue;
                }
            }
            collapsed = false;
            if segments.len() >= request.limits.max_output_segments {
                return Err(GeomError::BudgetExceeded {
                    resource: "hidden-line output segments",
                }
                .into());
            }
            segments.try_reserve(1).map_err(|_| memory())?;
            segments.push(HiddenLineSegment::new(
                source, edge.kind, visibility, from, to,
            ));
        }
    }
    for segment in &segments {
        match segment.visibility {
            Visibility::Hidden => evidence.hidden_segments += 1,
            _ => evidence.visible_segments += 1,
        }
    }
    evidence.candidate_pairs = work.candidate_pairs;
    evidence.exact_escalations = work.escalations + escalations;
    Ok(HiddenLineOutcome::new(request.view, segments, evidence))
}

/// The scalar provider with a linear-scan broad phase: the baseline the
/// benchmark times the BVH against, and the reference its candidates must
/// equal. Its output is identical to [`ScalarHiddenLines`]; it is slower on
/// large scenes. Not part of the supported API.
#[doc(hidden)]
#[derive(Debug, Default, Clone, Copy)]
pub struct ScanHiddenLines;

impl Backend for ScanHiddenLines {
    fn descriptor(&self) -> BackendDescriptor {
        BackendDescriptor::new(
            BackendId::new("scalar-hidden-lines-scan"),
            ExecutionTarget::PortableCpu,
        )
    }
}

impl MeshHiddenLines for ScanHiddenLines {
    fn scratch_requirement(&self) -> ScratchRequirement {
        ScalarHiddenLines.scratch_requirement()
    }

    fn cancellation_granularity(&self) -> CancellationGranularity {
        CancellationGranularity::Incremental
    }

    fn supports_projection(&self, _projection: ViewProjection) -> bool {
        true
    }

    fn hidden_lines(
        &self,
        request: &HiddenLineRequest<'_>,
        options: &ExecutionOptions,
    ) -> Result<HiddenLineOutcome, HiddenLineError> {
        draw(BroadPhase::Scan, request, options)
    }
}
