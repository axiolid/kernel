//! The scalar exact visible-face provider.

use std::cmp::Ordering;

use axiolid_contracts::{
    Backend, BackendDescriptor, BackendId, CancellationGranularity, ExecutionOptions,
    ExecutionTarget, GeomError, ScratchRequirement,
};
use axiolid_exact_arrangement::{arrange, PairBudget, Segment};
use axiolid_mesh_hidden_lines_contract::{
    FaceGroups, HiddenLineError, MeshVisibleFaces, VisibleFaceEvidence, VisibleFaceOutcome,
    VisibleFaceRequest,
};

use super::intersect::intersection_segments;
use super::label::label_piece;
use super::points::Points;
use super::regions::{build, OutputLimits};
use super::segments::{clip_segments, crop_segments, edge_segments, PAIRS};
use super::{Faces, Work};
use crate::broad::{Broad, BroadPhase};
use crate::prepare::{prepare, Prepared, ResourceNames, SharedLimits};

/// Bytes of scratch per source triangle that the provider declares.
const SCRATCH_PER_ELEMENT: usize = 3072;

/// Portable exact orthographic visible-face provider.
///
/// Visibility is decided exactly (module docs of `faces`): the image is
/// cut by an exact arrangement of the visible mesh edges, the clip segments,
/// the cross-sections of interpenetrating triangles and the crop box; each
/// side of each piece is labelled with the sources frontmost just beside it;
/// and the pieces between different labels bound the regions. Nothing is
/// rounded before a decision; each region vertex is the exact point rounded
/// once, so it equals the hidden-line endpoint of the same point.
///
/// Where several sources are at equal depth over an area they all get it
/// and the evidence counts the tied sides. A triangle seen edge-on covers no
/// area and makes no region. A near plane that cuts a
/// [`MeshRole::Solid`](axiolid_mesh_hidden_lines_contract::MeshRole::Solid)
/// leaves a [`FaceSource::Cap`](axiolid_mesh_hidden_lines_contract::FaceSource::Cap)
/// region. A ring that collapses under rounding is `Degenerate` and names its
/// mesh.
#[derive(Debug, Default, Clone, Copy)]
pub struct ScalarVisibleFaces;

impl ScalarVisibleFaces {
    /// Stable identity of this provider.
    pub const ID: BackendId = BackendId::new("scalar-visible-faces");

    /// Construct the provider.
    #[must_use]
    pub const fn new() -> Self {
        Self
    }
}

impl Backend for ScalarVisibleFaces {
    fn descriptor(&self) -> BackendDescriptor {
        BackendDescriptor::new(Self::ID, ExecutionTarget::PortableCpu)
    }
}

impl MeshVisibleFaces for ScalarVisibleFaces {
    fn scratch_requirement(&self) -> ScratchRequirement {
        ScratchRequirement::PerElement {
            bytes_per_element: SCRATCH_PER_ELEMENT,
        }
    }

    fn cancellation_granularity(&self) -> CancellationGranularity {
        CancellationGranularity::Incremental
    }

    fn visible_faces(
        &self,
        request: &VisibleFaceRequest<'_>,
        options: &ExecutionOptions,
    ) -> Result<VisibleFaceOutcome, HiddenLineError> {
        run(BroadPhase::Tree, true, request, options)
    }
}

/// The provider with a linear-scan broad phase. Its outcome equals that of
/// [`ScalarVisibleFaces`]. Not part of the supported API.
#[doc(hidden)]
#[derive(Debug, Default, Clone, Copy)]
pub struct ScanVisibleFaces;

impl Backend for ScanVisibleFaces {
    fn descriptor(&self) -> BackendDescriptor {
        BackendDescriptor::new(
            BackendId::new("scalar-visible-faces-scan"),
            ExecutionTarget::PortableCpu,
        )
    }
}

impl MeshVisibleFaces for ScanVisibleFaces {
    fn scratch_requirement(&self) -> ScratchRequirement {
        ScalarVisibleFaces.scratch_requirement()
    }

    fn cancellation_granularity(&self) -> CancellationGranularity {
        CancellationGranularity::Incremental
    }

    fn visible_faces(
        &self,
        request: &VisibleFaceRequest<'_>,
        options: &ExecutionOptions,
    ) -> Result<VisibleFaceOutcome, HiddenLineError> {
        run(BroadPhase::Scan, true, request, options)
    }
}

/// The provider that keeps hidden pieces of edges in the arrangement. Its
/// regions equal those of [`ScalarVisibleFaces`]; it is the proof that the
/// pruning is exact and the baseline the benchmark times it against. Not part
/// of the supported API.
#[doc(hidden)]
#[derive(Debug, Default, Clone, Copy)]
pub struct UnprunedVisibleFaces;

impl Backend for UnprunedVisibleFaces {
    fn descriptor(&self) -> BackendDescriptor {
        BackendDescriptor::new(
            BackendId::new("scalar-visible-faces-unpruned"),
            ExecutionTarget::PortableCpu,
        )
    }
}

impl MeshVisibleFaces for UnprunedVisibleFaces {
    fn scratch_requirement(&self) -> ScratchRequirement {
        ScalarVisibleFaces.scratch_requirement()
    }

    fn cancellation_granularity(&self) -> CancellationGranularity {
        CancellationGranularity::Incremental
    }

    fn visible_faces(
        &self,
        request: &VisibleFaceRequest<'_>,
        options: &ExecutionOptions,
    ) -> Result<VisibleFaceOutcome, HiddenLineError> {
        run(BroadPhase::Tree, false, request, options)
    }
}

fn run(
    phase: BroadPhase,
    prune: bool,
    request: &VisibleFaceRequest<'_>,
    options: &ExecutionOptions,
) -> Result<VisibleFaceOutcome, HiddenLineError> {
    options.check_cancelled()?;
    for entry in request.meshes {
        if let FaceGroups::PerTriangle(groups) = entry.groups {
            if groups.len() != entry.mesh.triangle_count() {
                return Err(HiddenLineError::mesh(
                    entry.id,
                    GeomError::InvalidInput(format!(
                        "{} face groups for {} triangles",
                        groups.len(),
                        entry.mesh.triangle_count()
                    )),
                ));
            }
        }
    }
    let limits = request.limits;
    let prepared = prepare(
        request.meshes,
        &request.view,
        request.region,
        SharedLimits {
            max_meshes: limits.max_meshes,
            max_source_vertices: limits.max_source_vertices,
            max_source_triangles: limits.max_source_triangles,
        },
        &ResourceNames::VISIBLE_FACES,
        options,
    )?;
    let Prepared {
        verts,
        tris,
        total_triangles,
        escalations,
        ..
    } = &prepared;
    let broad = Broad::build(tris, phase);
    let mut work = Work {
        budget: PairBudget::new(limits.max_candidate_pairs, PAIRS),
        escalations: *escalations,
    };
    let faces = Faces::new(
        &prepared,
        &broad,
        request.meshes,
        (request.view.near, request.view.far, request.region),
        &mut work,
        options,
    )?;

    let mut points = Points::new(verts.len());
    let mut segments: Vec<Segment> = Vec::new();
    edge_segments(
        &faces,
        &mut points,
        prune,
        &mut work,
        options,
        &mut segments,
    )?;
    clip_segments(&faces, &mut points, &mut work, options, &mut segments)?;
    crop_segments(&faces, &mut points, &mut segments)?;
    let outline_segments = segments.len();
    let per_element = SCRATCH_PER_ELEMENT.saturating_mul(*total_triangles);
    intersection_segments(
        &faces,
        &mut points,
        &mut work,
        per_element,
        options,
        &mut segments,
    )?;
    let arrangement = arrange(
        &mut points.store,
        &segments,
        &mut work.budget,
        per_element,
        options,
    )?;

    let mut labels = Vec::with_capacity(arrangement.pieces.len());
    let mut evidence = VisibleFaceEvidence::input_mesh(request.meshes.len(), *total_triangles);
    for (index, piece) in arrangement.pieces.iter().enumerate() {
        if index % 256 == 0 {
            options.check_cancelled()?;
        }
        let ends = piece.ends.map(|v| arrangement.vertices[v.0].point);
        let label = label_piece(&faces, &points.store, ends, &mut work)?;
        evidence.tied_sides += label.iter().filter(|l| l.keys.len() > 1).count();
        if label[0].keys != label[1].keys && piece.segments[0] < outline_segments {
            evidence.boundary_segments += 1;
        }
        labels.push(label);
    }

    let regions = build(
        &faces,
        &points.store,
        &arrangement,
        &labels,
        OutputLimits {
            max_regions: limits.max_regions,
            max_region_vertices: limits.max_region_vertices,
        },
        &mut work,
        options,
    )?;

    evidence.edge_on_triangles = tris.iter().filter(|t| t.sigma == Ordering::Equal).count();
    evidence.removed_triangles = faces.removed.iter().filter(|&&r| r).count();
    evidence.intersection_segments = segments.len() - outline_segments;
    evidence.crossings = arrangement.crossings;
    evidence.pieces = arrangement.pieces.len();
    evidence.candidate_pairs = work.budget.used;
    evidence.exact_escalations = work.escalations + points.store.escalations();
    evidence.regions = regions.len();
    evidence.cap_regions = regions
        .iter()
        .filter(|r| r.source == axiolid_mesh_hidden_lines_contract::FaceSource::Cap)
        .count();
    Ok(VisibleFaceOutcome::new(request.view, regions, evidence))
}
