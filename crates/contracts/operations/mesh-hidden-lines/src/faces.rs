//! Portable visible-face-region provider contract (ADR 0094).
//!
//! The request is a hidden-line request without edge selection: the same
//! meshes, view, crop and error type, plus per-triangle [`FaceGroup`] ids.
//! The outcome lists, per mesh and source, the regions of the image where a
//! face of that source is frontmost. Every decision is exact and no tolerance
//! enters; the contract fixes what the decision means:
//!
//! - a point of the image belongs to source `(mesh, source)` when some kept
//!   surface of that source covers it and no kept surface lies strictly in
//!   front of that surface there;
//! - a kept surface is a triangle of nonzero image area inside
//!   `far <= W <= near` and the crop box, or the cap of a [`MeshRole::Solid`]
//!   cut by `near` ([`FaceSource::Cap`], a surface at depth `near`);
//! - equal depth never occludes, so sources tied over an area each get it
//!   and the evidence counts the tied sides;
//! - a triangle seen edge-on covers no area: it makes no region and never
//!   occludes a region, and is counted in the evidence.
//!
//! A region is a maximal edge-connected open set of one source, given by its
//! closure. Its outer ring runs counter-clockwise and each hole clockwise;
//! every ring starts at its lexicographically smallest point and does not
//! repeat it at the end. A ring vertex is kept unless exactly two boundary
//! pieces meet there in a straight line, so neighbouring regions share
//! bit-identical vertices. Each vertex is an exact point rounded once to the
//! nearest binary64, as for [`HiddenLineSegment`](crate::HiddenLineSegment)
//! endpoints.

use axiolid_contracts::{Backend, CancellationGranularity, ExecutionOptions, ScratchRequirement};
use axiolid_core::{Aabb2, Point2};

use crate::{
    HiddenLineError, HiddenLineMesh, HiddenLineSource, MeshId, MeshRole, OrthographicView,
};

/// Caller-chosen id of a set of triangles that share one surface pattern.
///
/// A material id and a planar face are both encodings the caller picks; the
/// kernel does not infer planes.
#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub struct FaceGroup(pub u32);

/// The face groups of one mesh.
#[non_exhaustive]
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum FaceGroups<'a> {
    /// Every triangle is in `FaceGroup(0)`.
    Single,
    /// One id per triangle; the length equals the triangle count.
    PerTriangle(&'a [FaceGroup]),
}

impl FaceGroups<'_> {
    /// The group of triangle `triangle`.
    ///
    /// # Panics
    ///
    /// Panics for `PerTriangle` when `triangle` is out of range; providers
    /// check the length first.
    pub fn of(&self, triangle: usize) -> FaceGroup {
        match self {
            Self::Single => FaceGroup(0),
            Self::PerTriangle(groups) => groups[triangle],
        }
    }
}

/// What a visible region is a region of.
///
/// The derived order puts every `Group` before `Cap`.
#[non_exhaustive]
#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub enum FaceSource {
    /// The triangles of one group.
    Group(FaceGroup),
    /// The cap of a [`MeshRole::Solid`] cut by the near plane: a surface at
    /// depth `near` covering the image points whose point in that plane lies
    /// inside the solid. A drawing fills it with the cut pattern.
    Cap,
}

/// Hard bounds for one visible-face request.
#[non_exhaustive]
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct VisibleFaceLimits {
    /// Maximum input meshes.
    pub max_meshes: usize,
    /// Maximum source positions, summed over all meshes.
    pub max_source_vertices: usize,
    /// Maximum source triangles, summed over all meshes.
    pub max_source_triangles: usize,
    /// Maximum candidate pairs examined, over all pair work.
    pub max_candidate_pairs: usize,
    /// Maximum output regions.
    pub max_regions: usize,
    /// Maximum ring vertices, summed over all output regions.
    pub max_region_vertices: usize,
}

impl VisibleFaceLimits {
    /// Construct explicit input and output work limits.
    pub const fn new(
        max_meshes: usize,
        max_source_vertices: usize,
        max_source_triangles: usize,
        max_candidate_pairs: usize,
        max_regions: usize,
        max_region_vertices: usize,
    ) -> Self {
        Self {
            max_meshes,
            max_source_vertices,
            max_source_triangles,
            max_candidate_pairs,
            max_regions,
            max_region_vertices,
        }
    }
}

/// One visible-face request.
///
/// All meshes occlude each other and themselves. A `region` crops the image
/// (closed box in image coordinates); geometry outside it neither makes a
/// region nor, being outside the image, occludes anything inside.
#[non_exhaustive]
#[derive(Debug, Clone)]
pub struct VisibleFaceRequest<'a> {
    /// Input meshes; their order orders the output.
    pub meshes: &'a [HiddenLineMesh<'a>],
    /// The orthographic view.
    pub view: OrthographicView,
    /// Optional image crop.
    pub region: Option<Aabb2>,
    /// Work limits.
    pub limits: VisibleFaceLimits,
}

impl<'a> VisibleFaceRequest<'a> {
    /// Find the visible faces of `meshes` in `view` without a crop.
    pub const fn new(
        meshes: &'a [HiddenLineMesh<'a>],
        view: OrthographicView,
        limits: VisibleFaceLimits,
    ) -> Self {
        Self {
            meshes,
            view,
            region: None,
            limits,
        }
    }

    /// Crop the image to a closed box.
    pub const fn with_region(mut self, region: Aabb2) -> Self {
        self.region = Some(region);
        self
    }
}

/// A maximal region of the image where one source is frontmost.
///
/// Rings are counter-clockwise (outer) and clockwise (holes), start at their
/// lexicographically smallest point and do not repeat it.
#[non_exhaustive]
#[derive(Debug, Clone, PartialEq)]
pub struct VisibleRegion {
    /// The mesh the region belongs to.
    pub mesh: MeshId,
    /// The group or cap the region belongs to.
    pub source: FaceSource,
    /// Outer boundary.
    pub outer: Vec<Point2>,
    /// Holes, sorted lexicographically by their points.
    pub holes: Vec<Vec<Point2>>,
    /// Indices of the mesh's triangles that are frontmost somewhere in the
    /// region, ascending. Empty for [`FaceSource::Cap`].
    pub triangles: Vec<u32>,
}

impl VisibleRegion {
    /// Construct a region of `source` of `mesh`.
    pub fn new(
        mesh: MeshId,
        source: FaceSource,
        outer: Vec<Point2>,
        holes: Vec<Vec<Point2>>,
        triangles: Vec<u32>,
    ) -> Self {
        Self {
            mesh,
            source,
            outer,
            holes,
            triangles,
        }
    }
}

/// Auditable counts for one visible-face operation.
#[non_exhaustive]
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct VisibleFaceEvidence {
    /// Approximation source.
    pub source: HiddenLineSource,
    /// Input meshes.
    pub source_meshes: usize,
    /// Input triangles, summed over all meshes.
    pub source_triangles: usize,
    /// Triangles seen edge-on: they cover no area and make no region.
    pub edge_on_triangles: usize,
    /// Triangles removed entirely by depth clipping or the crop.
    pub removed_triangles: usize,
    /// Boundary pieces of the arrangement that come from triangle edges and
    /// cap outlines.
    pub boundary_segments: usize,
    /// Input segments of the arrangement where a triangle of one surface cuts
    /// another, before the arrangement cuts them into pieces.
    pub intersection_segments: usize,
    /// Crossings of segments in the image.
    pub crossings: usize,
    /// Pieces of the arrangement.
    pub pieces: usize,
    /// Sides of arrangement pieces that two or more sources share by equal
    /// depth, interior pieces included. The count depends on how the
    /// arrangement cuts the plane, so only whether it is zero is a property
    /// of the scene.
    pub tied_sides: usize,
    /// Candidate pairs examined, over all pair work.
    pub candidate_pairs: usize,
    /// Decisions the fast filter could not settle and the exact tier did.
    /// This count may differ between two otherwise identical requests; it is
    /// not part of the result.
    pub exact_escalations: usize,
    /// Output regions.
    pub regions: usize,
    /// Output regions whose source is [`FaceSource::Cap`].
    pub cap_regions: usize,
}

impl VisibleFaceEvidence {
    /// Evidence for a result derived from the input meshes; the remaining
    /// counts start at zero for the provider to fill in.
    pub const fn input_mesh(source_meshes: usize, source_triangles: usize) -> Self {
        Self {
            source: HiddenLineSource::InputMesh,
            source_meshes,
            source_triangles,
            edge_on_triangles: 0,
            removed_triangles: 0,
            boundary_segments: 0,
            intersection_segments: 0,
            crossings: 0,
            pieces: 0,
            tied_sides: 0,
            candidate_pairs: 0,
            exact_escalations: 0,
            regions: 0,
            cap_regions: 0,
        }
    }

    /// Whether this result came from the discrete input meshes.
    pub const fn is_derived_from_input_mesh(self) -> bool {
        matches!(self.source, HiddenLineSource::InputMesh)
    }
}

/// Visible regions plus evidence.
///
/// Regions are ordered by the mesh's position in the request, then by
/// `source`, then by the outer ring's points, lexicographically.
#[non_exhaustive]
#[derive(Debug, Clone, PartialEq)]
pub struct VisibleFaceOutcome {
    /// The view the regions are expressed in, unchanged.
    pub view: OrthographicView,
    /// The regions.
    pub regions: Vec<VisibleRegion>,
    /// Counts.
    pub evidence: VisibleFaceEvidence,
}

impl VisibleFaceOutcome {
    /// Construct a provider result. Registries validate it before returning it.
    pub fn new(
        view: OrthographicView,
        regions: Vec<VisibleRegion>,
        evidence: VisibleFaceEvidence,
    ) -> Self {
        Self {
            view,
            regions,
            evidence,
        }
    }
}

/// Provider for deterministic visible-face regions of triangle meshes.
///
/// Provider obligations beyond the types:
///
/// - a [`FaceSource::Cap`] region exists only for a [`MeshRole::Solid`] that
///   `view.near` cuts, and has no `triangles`;
/// - regions of one source never overlap; regions of different sources
///   overlap only where they tie in depth;
/// - flipping every triangle gives identical regions and identical evidence
///   except `exact_escalations`;
/// - ring vertices are the exact points rounded once to nearest;
/// - the tolerance of `options` applies to mesh audit and frame validation
///   only, never to a visibility decision;
/// - the same input refusals as [`MeshHiddenLines`](crate::MeshHiddenLines),
///   and a `PerTriangle` group slice whose length differs from the triangle
///   count is `InvalidInput` naming the mesh; no substitute result is
///   returned;
/// - output storage is outside the scratch budget.
pub trait MeshVisibleFaces: Backend {
    /// Scratch needed beyond inputs and result.
    fn scratch_requirement(&self) -> ScratchRequirement {
        ScratchRequirement::Unbounded
    }

    /// How finely the provider polls cancellation.
    fn cancellation_granularity(&self) -> CancellationGranularity {
        CancellationGranularity::None
    }

    /// Find the visible face regions of `request.meshes` in `request.view`.
    /// Implementations must enforce `request.limits` before growing output.
    fn visible_faces(
        &self,
        request: &VisibleFaceRequest<'_>,
        options: &ExecutionOptions,
    ) -> Result<VisibleFaceOutcome, HiddenLineError>;
}

#[cfg(test)]
mod tests {
    use super::*;
    use axiolid_mesh::TriMesh;

    #[test]
    fn single_group_is_zero_for_every_triangle() {
        assert_eq!(FaceGroups::Single.of(0), FaceGroup(0));
        assert_eq!(FaceGroups::Single.of(41), FaceGroup(0));
    }

    #[test]
    fn per_triangle_groups_are_indexed() {
        let ids = [FaceGroup(5), FaceGroup(2)];
        let groups = FaceGroups::PerTriangle(&ids);
        assert_eq!(groups.of(0), FaceGroup(5));
        assert_eq!(groups.of(1), FaceGroup(2));
    }

    #[test]
    fn mesh_defaults_to_single_and_with_groups_replaces_it() {
        let mesh = TriMesh::default();
        let m = HiddenLineMesh::new(MeshId(1), MeshRole::Sheet, &mesh);
        assert_eq!(m.groups, FaceGroups::Single);
        let ids = [FaceGroup(3)];
        let m = m.with_groups(&ids);
        assert_eq!(m.groups, FaceGroups::PerTriangle(&ids));
    }

    #[test]
    fn groups_order_before_the_cap() {
        let mut sources = [
            FaceSource::Cap,
            FaceSource::Group(FaceGroup(9)),
            FaceSource::Group(FaceGroup(1)),
        ];
        sources.sort();
        assert_eq!(
            sources,
            [
                FaceSource::Group(FaceGroup(1)),
                FaceSource::Group(FaceGroup(9)),
                FaceSource::Cap
            ]
        );
    }

    #[test]
    fn input_mesh_evidence_starts_at_zero() {
        let e = VisibleFaceEvidence::input_mesh(2, 12);
        assert!(e.is_derived_from_input_mesh());
        assert_eq!((e.source_meshes, e.source_triangles), (2, 12));
        let counts = [
            e.edge_on_triangles,
            e.removed_triangles,
            e.boundary_segments,
            e.intersection_segments,
            e.crossings,
            e.pieces,
            e.tied_sides,
            e.candidate_pairs,
            e.exact_escalations,
            e.regions,
            e.cap_regions,
        ];
        assert!(counts.iter().all(|&c| c == 0));
    }
}
