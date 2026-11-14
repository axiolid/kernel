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
//!   `far <= W <= near` and the crop box, or the cap of a [`MeshRole::Solid`](crate::MeshRole::Solid)
//!   cut by `near` ([`FaceSource::Cap`], a surface at depth `near`);
//! - equal depth never occludes, so sources tied over an area each get it
//!   and the evidence counts the tied sides;
//! - a triangle seen edge-on covers no area: it makes no region and never
//!   occludes a region, and is counted in the evidence;
//! - a triangle flagged [`Occlusion::SeeThrough`](crate::Occlusion::SeeThrough)
//!   is a kept surface like any other and makes regions where it is frontmost,
//!   but those regions list what is seen through it (ADR 0111).
//!
//! What is seen through a see-through surface is its stack. The kept surfaces
//! covering an image point are grouped into depth layers by exactly equal
//! depth, layer 0 frontmost. The seen stack is every layer up to and
//! including the first that holds an occluding surface (an opaque triangle or
//! a cap), or every layer if none does. A region of a see-through source
//! lists the seen stack without its own layer-0 entry in
//! [`VisibleRegion::behind`]; a region of an opaque source or a cap lists
//! nothing. Within a layer there is one entry per distinct `(mesh, source)`,
//! ordered by the mesh's position in the request and then by source; a source
//! crossed at two depths is listed once per depth. Every face group must be
//! all opaque or all see-through, else the request is refused naming the
//! mesh.
//!
//! A region is a maximal edge-connected open set of one source and one
//! `behind` list, given by its closure. Its outer ring runs counter-clockwise and each hole clockwise;
//! every ring starts at its lexicographically smallest point and does not
//! repeat it at the end.
//!
//! The regions are the exact regions snap-rounded to the binary64 grid (ADR
//! 0100). Where no exact edge passes through the rounding cell of another
//! ring vertex, every vertex is the exact point rounded once to the nearest
//! binary64, as for [`HiddenLineSegment`](crate::HiddenLineSegment)
//! endpoints. Otherwise the edge is routed through that vertex, and a part
//! of a region whose area vanishes under rounding is dropped and counted in
//! the evidence, never refused. Rings are weakly simple: they may touch at a
//! vertex, never cross and never run along each other, and no ring vertex
//! lies inside the edge of any ring. A ring vertex is kept unless exactly two
//! boundary pieces of the snapped boundary meet there in a straight line, so
//! neighbouring regions share bit-identical vertices.

use axiolid_contracts::{
    Backend, CancellationGranularity, DisplacementBound, ExecutionOptions, ScratchRequirement,
};
use axiolid_core::{Aabb2, Point2};

use crate::{
    HiddenLineError, HiddenLineMesh, HiddenLineSource, MeshId, Occlusion, View, ViewProjection,
};

/// Caller-chosen id of a set of triangles that share one surface pattern.
///
/// A material id and a planar face are both encodings the caller picks; the
/// kernel infers planes only on request, as [`FaceGroups::Planar`].
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
    /// The planar patches of the mesh (ADR 0106). Two triangles are in one
    /// group when a chain of triangles joins them in which every two
    /// consecutive triangles share an edge that has exactly these two
    /// triangles and whose model-space normals are parallel and point the
    /// same way. Parallel and same-facing is decided exactly as
    /// `n1 × n2 = 0` in all three components and `n1 · n2 > 0`; no tolerance
    /// enters. A group's id is `FaceGroup(i)` for the smallest triangle
    /// index `i` in it, so it does not depend on the view or on the winding
    /// of the whole mesh. Providers resolve the groups before use.
    Planar,
}

impl FaceGroups<'_> {
    /// The group of triangle `triangle`.
    ///
    /// # Panics
    ///
    /// Panics for `PerTriangle` when `triangle` is out of range; providers
    /// check the length first. Panics for `Planar`, whose groups depend on
    /// the geometry; providers resolve it first.
    pub fn of(&self, triangle: usize) -> FaceGroup {
        match self {
            Self::Single => FaceGroup(0),
            Self::PerTriangle(groups) => groups[triangle],
            Self::Planar => panic!("planar face groups are resolved from the mesh"),
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
    /// The cap of a [`MeshRole::Solid`](crate::MeshRole::Solid) cut by the near plane: a surface at
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
    /// Maximum [`BehindSurface`] entries, summed over all output regions
    /// (ADR 0111). Zero refuses every result with an entry, so a request with
    /// see-through triangles needs it set.
    pub max_behind_surfaces: usize,
}

impl VisibleFaceLimits {
    /// Construct explicit input and output work limits, with no
    /// [`BehindSurface`] entry allowed.
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
            max_behind_surfaces: 0,
        }
    }

    /// Allow up to `max` [`BehindSurface`] entries over all output regions.
    pub const fn with_max_behind_surfaces(mut self, max: usize) -> Self {
        self.max_behind_surfaces = max;
        self
    }
}

/// One visible-face request.
///
/// All meshes occlude each other and themselves, except through their
/// see-through triangles, which list what lies behind them. A `region` crops
/// the image (closed box in image coordinates); geometry outside it neither
/// makes a
/// region nor, being outside the image, occludes anything inside.
#[non_exhaustive]
#[derive(Debug, Clone)]
pub struct VisibleFaceRequest<'a> {
    /// Input meshes; their order orders the output.
    pub meshes: &'a [HiddenLineMesh<'a>],
    /// The view.
    pub view: View,
    /// Optional image crop.
    pub region: Option<Aabb2>,
    /// Work limits.
    pub limits: VisibleFaceLimits,
}

impl<'a> VisibleFaceRequest<'a> {
    /// Find the visible faces of `meshes` in `view` without a crop.
    pub const fn new(
        meshes: &'a [HiddenLineMesh<'a>],
        view: View,
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

/// One surface seen through a see-through region (ADR 0111).
///
/// The derived equality compares all three fields; the order of a `behind`
/// list is by `layer`, then by the mesh's position in the request, then by
/// `source`, so the type has no `Ord` of its own.
#[non_exhaustive]
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub struct BehindSurface {
    /// The mesh the surface belongs to.
    pub mesh: MeshId,
    /// The group or cap of that mesh.
    pub source: FaceSource,
    /// The depth layer: 0 is the region's own depth, where sources tie with
    /// it, and each further layer is strictly behind the one before.
    pub layer: u32,
}

impl BehindSurface {
    /// Construct an entry for `source` of `mesh` in depth layer `layer`.
    pub const fn new(mesh: MeshId, source: FaceSource, layer: u32) -> Self {
        Self {
            mesh,
            source,
            layer,
        }
    }
}

/// A maximal region of the image where one source is frontmost and the same
/// surfaces are seen through it.
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
    /// Indices of the mesh's triangles that are frontmost somewhere in an
    /// exact region this region was snapped from, ascending. A region that
    /// snapping split at a pinch lists the whole set in each part. Empty for
    /// [`FaceSource::Cap`].
    pub triangles: Vec<u32>,
    /// What is seen through the region, front to back (ADR 0111): empty
    /// unless the source is see-through, and then the seen stack without the
    /// source's own layer-0 entry.
    pub behind: Vec<BehindSurface>,
}

impl VisibleRegion {
    /// Construct a region of `source` of `mesh` with nothing behind it.
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
            behind: Vec::new(),
        }
    }

    /// Set what is seen through the region.
    pub fn with_behind(mut self, behind: Vec<BehindSurface>) -> Self {
        self.behind = behind;
        self
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
    /// Bound, in the maximum norm of the image, on how far snap rounding
    /// moved region boundaries (ADR 0100): an image point farther than this
    /// from the exact boundary of a source lies in that source's regions
    /// exactly when it lies in its exact regions. Zero when every ring
    /// vertex is a binary64 point and no edge was rerouted.
    pub snap_displacement: DisplacementBound,
    /// Exact ring edges routed through the rounded point of another ring
    /// vertex whose rounding cell they pass through.
    pub rerouted_edges: usize,
    /// Exact regions whose area collapsed under snap rounding; no region is
    /// returned for them.
    pub collapsed_regions: usize,
    /// Holes of returned regions that collapsed under snap rounding.
    pub collapsed_holes: usize,
    /// [`Occlusion::SeeThrough`] flags in the request, summed over all meshes.
    pub see_through_triangles: usize,
    /// [`BehindSurface`] entries, summed over all output regions.
    pub behind_surfaces: usize,
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
            snap_displacement: DisplacementBound::ZERO,
            rerouted_edges: 0,
            collapsed_regions: 0,
            collapsed_holes: 0,
            see_through_triangles: 0,
            behind_surfaces: 0,
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
/// `source`, then by the outer ring's points, lexicographically, then by
/// `behind`, entry by entry, by `layer`, the mesh's position and `source`.
#[non_exhaustive]
#[derive(Debug, Clone, PartialEq)]
pub struct VisibleFaceOutcome {
    /// The view the regions are expressed in, unchanged.
    pub view: View,
    /// The regions.
    pub regions: Vec<VisibleRegion>,
    /// Counts.
    pub evidence: VisibleFaceEvidence,
}

impl VisibleFaceOutcome {
    /// Construct a provider result. Registries validate it before returning it.
    pub fn new(view: View, regions: Vec<VisibleRegion>, evidence: VisibleFaceEvidence) -> Self {
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
/// - a [`FaceSource::Cap`] region exists only for a [`MeshRole::Solid`](crate::MeshRole::Solid) that
///   `view.near` cuts, and has no `triangles`;
/// - regions of one source never overlap; regions of different sources
///   overlap only where they tie in depth, and a drawing that composites the
///   `behind` lists paints one of the tied regions only;
/// - flipping every triangle gives identical regions and identical evidence
///   except `exact_escalations`;
/// - the regions are the exact regions snap-rounded to the binary64 grid (ADR
///   0100), with the rings of [`VisibleRegion`] and the evidence of
///   [`VisibleFaceEvidence`]; a region is never refused because rounding
///   collapses it, and where no edge was rerouted the ring vertices are the
///   exact points rounded once to nearest;
/// - the tolerance of `options` applies to mesh audit and frame validation
///   only, never to a visibility decision;
/// - the same input refusals as [`MeshHiddenLines`](crate::MeshHiddenLines),
///   and a `PerTriangle` group slice whose length differs from the triangle
///   count is `InvalidInput` naming the mesh, as is a face group that holds
///   both opaque and see-through triangles; no substitute result is
///   returned;
/// - the [`BehindSurface`] entries are bounded by
///   [`VisibleFaceLimits::max_behind_surfaces`];
/// - [`FaceGroups::Planar`] is resolved by the provider as documented there,
///   exactly, with polled cancellation and a memory-budget refusal, and the
///   result equals the `PerTriangle` request with the resolved groups, except
///   `exact_escalations`;
/// - output storage is outside the scratch budget.
pub trait MeshVisibleFaces: Backend {
    /// Whether the provider computes regions for views of this projection.
    /// The default is `true` for [`ViewProjection::Orthographic`] only; see
    /// [`MeshHiddenLines::supports_projection`](crate::MeshHiddenLines::supports_projection)
    /// (ADR 0107).
    fn supports_projection(&self, projection: ViewProjection) -> bool {
        projection == ViewProjection::Orthographic
    }

    /// Whether the provider honours triangles of this occlusion. The default
    /// is `true` for [`Occlusion::Opaque`] only; see
    /// [`MeshHiddenLines::supports_occlusion`](crate::MeshHiddenLines::supports_occlusion)
    /// (ADR 0111). A provider that returns `true` for
    /// [`Occlusion::SeeThrough`] must pass `SeeThroughConformanceSuite::run_faces`.
    fn supports_occlusion(&self, occlusion: Occlusion) -> bool {
        occlusion == Occlusion::Opaque
    }

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
    use crate::MeshRole;
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
    fn with_planar_groups_replaces_the_groups() {
        let mesh = TriMesh::default();
        let m = HiddenLineMesh::new(MeshId(1), MeshRole::Sheet, &mesh).with_planar_groups();
        assert_eq!(m.groups, FaceGroups::Planar);
    }

    #[test]
    #[should_panic(expected = "planar face groups")]
    fn planar_groups_are_not_indexable() {
        let _ = FaceGroups::Planar.of(0);
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
    fn a_region_has_nothing_behind_it_unless_set() {
        let region = VisibleRegion::new(
            MeshId(2),
            FaceSource::Group(FaceGroup(0)),
            Vec::new(),
            Vec::new(),
            vec![0],
        );
        assert!(region.behind.is_empty());
        let entry = BehindSurface::new(MeshId(1), FaceSource::Cap, 2);
        assert_eq!(
            (entry.mesh, entry.source, entry.layer),
            (MeshId(1), FaceSource::Cap, 2)
        );
        let region = region.with_behind(vec![entry]);
        assert_eq!(region.behind, [entry]);
        assert_eq!(region.triangles, [0]);
    }

    #[test]
    fn limits_allow_no_behind_surface_unless_set() {
        let limits = VisibleFaceLimits::new(1, 2, 3, 4, 5, 6);
        assert_eq!(limits.max_behind_surfaces, 0);
        let limits = limits.with_max_behind_surfaces(7);
        assert_eq!(limits.max_behind_surfaces, 7);
        assert_eq!(limits.max_region_vertices, 6);
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
            e.rerouted_edges,
            e.collapsed_regions,
            e.collapsed_holes,
            e.see_through_triangles,
            e.behind_surfaces,
        ];
        assert!(counts.iter().all(|&c| c == 0));
        assert_eq!(e.snap_displacement, DisplacementBound::ZERO);
    }
}
