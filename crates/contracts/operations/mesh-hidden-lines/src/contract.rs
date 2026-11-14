//! Portable hidden-line provider contract: the views, the request and the
//! drawing a provider returns.

use axiolid_contracts::{
    Backend, CancellationGranularity, ExecutionOptions, GeomError, ScratchRequirement,
};
use axiolid_core::{Aabb2, Frame3, Point2};
use axiolid_mesh::TriMesh;

use crate::{FaceGroup, FaceGroups, HiddenLineError};

/// Caller-chosen identity of one input mesh, echoed in every segment and
/// refusal that concerns it. Ids must be unique within a request.
#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub struct MeshId(pub u64);

/// How an input mesh occludes and whether it is capped by a near plane.
#[non_exhaustive]
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum MeshRole {
    /// Closed, consistently wound two-manifold. A near clipping plane that
    /// cuts it (some vertex beyond the plane, `W > near`, and some not)
    /// leaves a cap that occludes: a point with `W < near`, farther from the
    /// viewer than the plane, is hidden by the cap when the point of the plane
    /// above it is inside the solid. A face of the solid in the near plane is
    /// no special case: it hides what lies below it itself.
    Solid,
    /// Edge-manifold with consistent winding, open or closed. Never capped.
    Sheet,
}

/// Whether a triangle hides what lies behind it (ADR 0111).
#[non_exhaustive]
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash, Default)]
pub enum Occlusion {
    /// The triangle hides what lies strictly behind it.
    #[default]
    Opaque,
    /// The triangle hides nothing. It is still drawn: its edges are selected
    /// and clipped like any other, a hidden-line piece behind it counts it in
    /// [`HiddenLineSegment::seen_through`], and a visible-face region of it
    /// lists what lies behind it in `VisibleRegion::behind`.
    SeeThrough,
}

/// The [`Occlusion`] of every triangle of one mesh.
#[non_exhaustive]
#[derive(Debug, Clone, Copy, PartialEq, Eq, Default)]
pub enum TriangleOcclusion<'a> {
    /// Every triangle is [`Occlusion::Opaque`].
    #[default]
    AllOpaque,
    /// One flag per triangle; the length equals the triangle count, which a
    /// provider checks.
    PerTriangle(&'a [Occlusion]),
}

impl TriangleOcclusion<'_> {
    /// Whether any triangle is [`Occlusion::SeeThrough`]. The length is not
    /// checked here.
    pub fn has_see_through(self) -> bool {
        self.see_through_count() > 0
    }

    /// The number of [`Occlusion::SeeThrough`] flags. The length is not
    /// checked here.
    pub fn see_through_count(self) -> usize {
        match self {
            Self::AllOpaque => 0,
            Self::PerTriangle(flags) => flags
                .iter()
                .filter(|&&flag| flag == Occlusion::SeeThrough)
                .count(),
        }
    }
}

/// One input mesh with its identity and role.
#[non_exhaustive]
#[derive(Debug, Clone, Copy)]
pub struct HiddenLineMesh<'a> {
    /// Identity echoed in the outcome.
    pub id: MeshId,
    /// Occlusion role.
    pub role: MeshRole,
    /// Geometry. Winding does not change the outcome.
    pub mesh: &'a TriMesh,
    /// Face groups, read by `MeshVisibleFaces` and, when
    /// [`EdgeSelection::group_boundary`] is set, by `MeshHiddenLines`.
    pub groups: FaceGroups<'a>,
    /// Which triangles hide what lies behind them (ADR 0111). With no
    /// [`Occlusion::SeeThrough`] flag the outcome is exactly that of
    /// [`TriangleOcclusion::AllOpaque`].
    pub occlusion: TriangleOcclusion<'a>,
}

impl<'a> HiddenLineMesh<'a> {
    /// Bind a mesh to an id and a role; every triangle is in
    /// [`FaceGroup`]`(0)` and opaque.
    pub const fn new(id: MeshId, role: MeshRole, mesh: &'a TriMesh) -> Self {
        Self {
            id,
            role,
            mesh,
            groups: FaceGroups::Single,
            occlusion: TriangleOcclusion::AllOpaque,
        }
    }

    /// Give every triangle its own [`Occlusion`]; `occlusion.len()` must equal
    /// the mesh's triangle count, which a provider checks.
    pub const fn with_occlusion(mut self, occlusion: &'a [Occlusion]) -> Self {
        self.occlusion = TriangleOcclusion::PerTriangle(occlusion);
        self
    }

    /// Give every triangle its own group id; `groups.len()` must equal the
    /// mesh's triangle count, which a provider checks.
    pub const fn with_groups(mut self, groups: &'a [FaceGroup]) -> Self {
        self.groups = FaceGroups::PerTriangle(groups);
        self
    }

    /// Group the triangles into planar patches ([`FaceGroups::Planar`]).
    pub const fn with_planar_groups(mut self) -> Self {
        self.groups = FaceGroups::Planar;
        self
    }
}

/// How a view maps model points to the image.
///
/// [`Orthographic`](Self::Orthographic) is the default. The in-tree hidden-line
/// providers draw both projections; the visible-face providers draw only the
/// orthographic one. A provider says which it draws with
/// `supports_projection` on [`MeshHiddenLines`] and [`crate::MeshVisibleFaces`]
/// (ADR 0107), and a registry never hands it a view it does not support.
#[non_exhaustive]
#[derive(Debug, Clone, Copy, PartialEq, Eq, Default)]
pub enum ViewProjection {
    /// Parallel projection along the frame's `-z` axis.
    #[default]
    Orthographic,
    /// Central projection from the eye at the frame's origin onto the plane
    /// `W = -1`: the image of a point with view coordinates `(U, V, W)` and
    /// `W < 0` is `(U / -W, V / -W)`. The view requires `near < 0`, so every
    /// kept point is strictly in front of the eye.
    Perspective,
}

/// A view: orthographic by default, or perspective (see [`ViewProjection`]).
///
/// For a model point `p` the view coordinates are `((p - o)·x, (p - o)·y)`
/// and the depth `W = (p - o)·z`, where `o`, `x`, `y`, `z` are the frame's
/// origin and axes. The viewer sits on the `+z` side and looks along `-z`, so
/// a larger `W` is nearer. An orthographic view's image is the view
/// coordinates; a perspective view's is described at
/// [`ViewProjection::Perspective`]. Geometry is kept where `far <= W <= near`
/// (both closed); an absent bound is no bound, except that a perspective view
/// has a `near` below zero.
#[non_exhaustive]
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct View {
    /// Right-handed orthonormal view frame.
    pub frame: Frame3,
    /// Nearest kept depth.
    pub near: Option<f64>,
    /// Farthest kept depth.
    pub far: Option<f64>,
    /// How points are projected to the image.
    pub projection: ViewProjection,
}

/// The former name of [`View`], from when every view was orthographic (ADR
/// 0107).
#[deprecated(note = "use View")]
pub type OrthographicView = View;

impl View {
    /// An orthographic view without depth clipping.
    pub const fn new(frame: Frame3) -> Self {
        Self {
            frame,
            near: None,
            far: None,
            projection: ViewProjection::Orthographic,
        }
    }

    /// A perspective view from the eye at `frame.origin`, keeping depths
    /// `W <= near`. `near` must be finite and below zero; anything else is
    /// refused, since the central projection is undefined at and behind the
    /// eye.
    pub fn perspective(frame: Frame3, near: f64) -> Result<Self, HiddenLineError> {
        if !near.is_finite() || near >= 0.0 {
            return Err(HiddenLineError::request(GeomError::InvalidInput(
                "perspective near depth must be finite and below zero".into(),
            )));
        }
        Ok(Self {
            frame,
            near: Some(near),
            far: None,
            projection: ViewProjection::Perspective,
        })
    }

    /// Keep only depths `W <= near`.
    pub const fn with_near(mut self, near: f64) -> Self {
        self.near = Some(near);
        self
    }

    /// Keep only depths `W >= far`.
    pub const fn with_far(mut self, far: f64) -> Self {
        self.far = Some(far);
        self
    }
}

/// Which edges to draw.
///
/// An edge with one incident triangle is a boundary edge. An edge with two
/// is a silhouette when the two triangles' image orientations differ (a
/// triangle seen edge-on is its own orientation, so the test is symmetric
/// and independent of winding), and otherwise a crease when the cosine of
/// the angle between their normals is strictly below `crease_cosine`, and
/// otherwise, when `group_boundary` is set, a group boundary when the two
/// triangles are in different [`FaceGroup`]s (ADR 0106). The kinds are
/// exclusive in that order: a silhouette is never a crease or a group
/// boundary, and an edge that is a crease at `crease_cosine` is never a group
/// boundary. With `group_boundary` unset the face groups are ignored.
#[non_exhaustive]
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct EdgeSelection {
    /// Crease threshold in `[-1, 1]`; `None` draws no creases.
    pub crease_cosine: Option<f64>,
    /// Draw boundary edges.
    pub boundary: bool,
    /// Draw silhouette edges.
    pub silhouette: bool,
    /// Draw edges between triangles of different face groups that are
    /// neither silhouettes nor selected creases. The groups are resolved as
    /// for `MeshVisibleFaces` ([`FaceGroups`]), and a `PerTriangle` slice of
    /// the wrong length is refused.
    pub group_boundary: bool,
    /// Also report the hidden pieces.
    pub hidden: bool,
}

impl EdgeSelection {
    /// Every edge kind, hidden pieces included.
    pub const fn new(crease_cosine: Option<f64>) -> Self {
        Self {
            crease_cosine,
            boundary: true,
            silhouette: true,
            group_boundary: false,
            hidden: true,
        }
    }

    /// Select or deselect boundary edges.
    pub const fn with_boundary(mut self, boundary: bool) -> Self {
        self.boundary = boundary;
        self
    }

    /// Select or deselect silhouette edges.
    pub const fn with_silhouette(mut self, silhouette: bool) -> Self {
        self.silhouette = silhouette;
        self
    }

    /// Select or deselect group-boundary edges.
    pub const fn with_group_boundary(mut self, group_boundary: bool) -> Self {
        self.group_boundary = group_boundary;
        self
    }

    /// Report or omit hidden pieces.
    pub const fn with_hidden(mut self, hidden: bool) -> Self {
        self.hidden = hidden;
        self
    }
}

/// Hard bounds for one hidden-line request.
#[non_exhaustive]
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct HiddenLineLimits {
    /// Maximum input meshes.
    pub max_meshes: usize,
    /// Maximum source positions, summed over all meshes.
    pub max_source_vertices: usize,
    /// Maximum source triangles, summed over all meshes.
    pub max_source_triangles: usize,
    /// Maximum (edge, triangle) pairs examined for occlusion.
    pub max_candidate_pairs: usize,
    /// Maximum output segments.
    pub max_output_segments: usize,
}

impl HiddenLineLimits {
    /// Construct explicit input and output work limits.
    pub const fn new(
        max_meshes: usize,
        max_source_vertices: usize,
        max_source_triangles: usize,
        max_candidate_pairs: usize,
        max_output_segments: usize,
    ) -> Self {
        Self {
            max_meshes,
            max_source_vertices,
            max_source_triangles,
            max_candidate_pairs,
            max_output_segments,
        }
    }
}

/// One hidden-line request.
///
/// All meshes occlude each other and themselves, except through their
/// [`Occlusion::SeeThrough`] triangles. A `region` crops the image
/// (closed box in image coordinates); geometry outside it neither is drawn
/// nor, being outside the image, occludes anything inside.
#[non_exhaustive]
#[derive(Debug, Clone)]
pub struct HiddenLineRequest<'a> {
    /// Input meshes; their order orders the output.
    pub meshes: &'a [HiddenLineMesh<'a>],
    /// The view.
    pub view: View,
    /// Edge kinds to draw.
    pub edges: EdgeSelection,
    /// Optional image crop.
    pub region: Option<Aabb2>,
    /// Work limits.
    pub limits: HiddenLineLimits,
}

impl<'a> HiddenLineRequest<'a> {
    /// Draw `meshes` in `view` without a crop.
    pub const fn new(
        meshes: &'a [HiddenLineMesh<'a>],
        view: View,
        edges: EdgeSelection,
        limits: HiddenLineLimits,
    ) -> Self {
        Self {
            meshes,
            view,
            edges,
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

/// Why an edge is drawn.
#[non_exhaustive]
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum EdgeKind {
    /// One incident triangle.
    Boundary,
    /// Two incident triangles seen with different image orientations.
    Silhouette,
    /// Two incident triangles meeting at a sharp angle.
    Crease,
    /// Two incident triangles of different face groups that are neither a
    /// silhouette nor a selected crease (ADR 0106).
    GroupBoundary,
}

/// Whether a piece of an edge can be seen.
#[non_exhaustive]
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Visibility {
    /// Nothing kept that occludes lies strictly in front of it. See-through
    /// triangles in front are counted in [`HiddenLineSegment::seen_through`].
    Visible,
    /// A kept opaque triangle point, or a cap, lies strictly in front of it.
    Hidden,
}

/// A source mesh edge.
#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub struct SourceEdge {
    /// The mesh the edge belongs to.
    pub mesh: MeshId,
    /// Position indices, ascending.
    pub vertices: [u32; 2],
}

/// A maximal positive-length piece of one source edge with one visibility
/// and one [`seen_through`](Self::seen_through) count.
///
/// Each endpoint is the exact point of the edge rounded once: every image
/// coordinate of the (rational) point is evaluated exactly and rounded to
/// the nearest binary64. An endpoint is therefore a function of the exact
/// point alone, whichever decision produced it, and equal points from
/// different decisions are bit-identical. All pieces of one source edge
/// have the same `kind`.
#[non_exhaustive]
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct HiddenLineSegment {
    /// The edge this piece belongs to.
    pub source: SourceEdge,
    /// Why the edge is drawn.
    pub kind: EdgeKind,
    /// Visibility of the whole piece.
    pub visibility: Visibility,
    /// Image point nearer `source.vertices[0]`.
    pub start: Point2,
    /// Image point nearer `source.vertices[1]`.
    pub end: Point2,
    /// For a visible piece, the number of kept [`Occlusion::SeeThrough`]
    /// triangles of nonzero image area that lie strictly in front of it
    /// (ADR 0111). A triangle counts at a point of the piece when its plane is
    /// strictly in front there and its image contains the point shifted by
    /// `(ε, ε²)` for an infinitesimal `ε > 0`, so two triangles sharing an
    /// edge above the piece count once, and a triangle at equal depth or seen
    /// edge-on never counts. Zero for a hidden piece.
    pub seen_through: u32,
}

impl HiddenLineSegment {
    /// Construct a piece of `source` with nothing seen through.
    pub const fn new(
        source: SourceEdge,
        kind: EdgeKind,
        visibility: Visibility,
        start: Point2,
        end: Point2,
    ) -> Self {
        Self {
            source,
            kind,
            visibility,
            start,
            end,
            seen_through: 0,
        }
    }

    /// Set the number of see-through triangles in front of the piece.
    pub const fn with_seen_through(mut self, seen_through: u32) -> Self {
        self.seen_through = seen_through;
        self
    }
}

/// Provenance of a hidden-line drawing.
#[non_exhaustive]
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum HiddenLineSource {
    /// Segments were derived from the supplied triangle meshes.
    InputMesh,
}

/// Auditable counts for one hidden-line operation.
#[non_exhaustive]
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct HiddenLineEvidence {
    /// Approximation source.
    pub source: HiddenLineSource,
    /// Input meshes.
    pub source_meshes: usize,
    /// Input triangles, summed over all meshes.
    pub source_triangles: usize,
    /// Edges of the selected kinds, before clipping.
    pub selected_edges: usize,
    /// Selected edges whose image is a single point. Such an edge is counted
    /// here only, before clipping: it is neither drawn nor counted in
    /// `clipped_edges`, whatever the depth range or the region.
    pub point_edges: usize,
    /// Selected edges removed entirely by depth clipping or the region.
    pub clipped_edges: usize,
    /// (edge, triangle) pairs examined for occlusion.
    pub candidate_pairs: usize,
    /// Decisions the fast filter could not settle and the exact tier did.
    /// This count may differ between two otherwise identical requests (for
    /// example with the triangles reversed); it is not part of the drawing.
    pub exact_escalations: usize,
    /// Pieces dropped because their rounded endpoints coincide.
    pub collapsed_segments: usize,
    /// Visible output segments.
    pub visible_segments: usize,
    /// Hidden output segments.
    pub hidden_segments: usize,
    /// [`Occlusion::SeeThrough`] flags in the request, summed over all meshes.
    pub see_through_triangles: usize,
    /// Output segments with a nonzero [`HiddenLineSegment::seen_through`].
    pub seen_through_segments: usize,
}

impl HiddenLineEvidence {
    /// Evidence for a drawing derived from the input meshes; the remaining
    /// counts start at zero for the provider to fill in.
    pub const fn input_mesh(source_meshes: usize, source_triangles: usize) -> Self {
        Self {
            source: HiddenLineSource::InputMesh,
            source_meshes,
            source_triangles,
            selected_edges: 0,
            point_edges: 0,
            clipped_edges: 0,
            candidate_pairs: 0,
            exact_escalations: 0,
            collapsed_segments: 0,
            visible_segments: 0,
            hidden_segments: 0,
            see_through_triangles: 0,
            seen_through_segments: 0,
        }
    }

    /// Whether this result came from the discrete input meshes.
    pub const fn is_derived_from_input_mesh(self) -> bool {
        matches!(self.source, HiddenLineSource::InputMesh)
    }
}

/// Image-space segments plus evidence.
///
/// Segments are ordered by the mesh's position in the request, then by
/// `source.vertices`, then by parameter along the edge.
#[derive(Debug, Clone, PartialEq)]
pub struct HiddenLineOutcome {
    /// The view the segments are expressed in, unchanged.
    pub view: View,
    /// Visible and hidden pieces.
    pub segments: Vec<HiddenLineSegment>,
    /// Counts.
    pub evidence: HiddenLineEvidence,
}

impl HiddenLineOutcome {
    /// Construct a provider result. Registries validate it before returning it.
    pub fn new(view: View, segments: Vec<HiddenLineSegment>, evidence: HiddenLineEvidence) -> Self {
        Self {
            view,
            segments,
            evidence,
        }
    }
}

/// Provider for deterministic hidden-line drawings of triangle meshes.
///
/// Provider obligations beyond the types:
///
/// - each segment is a maximal positive-length piece of one edge with one
///   visibility, finite, with `start != end`;
/// - a selected edge whose image is a point is counted in `point_edges`
///   and not drawn;
/// - with `edges.hidden == false` no hidden segment is returned;
/// - flipping every triangle gives identical segments and identical
///   evidence except `exact_escalations`;
/// - with `edges.group_boundary` the groups are resolved exactly as for
///   `MeshVisibleFaces`, and a `PerTriangle` group slice whose length differs
///   from the triangle count is `InvalidInput` naming the mesh; without it
///   the groups are not read;
/// - a [`TriangleOcclusion::PerTriangle`] slice whose length differs from
///   the triangle count is `InvalidInput` naming the mesh; a request whose
///   flags are all [`Occlusion::Opaque`] is drawn exactly as one without
///   flags, evidence included;
/// - endpoints are the exact points rounded once to nearest (see
///   [`HiddenLineSegment`]);
/// - the tolerance of `options` applies to mesh audit and frame validation
///   only, never to a visibility decision;
/// - input that cannot be decided is refused with a typed error naming the
///   mesh ([`HiddenLineError::mesh`]) where one mesh is at fault; no
///   substitute drawing is returned;
/// - output storage is outside the scratch budget.
pub trait MeshHiddenLines: Backend {
    /// Whether the provider draws views of this projection. The default is
    /// `true` for [`ViewProjection::Orthographic`] only, so a provider written
    /// before [`ViewProjection`] existed, which ignores the field, is never given
    /// a perspective view by a registry (ADR 0107). A provider that returns
    /// `true` for a projection must draw it exactly and pass that
    /// projection's conformance suite.
    fn supports_projection(&self, projection: ViewProjection) -> bool {
        projection == ViewProjection::Orthographic
    }

    /// Whether the provider honours triangles of this occlusion. The default
    /// is `true` for [`Occlusion::Opaque`] only, so a provider written before
    /// [`Occlusion`] existed, which would let every triangle hide, is never
    /// given a see-through triangle by a registry (ADR 0111). A provider that
    /// returns `true` for [`Occlusion::SeeThrough`] must pass
    /// `SeeThroughConformanceSuite::run`.
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

    /// Draw the selected edges of `request.meshes` in `request.view`.
    /// Implementations must enforce `request.limits` before growing output.
    fn hidden_lines(
        &self,
        request: &HiddenLineRequest<'_>,
        options: &ExecutionOptions,
    ) -> Result<HiddenLineOutcome, HiddenLineError>;
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::fixtures::plan_frame;

    #[test]
    fn a_view_is_orthographic_unless_built_as_perspective() {
        let view = View::new(plan_frame());
        assert_eq!(view.projection, ViewProjection::Orthographic);
        assert_eq!(ViewProjection::default(), ViewProjection::Orthographic);
        assert_eq!(
            view.with_near(-1.0).projection,
            ViewProjection::Orthographic
        );
    }

    #[test]
    fn occlusion_defaults_to_opaque_for_every_triangle() {
        let mesh = TriMesh::default();
        let m = HiddenLineMesh::new(MeshId(1), MeshRole::Sheet, &mesh);
        assert_eq!(m.occlusion, TriangleOcclusion::AllOpaque);
        assert_eq!(TriangleOcclusion::default(), TriangleOcclusion::AllOpaque);
        assert_eq!(Occlusion::default(), Occlusion::Opaque);
        assert!(!m.occlusion.has_see_through());
        assert_eq!(m.occlusion.see_through_count(), 0);
        let flags = [
            Occlusion::Opaque,
            Occlusion::SeeThrough,
            Occlusion::SeeThrough,
        ];
        let m = m.with_occlusion(&flags);
        assert_eq!(m.occlusion, TriangleOcclusion::PerTriangle(&flags));
        assert!(m.occlusion.has_see_through());
        assert_eq!(m.occlusion.see_through_count(), 2);
        let opaque = [Occlusion::Opaque; 4];
        assert!(!TriangleOcclusion::PerTriangle(&opaque).has_see_through());
        assert!(!TriangleOcclusion::PerTriangle(&[]).has_see_through());
    }

    #[test]
    fn a_segment_sees_through_nothing_unless_set() {
        let source = SourceEdge {
            mesh: MeshId(1),
            vertices: [0, 1],
        };
        let p = Point2::new(0.0, 0.0);
        let q = Point2::new(1.0, 0.0);
        let s = HiddenLineSegment::new(source, EdgeKind::Boundary, Visibility::Visible, p, q);
        assert_eq!(s.seen_through, 0);
        let s = s.with_seen_through(3);
        assert_eq!(s.seen_through, 3);
        assert_eq!((s.start, s.end, s.visibility), (p, q, Visibility::Visible));
    }

    #[test]
    fn evidence_starts_without_see_through_counts() {
        let e = HiddenLineEvidence::input_mesh(2, 12);
        assert_eq!((e.see_through_triangles, e.seen_through_segments), (0, 0));
    }

    #[test]
    #[allow(deprecated)]
    fn the_old_name_is_the_same_type() {
        let old: crate::OrthographicView = crate::OrthographicView::new(plan_frame());
        let new: View = old;
        assert_eq!(new, View::new(plan_frame()));
    }

    #[test]
    fn a_perspective_view_needs_a_finite_negative_near() {
        let view = View::perspective(plan_frame(), -1.0).expect("valid");
        assert_eq!(view.projection, ViewProjection::Perspective);
        assert_eq!(view.near, Some(-1.0));
        assert_eq!(view.far, None);
        for near in [0.0, -0.0, 1.0, f64::NAN, f64::INFINITY, f64::NEG_INFINITY] {
            let error = View::perspective(plan_frame(), near).expect_err("refused");
            assert_eq!(error.mesh, None);
            assert!(matches!(error.error, GeomError::InvalidInput(_)), "{near}");
        }
    }
}
