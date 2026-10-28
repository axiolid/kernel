//! Portable orthographic hidden-line provider contract.

use axiolid_contracts::{Backend, CancellationGranularity, ExecutionOptions, ScratchRequirement};
use axiolid_core::{Aabb2, Frame3, Point2};
use axiolid_mesh::TriMesh;

use crate::HiddenLineError;

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
}

impl<'a> HiddenLineMesh<'a> {
    /// Bind a mesh to an id and a role.
    pub const fn new(id: MeshId, role: MeshRole, mesh: &'a TriMesh) -> Self {
        Self { id, role, mesh }
    }
}

/// An orthographic view.
///
/// For a model point `p` the image is `((p - o)·x, (p - o)·y)` and the depth
/// is `W = (p - o)·z`, where `o`, `x`, `y`, `z` are the frame's origin and
/// axes. The viewer sits on the `+z` side and looks along `-z`, so a larger
/// `W` is nearer. Geometry is kept where `far <= W <= near` (both closed);
/// an absent bound is no bound.
#[non_exhaustive]
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct OrthographicView {
    /// Right-handed orthonormal view frame.
    pub frame: Frame3,
    /// Nearest kept depth.
    pub near: Option<f64>,
    /// Farthest kept depth.
    pub far: Option<f64>,
}

impl OrthographicView {
    /// A view without depth clipping.
    pub const fn new(frame: Frame3) -> Self {
        Self {
            frame,
            near: None,
            far: None,
        }
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
/// the angle between their normals is strictly below `crease_cosine`.
#[non_exhaustive]
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct EdgeSelection {
    /// Crease threshold in `[-1, 1]`; `None` draws no creases.
    pub crease_cosine: Option<f64>,
    /// Draw boundary edges.
    pub boundary: bool,
    /// Draw silhouette edges.
    pub silhouette: bool,
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
/// All meshes occlude each other and themselves. A `region` crops the image
/// (closed box in image coordinates); geometry outside it neither is drawn
/// nor, being outside the image, occludes anything inside.
#[non_exhaustive]
#[derive(Debug, Clone)]
pub struct HiddenLineRequest<'a> {
    /// Input meshes; their order orders the output.
    pub meshes: &'a [HiddenLineMesh<'a>],
    /// The orthographic view.
    pub view: OrthographicView,
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
        view: OrthographicView,
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
}

/// Whether a piece of an edge can be seen.
#[non_exhaustive]
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Visibility {
    /// Nothing kept lies strictly in front of it.
    Visible,
    /// A kept triangle point lies strictly in front of it.
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

/// A maximal positive-length piece of one source edge with one visibility.
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
}

impl HiddenLineSegment {
    /// Construct a piece of `source`.
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
        }
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
    pub view: OrthographicView,
    /// Visible and hidden pieces.
    pub segments: Vec<HiddenLineSegment>,
    /// Counts.
    pub evidence: HiddenLineEvidence,
}

impl HiddenLineOutcome {
    /// Construct a provider result. Registries validate it before returning it.
    pub fn new(
        view: OrthographicView,
        segments: Vec<HiddenLineSegment>,
        evidence: HiddenLineEvidence,
    ) -> Self {
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
/// - endpoints are the exact points rounded once to nearest (see
///   [`HiddenLineSegment`]);
/// - the tolerance of `options` applies to mesh audit and frame validation
///   only, never to a visibility decision;
/// - input that cannot be decided is refused with a typed error naming the
///   mesh ([`HiddenLineError::mesh`]) where one mesh is at fault; no
///   substitute drawing is returned;
/// - output storage is outside the scratch budget.
pub trait MeshHiddenLines: Backend {
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
