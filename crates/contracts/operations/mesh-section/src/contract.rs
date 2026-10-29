//! Portable mesh plane-section provider contract.

use axiolid_contracts::{
    Backend, CancellationGranularity, ExecutionOptions, GeomError, GeomResult, Operation,
    ScratchRequirement,
};
use axiolid_core::{Frame3, Point2};
use axiolid_mesh::TriMesh;

use crate::request::{PlaneSection, SectionRequest};

/// Hard bounds for one mesh section request.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct SectionLimits {
    /// Maximum source positions inspected.
    pub max_source_vertices: usize,
    /// Maximum source triangles inspected.
    pub max_source_triangles: usize,
    /// Maximum vertices across all output contours.
    pub max_output_vertices: usize,
    /// Maximum output contours.
    pub max_contours: usize,
}

impl SectionLimits {
    /// Construct explicit source and output work limits.
    pub const fn new(
        max_source_vertices: usize,
        max_source_triangles: usize,
        max_output_vertices: usize,
        max_contours: usize,
    ) -> Self {
        Self {
            max_source_vertices,
            max_source_triangles,
            max_output_vertices,
            max_contours,
        }
    }
}

/// One closed plane-local polyline.
///
/// The terminal point is implicit and must not duplicate the first point.
#[derive(Debug, Clone, PartialEq)]
pub struct SectionContour {
    /// Plane-local points in traversal order.
    pub points: Vec<Point2>,
    /// Private invariant marker: every constructible contour is closed.
    closed: bool,
}

impl SectionContour {
    /// Construct a closed contour. Registry validation checks cardinality and
    /// finite, non-duplicated coordinates after provider dispatch.
    pub fn new(points: Vec<Point2>) -> Self {
        Self {
            points,
            closed: true,
        }
    }

    /// Whether this polyline closes from its last point back to its first.
    pub const fn is_closed(&self) -> bool {
        self.closed
    }
}

/// Provenance of a plane-section approximation.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
#[non_exhaustive]
pub enum SectionSource {
    /// Contours were intersected with the supplied triangle mesh.
    InputMesh,
}

/// Auditable counts for one section operation.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct SectionEvidence {
    /// Approximation source.
    pub source: SectionSource,
    /// Source triangles inspected.
    pub source_triangles: usize,
    /// Output contour vertices.
    pub output_vertices: usize,
    /// Output contours.
    pub output_contours: usize,
}

impl SectionEvidence {
    /// Record evidence for a section derived from the input mesh.
    pub const fn input_mesh(
        source_triangles: usize,
        output_vertices: usize,
        output_contours: usize,
    ) -> Self {
        Self {
            source: SectionSource::InputMesh,
            source_triangles,
            output_vertices,
            output_contours,
        }
    }

    /// Whether this result came from the discrete input mesh.
    pub const fn is_derived_from_input_mesh(self) -> bool {
        matches!(self.source, SectionSource::InputMesh)
    }
}

/// Plane-local closed section contours plus approximation evidence.
#[derive(Debug, Clone, PartialEq)]
pub struct SectionOutcome {
    /// Right-handed orthonormal mapping from local `(x,y,0)` to model space.
    pub frame: Frame3,
    /// Closed section contours. Empty means the plane misses the solid.
    pub contours: Vec<SectionContour>,
    /// Source and output counts.
    pub evidence: SectionEvidence,
}

impl SectionOutcome {
    /// Construct a provider result. Registries validate it before returning it.
    pub fn new(frame: Frame3, contours: Vec<SectionContour>, evidence: SectionEvidence) -> Self {
        Self {
            frame,
            contours,
            evidence,
        }
    }
}

/// Provider for deterministic sections of closed oriented triangle solids.
pub trait MeshPlaneSection: Backend {
    /// Scratch needed beyond inputs and result.
    fn scratch_requirement(&self) -> ScratchRequirement {
        ScratchRequirement::Unbounded
    }

    /// How finely the provider polls cancellation.
    fn cancellation_granularity(&self) -> CancellationGranularity {
        CancellationGranularity::None
    }

    /// Intersect one validated closed oriented mesh with `frame`'s local XY
    /// plane. Implementations must enforce `limits` before growing output.
    fn section(
        &self,
        mesh: &TriMesh,
        frame: Frame3,
        limits: SectionLimits,
        options: &ExecutionOptions,
    ) -> GeomResult<SectionOutcome>;

    /// Scratch needed beyond inputs and result by
    /// [`section_request`](Self::section_request). Defaults to
    /// [`scratch_requirement`](Self::scratch_requirement).
    fn request_scratch_requirement(&self) -> ScratchRequirement {
        self.scratch_requirement()
    }

    /// Intersect one validated mesh with the request frame's local XY plane,
    /// reporting faces lying in the plane as regions (ADR 0091). Providers
    /// that do not implement it refuse with `Unsupported`, which lets the
    /// registry fall back to the next provider.
    ///
    /// A provider may implement `SectionRole::Solid` only: it then refuses a
    /// `SectionRole::Sheet` request with `Unsupported`, and the region
    /// conformance suite reports `sheets_claimed == false` instead of failing.
    ///
    /// A provider that does not trace cut regions refuses a request that has
    /// them with `Unsupported`; the suite then reports
    /// `cut_regions_claimed == false`. A provider written before cut regions
    /// existed ignores the request and answers with `cut_regions == None`: the
    /// registry treats that answer as `Unsupported` and tries the next
    /// provider, and the suite treats it on the first solid case as not
    /// claimed. A provider that returns cut regions on one case and drops them
    /// on another fails the suite.
    fn section_request(
        &self,
        mesh: &TriMesh,
        request: &SectionRequest,
        options: &ExecutionOptions,
    ) -> GeomResult<PlaneSection> {
        let _ = (mesh, request, options);
        Err(GeomError::Unsupported {
            backend: self.descriptor().id,
            operation: Operation::MeshPlaneSection,
        })
    }
}
