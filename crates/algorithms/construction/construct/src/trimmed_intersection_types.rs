use axiolid_brep::{Curve2Id, ExactBRep};
use axiolid_core::{Interval, Scalar};
use axiolid_nurbs::{
    CertifiedSurfaceSurfaceIntersection3, CertifiedSurfaceSurfaceIntersectionOptions,
    TransverseSurfaceSurfaceTrace3,
};
use axiolid_topology::{EdgeId, FaceId};

/// Explicit resource and residual policy for certified topology integration.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct CertifiedSurfacePairSplitOptions {
    /// Policy for the underlying certified surface/surface query.
    intersection: CertifiedSurfaceSurfaceIntersectionOptions,
    /// Largest certified surface/carrier residual accepted by the B-rep handoff.
    max_surface_residual: Scalar,
}

impl CertifiedSurfacePairSplitOptions {
    /// Validate and construct a split policy.
    pub fn new(
        intersection: CertifiedSurfaceSurfaceIntersectionOptions,
        max_surface_residual: Scalar,
    ) -> Result<Self, axiolid_contracts::GeomError> {
        if !max_surface_residual.is_finite() || max_surface_residual <= 0.0 {
            return Err(axiolid_contracts::GeomError::InvalidInput(
                "surface-pair split residual must be finite and positive".into(),
            ));
        }
        Ok(Self {
            intersection,
            max_surface_residual,
        })
    }

    pub(super) fn intersection_options(self) -> CertifiedSurfaceSurfaceIntersectionOptions {
        self.intersection
    }

    pub(super) fn max_surface_residual(self) -> Scalar {
        self.max_surface_residual
    }
}

impl Default for CertifiedSurfacePairSplitOptions {
    fn default() -> Self {
        Self {
            intersection: CertifiedSurfaceSurfaceIntersectionOptions::default(),
            max_surface_residual: 1.0e-7,
        }
    }
}

/// Which input surface is partitioned by the certified chord.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum SurfacePairMember {
    /// The first function argument.
    First,
    /// The second function argument.
    Second,
}

/// A certified edge embedded in a face without pretending it is a closed trim.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct EmbeddedFaceCurve {
    /// Unsplit face that contains the intersection chord.
    pub face: FaceId,
    /// Shared model-space intersection edge.
    pub edge: EdgeId,
    /// Surface-parameter image from edge start to edge end.
    pub pcurve: Curve2Id,
    /// Native pcurve interval, oriented from edge start to edge end.
    pub interval: Interval,
}

/// Validated analytic B-rep arrangement for one certified finite trace.
#[derive(Debug, Clone, PartialEq)]
pub struct CertifiedTrimmedSurfacePair3 {
    /// Strict analytic B-rep containing three closed faces.
    pub brep: ExactBRep,
    /// The same edge used by the two split-face loops and embedded in the other face.
    pub intersection_edge: EdgeId,
    /// Input member whose rectangular face was split.
    pub split_surface: SurfacePairMember,
    /// Two deterministic trimmed faces on `split_surface`.
    pub split_faces: [FaceId; 2],
    /// Rectangular face that the finite chord does not partition.
    pub unsplit_face: FaceId,
    /// Explicit interior attachment on `unsplit_face`.
    pub embedded_curve: EmbeddedFaceCurve,
    /// Original certified bounded trace; no endpoint is widened by construction.
    pub trace: TransverseSurfaceSurfaceTrace3,
    /// Conservative global carrier-to-surface residual bound.
    pub max_surface_residual_upper_bound: Scalar,
    /// Certified patch pairs processed by the intersection query.
    pub visited_patch_pairs: u32,
    /// Bounded boundary queries used by the intersection query.
    pub boundary_queries: u32,
}

/// Why valid geometry could not be promoted to a closed trimmed arrangement.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum SurfacePairSplitUnresolvedReason {
    /// The underlying certified query retained uncertainty.
    IntersectionUnresolved,
    /// The query did not produce exactly one finite trace.
    UnsupportedTraceCount,
    /// Endpoint ownership did not identify exactly one partitioned rectangle.
    ///
    /// The arrangement may still be constructible; this crate has not
    /// implemented it. Distinct from [`Self::NoPartitionExists`], which is a
    /// proof that no arrangement exists at all.
    UnsupportedEndpointOwnership,
    /// Proven: this trace cannot partition either patch, so no split exists.
    ///
    /// A trace that stops strictly inside a patch is a slit, not a partition:
    /// the face stays simply connected, exactly as a slot cut halfway into a
    /// sheet leaves one piece. When that holds on *both* patches, neither can
    /// become two closed trimmed faces and retrying with different options or
    /// tighter tolerances cannot change the answer.
    ///
    /// This is terminal. Callers should stop rather than escalate.
    NoPartitionExists,
    /// The certified carrier residual exceeded explicit policy.
    ResidualExceedsPolicy,
    /// Representative parameters or carrier were finite but degenerate for topology.
    DegenerateRepresentative,
}

/// Certified topology-integration outcome.
///
/// The success payload stays inline deliberately: boxing it would add an
/// infallible allocation after all certified construction allocations have
/// already been made fallible and bounded.
#[allow(clippy::large_enum_variant)]
#[derive(Debug, Clone, PartialEq)]
#[non_exhaustive]
pub enum CertifiedSurfacePairSplit3 {
    /// Conservative proof that the bounded patches do not intersect.
    Empty {
        /// Certified patch pairs processed.
        visited_patch_pairs: u32,
        /// Bounded boundary queries used.
        boundary_queries: u32,
    },
    /// One finite trace integrated into a strict trimmed arrangement.
    Split(CertifiedTrimmedSurfacePair3),
    /// Geometry remains usable, but no B-rep split was invented.
    Unresolved {
        /// Original intersection evidence retained without widening.
        intersection: CertifiedSurfaceSurfaceIntersection3,
        /// Topology-specific refusal reason.
        reason: SurfacePairSplitUnresolvedReason,
    },
}
