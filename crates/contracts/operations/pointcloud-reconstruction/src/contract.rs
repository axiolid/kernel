//! The provider trait: what implementing reconstruction obliges you to.

use axiolid_contracts::{
    Backend, CancellationGranularity, Determinism, ExecutionOptions, GeomResult, ScratchRequirement,
};
use axiolid_core::Scalar;
use axiolid_pointcloud::PointCloud;

use crate::{ReconstructionOutcome, ReconstructionRequest};

/// Why a reconstruction could not be performed.
///
/// Every variant is a *typed refusal*: the provider knows why it cannot
/// answer and says so. A provider must never return an empty mesh in place
/// of one of these — an empty result and "I cannot do this" are different
/// answers, and a caller that cannot tell them apart will treat a failure as
/// an object that is not there.
#[derive(Debug, Clone, PartialEq)]
#[non_exhaustive]
pub enum ReconstructionRefusal {
    /// Fewer points than the method needs.
    TooFewPoints {
        /// Points supplied.
        supplied: usize,
        /// Minimum this provider needs.
        required: usize,
    },
    /// The points are degenerate: all coincident, collinear, or coplanar.
    ///
    /// A surface cannot be reconstructed from a set with no volume, and
    /// producing a zero-thickness sheet would misrepresent the input.
    DegenerateExtent {
        /// What the provider measured.
        detail: String,
    },
    /// The request demanded a closed result the data cannot support.
    ///
    /// Raised instead of inventing surface across a gap the capture never
    /// covered.
    CannotClose {
        /// Why closure is impossible for this input.
        detail: String,
    },
    /// The requested resolution is finer than the samples resolve.
    ///
    /// Honouring it would interpolate detail that was never measured and
    /// present it as geometry.
    ResolutionExceedsData {
        /// Edge length the caller asked for.
        requested: Scalar,
        /// Spacing the samples actually resolve.
        sample_spacing: Scalar,
    },
    /// The provider does not implement some part of the request.
    Unsupported {
        /// What is unsupported.
        detail: String,
    },
}

impl core::fmt::Display for ReconstructionRefusal {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        match self {
            Self::TooFewPoints { supplied, required } => write!(
                f,
                "reconstruction needs at least {required} points, got {supplied}"
            ),
            Self::DegenerateExtent { detail } => {
                write!(f, "point set has no reconstructable extent: {detail}")
            }
            Self::CannotClose { detail } => {
                write!(f, "a closed result was required but is not supported by the data: {detail}")
            }
            Self::ResolutionExceedsData {
                requested,
                sample_spacing,
            } => write!(
                f,
                "requested edge length {requested} is finer than the sample spacing {sample_spacing}"
            ),
            Self::Unsupported { detail } => write!(f, "unsupported request: {detail}"),
        }
    }
}

impl core::error::Error for ReconstructionRefusal {}

/// Either a reconstruction or a typed reason there is none.
///
/// Deliberately not `Option<TriMesh>`: an absent surface always has a
/// reason, and this type makes the reason impossible to drop.
#[derive(Debug, Clone, PartialEq)]
pub enum Reconstruction {
    /// A surface was produced.
    Surface(Box<ReconstructionOutcome>),
    /// No surface was produced, for this stated reason.
    Refused(ReconstructionRefusal),
}

impl Reconstruction {
    /// The outcome, if one was produced.
    pub fn outcome(&self) -> Option<&ReconstructionOutcome> {
        match self {
            Self::Surface(outcome) => Some(outcome),
            Self::Refused(_) => None,
        }
    }

    /// The refusal, if there was one.
    pub fn refusal(&self) -> Option<&ReconstructionRefusal> {
        match self {
            Self::Surface(_) => None,
            Self::Refused(reason) => Some(reason),
        }
    }

    /// Whether a surface was produced.
    pub fn is_surface(&self) -> bool {
        matches!(self, Self::Surface(_))
    }
}

/// Pointcloud-to-surface reconstruction provider.
///
/// Implementing this trait is the capability declaration: a type that
/// implements it claims it can turn an unstructured point set into a
/// surface. Providers that cannot must not implement it.
pub trait PointcloudReconstruction: Backend {
    /// Scratch this provider needs beyond its inputs and result.
    ///
    /// Defaults to [`ScratchRequirement::Unbounded`] so an unaudited
    /// provider is treated as unbudgetable rather than assumed cheap.
    /// Reconstruction is memory-hungry — an octree over a hundred million
    /// points is not a detail a caller should discover by being killed.
    fn scratch_requirement(&self) -> ScratchRequirement {
        ScratchRequirement::Unbounded
    }

    /// Reproducibility this provider guarantees.
    ///
    /// Defaults to the weakest level so a provider that has not audited
    /// itself cannot silently satisfy a stronger request.
    fn determinism(&self) -> Determinism {
        Determinism::BestEffort
    }

    /// How finely this provider observes cancellation.
    ///
    /// Defaults to [`CancellationGranularity::None`]: a provider that has
    /// not audited its own polling must not advertise responsiveness it
    /// cannot deliver. Reconstruction runs long, so overstating this is
    /// exactly the case a caller would feel.
    fn cancellation_granularity(&self) -> CancellationGranularity {
        CancellationGranularity::None
    }

    /// Minimum points this provider needs to attempt a reconstruction.
    ///
    /// Callers can check before paying for a large request. A provider
    /// still refuses by name if given fewer.
    fn minimum_points(&self) -> usize {
        4
    }

    /// Reconstruct a surface from a point set.
    ///
    /// Returns [`Reconstruction::Refused`] with a typed reason when the
    /// request cannot be honoured. An `Err` is reserved for a failure of
    /// the machinery itself — a refusal is an answer, not an error.
    fn reconstruct(
        &self,
        cloud: &PointCloud,
        request: &ReconstructionRequest,
        options: &ExecutionOptions,
    ) -> GeomResult<Reconstruction>;
}
