//! What a reconstruction was asked to do, and what it actually did.
//!
//! A reconstruction is an *estimate*: a point set does not determine a unique
//! surface, so every result carries the assumptions that produced it. The
//! request records what the caller asked for; the evidence records what the
//! provider did, including where it had to guess.

use axiolid_core::Scalar;
use axiolid_mesh::TriMesh;

/// How densely the surface should be reconstructed.
///
/// Not an algorithm choice: a caller states the resolution it needs in model
/// units and the provider maps that onto whatever its method uses (octree
/// depth, ball radius, alpha value). Stating it in units keeps the request
/// portable across providers that share no parameters.
#[derive(Debug, Clone, Copy, PartialEq)]
pub enum Resolution {
    /// Target edge length of the reconstructed surface, in model units.
    ///
    /// The provider is expected to approach this, not to guarantee it: the
    /// achieved value is reported in [`ReconstructionEvidence`].
    TargetEdgeLength(Scalar),
    /// Let the provider choose from the sample spacing it measures.
    ///
    /// Honest default for a caller with no resolution requirement. What was
    /// chosen is reported rather than left implicit.
    FromSampleSpacing,
}

/// What the caller wants reconstructed.
#[derive(Debug, Clone, PartialEq)]
#[non_exhaustive]
pub struct ReconstructionRequest {
    /// Surface density to aim for.
    pub resolution: Resolution,
    /// Whether the result must be a closed solid.
    ///
    /// A scan of a room's interior has no back faces, so demanding closure
    /// forces the provider to invent surface where it has no data. Set this
    /// only when the capture genuinely covers the whole object; a provider
    /// that cannot honour it must refuse rather than fabricate.
    pub require_closed: bool,
    /// Whether per-point normals may be used when present.
    ///
    /// Normals dramatically improve most methods, but a capture's normals
    /// can be wrong. Turning this off asks the provider to work from
    /// positions alone.
    pub use_normals: bool,
}

impl Default for ReconstructionRequest {
    fn default() -> Self {
        Self {
            resolution: Resolution::FromSampleSpacing,
            require_closed: false,
            use_normals: true,
        }
    }
}

impl ReconstructionRequest {
    /// A request at a stated resolution.
    pub fn at_edge_length(edge_length: Scalar) -> Self {
        Self {
            resolution: Resolution::TargetEdgeLength(edge_length),
            ..Self::default()
        }
    }

    /// Demand a closed result, refusing rather than inventing surface.
    pub fn requiring_closed(mut self) -> Self {
        self.require_closed = true;
        self
    }

    /// Ignore per-point normals even when the cloud carries them.
    pub fn ignoring_normals(mut self) -> Self {
        self.use_normals = false;
        self
    }
}

/// Counters describing one reconstruction.
///
/// Every field is a fact about the computation, never a quality verdict:
/// whether 12% unsupported area is acceptable is the caller's decision, not
/// the kernel's.
#[derive(Debug, Clone, PartialEq, Default)]
#[non_exhaustive]
pub struct ReconstructionEvidence {
    /// Points supplied.
    pub input_points: usize,
    /// Points the provider actually used.
    ///
    /// Lower than `input_points` when the provider discarded outliers or
    /// duplicates. A large gap is worth a caller's attention.
    pub used_points: usize,
    /// Triangles in the result.
    pub output_triangles: usize,
    /// Connected components in the result.
    ///
    /// A capture with gaps often reconstructs into several shells. Reporting
    /// it lets a caller detect the split here rather than downstream.
    pub output_components: usize,
    /// Whether the result is a closed two-manifold solid.
    pub closed: bool,
    /// Edge length the provider achieved, in model units.
    ///
    /// Reported whether or not the caller stated one, so a
    /// `FromSampleSpacing` request still learns what it got.
    pub achieved_edge_length: Scalar,
    /// Median distance between neighbouring input samples.
    ///
    /// The scale the capture actually resolves. A result finer than this is
    /// interpolation, not measurement.
    pub sample_spacing: Scalar,
    /// Whether per-point normals were used.
    ///
    /// False when the cloud carried none, or the request declined them, or
    /// the provider could not use them. Distinguishes a normal-driven
    /// result from a positions-only one.
    pub used_normals: bool,
    /// Triangles whose vertices came from interpolation across a gap in the
    /// data rather than from measured samples.
    ///
    /// This is the honest measure of how much of the surface is invented.
    /// Zero means every triangle rests on real points.
    pub interpolated_triangles: usize,
}

impl ReconstructionEvidence {
    /// Start from a cleared record and fill in what the provider measured.
    ///
    /// `#[non_exhaustive]` keeps the struct additive, so a provider outside
    /// this crate cannot build one with a literal. This constructor is the
    /// supported route: take the default and set the fields you measured.
    /// A field left at its default is honestly "not measured" rather than a
    /// value invented to satisfy the type.
    pub fn measured() -> Self {
        Self::default()
    }
}

/// A reconstruction and the evidence that produced it.
#[derive(Debug, Clone, PartialEq)]
pub struct ReconstructionOutcome {
    /// The reconstructed surface. An empty mesh is a legitimate value only
    /// when the provider also explains it in the evidence.
    pub mesh: TriMesh,
    /// What the provider did.
    pub evidence: ReconstructionEvidence,
}

impl ReconstructionOutcome {
    /// Pair a mesh with its evidence.
    pub const fn new(mesh: TriMesh, evidence: ReconstructionEvidence) -> Self {
        Self { mesh, evidence }
    }
}
