//! Opt-in section request: faces lying in the cutting plane become regions.
//!
//! The legacy [`MeshPlaneSection::section`](crate::MeshPlaneSection::section)
//! refuses a mesh with a face in the cutting plane. This module adds the
//! request, the result and the refusal names for
//! [`MeshPlaneSection::section_request`](crate::MeshPlaneSection::section_request),
//! which reports such faces as [`InPlaneRegion`]s tagged with the side their
//! material lies on, kept apart from the cut curves (ADR 0091, #359).
//!
//! # Cut curves
//!
//! A [`SectionChain`] contains transversal crossings only. A mesh edge lying
//! in the plane yields a cut segment only when both of its triangles are off
//! the plane and their third vertices lie on opposite sides. An edge next to
//! an in-plane triangle is a region boundary, and an edge whose neighbours
//! all lie on one side only touches the plane; it is counted, not emitted.
//! A chain can therefore be open and end on a region boundary. Chain end
//! points and region vertices are projected by the same function, so they are
//! bit-identical.
//!
//! # Output invariants
//!
//! All points are finite and in the plane frame. A closed chain has at least
//! three points, runs counter-clockwise, starts at its lexicographically
//! smallest point (`x`, then `y`, by `total_cmp`) and has no adjacent
//! duplicates, wrap-around included. An open chain has at least two points,
//! no adjacent duplicates, and `points[0] < points[last]`. Chains are sorted
//! by their point sequences, a closed chain before an open one at a tie.
//! A region's outer ring is counter-clockwise and its holes are clockwise;
//! every ring has at least three points and starts at its lexicographically
//! smallest point. Holes are sorted by point sequence, regions by
//! `(side, outer)`. Rings may touch at vertices but never cross, and regions
//! of one side have disjoint interiors. Where a hole or a footprint touches the
//! outer boundary at a vertex, the region has one weakly simple outer ring that
//! visits that vertex twice, and no hole. A self-intersecting mesh is not
//! certified (ADR 0033). The registry does not check ring orientation: a
//! floating-point area sign is unreliable far from the origin (#274); the
//! conformance suite checks it on fixtures.

use axiolid_core::{Frame3, Point2};

use crate::contract::{SectionLimits, SectionSource};

/// Refusal detail: in-plane faces of one side overlap each other.
pub const REFUSED_OVERLAPPING_IN_PLANE_FACES: &str = "overlapping in-plane faces of one side";
/// Refusal detail: in-plane faces of opposite sides overlap each other.
pub const REFUSED_IN_PLANE_FACES_ON_BOTH_SIDES: &str = "in-plane faces of opposite sides overlap";
/// Refusal detail: more than two cut segments meet at one mesh vertex.
pub const REFUSED_SADDLE_VERTEX: &str = "section through a saddle vertex";

/// Hard bounds on the in-plane region work of one request.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
#[non_exhaustive]
pub struct RegionLimits {
    /// Maximum output regions.
    pub max_regions: usize,
    /// Maximum vertices across all output region rings.
    pub max_region_vertices: usize,
    /// Maximum candidate pairs of in-plane triangles compared for overlap.
    pub max_in_plane_pairs: usize,
}

impl RegionLimits {
    /// Construct explicit region limits.
    pub const fn new(
        max_regions: usize,
        max_region_vertices: usize,
        max_in_plane_pairs: usize,
    ) -> Self {
        Self {
            max_regions,
            max_region_vertices,
            max_in_plane_pairs,
        }
    }
}

/// One section request that reports in-plane faces as regions.
///
/// `limits.max_contours` bounds the chains and `limits.max_output_vertices`
/// bounds the chain points.
#[derive(Debug, Clone, Copy, PartialEq)]
#[non_exhaustive]
pub struct SectionRequest {
    /// Right-handed orthonormal mapping from local `(x,y,0)` to model space.
    pub frame: Frame3,
    /// Source, chain and chain-point bounds.
    pub limits: SectionLimits,
    /// Region bounds.
    pub region_limits: RegionLimits,
}

impl SectionRequest {
    /// Construct a request.
    pub const fn new(frame: Frame3, limits: SectionLimits, region_limits: RegionLimits) -> Self {
        Self {
            frame,
            limits,
            region_limits,
        }
    }
}

/// What an open chain ends on.
#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Hash)]
#[non_exhaustive]
pub enum ChainEnd {
    /// The end lies on the boundary of an in-plane region.
    InPlaneRegion,
}

/// One plane-local cut polyline: closed, or open and ending on a region.
///
/// A closed chain's terminal point is implicit and must not duplicate the
/// first point.
#[derive(Debug, Clone, PartialEq)]
#[non_exhaustive]
pub struct SectionChain {
    /// Plane-local points in traversal order.
    pub points: Vec<Point2>,
    /// What the first and last point end on; `None` for a closed chain.
    pub ends: Option<[ChainEnd; 2]>,
}

impl SectionChain {
    /// Construct a closed chain.
    pub fn closed(points: Vec<Point2>) -> Self {
        Self { points, ends: None }
    }

    /// Construct an open chain from what its first and last point end on.
    pub fn open(points: Vec<Point2>, start: ChainEnd, end: ChainEnd) -> Self {
        Self {
            points,
            ends: Some([start, end]),
        }
    }

    /// Whether the last point connects back to the first.
    pub const fn is_closed(&self) -> bool {
        self.ends.is_none()
    }
}

/// The side of the cutting plane a face's material lies on.
///
/// Declaration order is the sort order.
#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Hash)]
#[non_exhaustive]
pub enum InPlaneSide {
    /// The face normal agrees with the frame's `+z`, so the material (for a
    /// sheet, the back of the face) lies on the frame's `-z` side.
    Below,
    /// The face normal opposes the frame's `+z`; the material lies on `+z`.
    Above,
    /// Faces of both facings coincide. The first provider refuses this.
    Both,
}

/// A planar region of the section lying in the cutting plane.
#[derive(Debug, Clone, PartialEq)]
#[non_exhaustive]
pub struct InPlaneRegion {
    /// The side of the plane the material lies on.
    pub side: InPlaneSide,
    /// Counter-clockwise outer ring.
    pub outer: Vec<Point2>,
    /// Clockwise hole rings.
    pub holes: Vec<Vec<Point2>>,
}

impl InPlaneRegion {
    /// Construct a region.
    pub fn new(side: InPlaneSide, outer: Vec<Point2>, holes: Vec<Vec<Point2>>) -> Self {
        Self { side, outer, holes }
    }
}

/// Auditable counts for one request-based section.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
#[non_exhaustive]
pub struct PlaneSectionEvidence {
    /// Approximation source.
    pub source: SectionSource,
    /// Source triangles inspected.
    pub source_triangles: usize,
    /// Source triangles lying in the plane.
    pub in_plane_triangles: usize,
    /// Plane edges that only touch the plane: nothing is emitted for them.
    pub touching_edges: usize,
    /// Candidate pairs of in-plane triangles compared for overlap.
    pub in_plane_pairs: usize,
    /// Output closed chains.
    pub closed_chains: usize,
    /// Output open chains.
    pub open_chains: usize,
    /// Points across all output chains.
    pub chain_vertices: usize,
    /// Output regions.
    pub regions: usize,
    /// Points across all output region rings.
    pub region_vertices: usize,
}

impl PlaneSectionEvidence {
    /// Record evidence for a section derived from the input mesh; every
    /// other count starts at zero.
    pub const fn input_mesh(source_triangles: usize) -> Self {
        Self {
            source: SectionSource::InputMesh,
            source_triangles,
            in_plane_triangles: 0,
            touching_edges: 0,
            in_plane_pairs: 0,
            closed_chains: 0,
            open_chains: 0,
            chain_vertices: 0,
            regions: 0,
            region_vertices: 0,
        }
    }

    /// Record the counts taken from the source mesh.
    pub const fn with_source_counts(
        mut self,
        in_plane_triangles: usize,
        touching_edges: usize,
        in_plane_pairs: usize,
    ) -> Self {
        self.in_plane_triangles = in_plane_triangles;
        self.touching_edges = touching_edges;
        self.in_plane_pairs = in_plane_pairs;
        self
    }

    /// Record the counts of the output.
    pub fn with_output(mut self, chains: &[SectionChain], regions: &[InPlaneRegion]) -> Self {
        self.closed_chains = chains.iter().filter(|chain| chain.is_closed()).count();
        self.open_chains = chains.len() - self.closed_chains;
        self.chain_vertices = chains.iter().map(|chain| chain.points.len()).sum();
        self.regions = regions.len();
        self.region_vertices = regions
            .iter()
            .map(|region| region.outer.len() + region.holes.iter().map(Vec::len).sum::<usize>())
            .sum();
        self
    }

    /// Whether this result came from the discrete input mesh.
    pub const fn is_derived_from_input_mesh(self) -> bool {
        matches!(self.source, SectionSource::InputMesh)
    }
}

/// Cut chains and in-plane regions plus evidence.
#[derive(Debug, Clone, PartialEq)]
#[non_exhaustive]
pub struct PlaneSection {
    /// Right-handed orthonormal mapping from local `(x,y,0)` to model space.
    pub frame: Frame3,
    /// Transversal cut chains. Empty means no crossing.
    pub chains: Vec<SectionChain>,
    /// Regions of faces lying in the plane.
    pub regions: Vec<InPlaneRegion>,
    /// Source and output counts.
    pub evidence: PlaneSectionEvidence,
}

impl PlaneSection {
    /// Construct a provider result. Registries validate it before returning it.
    pub fn new(
        frame: Frame3,
        chains: Vec<SectionChain>,
        regions: Vec<InPlaneRegion>,
        evidence: PlaneSectionEvidence,
    ) -> Self {
        Self {
            frame,
            chains,
            regions,
            evidence,
        }
    }
}
