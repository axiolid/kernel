#![forbid(unsafe_code)]
//! Portable mesh plane-section operation contract and evidence.

pub mod conformance;
mod contract;
mod fixtures;
mod region_conformance;
mod request;

pub use contract::{
    MeshPlaneSection, SectionContour, SectionEvidence, SectionLimits, SectionOutcome, SectionSource,
};
pub use request::{
    ChainEnd, InPlaneRegion, InPlaneSide, PlaneSection, PlaneSectionEvidence, RegionLimits,
    SectionChain, SectionRequest, REFUSED_IN_PLANE_FACES_ON_BOTH_SIDES,
    REFUSED_OVERLAPPING_IN_PLANE_FACES, REFUSED_SADDLE_VERTEX,
};

pub const CAPABILITY_ID: axiolid_contracts::CapabilityId =
    axiolid_contracts::capability_ids::MESH_SECTION;
