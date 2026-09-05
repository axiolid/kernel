#![forbid(unsafe_code)]

//! Feature-gated facade for Axiolid geometry.
//!
//! The default build is intentionally small: core values, meshes, and the
//! portable CPU backend shell. Exact curves/surfaces/topology, algorithms,
//! parallel scheduling, and GPU adapters are opt-in. Leaf crates remain public
//! for consumers that want an even narrower dependency graph.

/// Versioned description of the compiled downstream surface.
#[cfg(feature = "integration")]
pub mod integration;

/// Supported provider-neutral application boundary.
#[cfg(feature = "application")]
pub mod application;

/// Always-available scalar, transform, and bounds vocabulary.
pub mod core {
    pub use axiolid_core::*;
}

pub use axiolid_core::{Aabb, Point2, Point3, Scalar, Tolerance, Transform3, Vec2, Vec3};

#[cfg(feature = "mesh")]
pub mod mesh {
    pub use axiolid_mesh::*;
}

/// Point-sampled geometry: scans and photogrammetry captures.
///
/// Representation only. Parsing LAS/LAZ/E57/PCD/COPC is deliberately out of
/// the kernel (ADR 0038); an ingestion crate converts a file into this value
/// and the kernel never sees the source format.
#[cfg(feature = "pointcloud")]
pub mod pointcloud {
    pub use axiolid_pointcloud::*;
}

/// Turning point-sampled geometry into a surface.
///
/// The contract and its types. A concrete provider arrives with the
/// `pointcloud-provider` feature, or a caller supplies its own.
#[cfg(feature = "pointcloud-reconstruction")]
pub mod pointcloud_reconstruction {
    pub use axiolid_pointcloud_reconstruction_contract::*;

    /// Reference provider: a signed-distance field extracted as a level set.
    #[cfg(feature = "pointcloud-provider")]
    pub use axiolid_pointcloud_reconstruction_sdf::SdfReconstruction;
}

#[cfg(feature = "linear")]
pub mod linear {
    pub use axiolid_linear::*;
}

#[cfg(feature = "predicates")]
pub mod predicates {
    pub use axiolid_predicates::*;
}

/// Certified linear intersection queries.
///
/// The facade adds one compilation unit by design. An application that wants
/// the smallest possible closure should depend on `axiolid-linear` and
/// `axiolid-linear-intersection` directly (ADR 0036).
#[cfg(feature = "linear-intersection")]
pub mod linear_intersection {
    pub use axiolid_linear_intersection::*;
}

#[cfg(feature = "profiles")]
pub mod profile {
    pub use axiolid_profile::*;
}

#[cfg(feature = "curves")]
pub mod curve {
    pub use axiolid_curve::*;
}

#[cfg(feature = "surfaces")]
pub mod surface {
    pub use axiolid_surface::*;
}

/// General NURBS analysis, inverse-query, and exact transformation algorithms.
/// Scalar evaluation for parametric curves and surfaces.
#[cfg(feature = "evaluate")]
pub mod evaluate {
    pub use axiolid_evaluate::*;
}

#[cfg(feature = "nurbs")]
pub mod nurbs {
    pub use axiolid_nurbs::*;
}

#[cfg(feature = "topology")]
pub mod topology {
    pub use axiolid_topology::*;
}

/// Strict exact B-rep results with typed analytic support catalogs and native
/// trim intervals. Tessellation remains an explicit, tolerance-bearing output.
#[cfg(feature = "brep")]
pub mod brep {
    pub use axiolid_brep::*;
}

#[cfg(feature = "model")]
pub mod model {
    pub use axiolid_model::*;
}

#[cfg(feature = "primitives")]
pub mod primitive {
    pub use axiolid_primitive::*;
}

#[cfg(feature = "tessellation")]
pub mod tessellation {
    pub use axiolid_tessellation_contract::*;
}

#[cfg(feature = "spatial")]
pub mod spatial {
    pub use axiolid_spatial::*;
}

/// Narrow-phase ray/triangle-mesh nearest-hit intersection.
///
/// Composes with the `spatial` broad phase: a caller can feed BVH candidate
/// keys into `nearest_hit_among`, or scan a whole mesh with `nearest_hit`.
#[cfg(feature = "ray-mesh")]
pub mod ray_mesh {
    pub use axiolid_ray_mesh::*;
}

#[cfg(feature = "measure")]
pub mod measure {
    pub use axiolid_measure::*;
}

#[cfg(feature = "overlay")]
pub mod overlay {
    pub use axiolid_overlay::*;
}

/// Planar projection of meshes: projected outlines and vertical prism
/// intersection.
#[cfg(feature = "project")]
pub mod project {
    pub use axiolid_project::*;
}

/// Exact planar shortest paths over a visibility graph, with typed
/// unreachable reasons.
///
/// The kernel reports that no route exists under a given envelope. It never
/// reports that a design is non-compliant: that reading belongs to the
/// consumer, not to geometry.
#[cfg(feature = "route")]
pub mod route {
    pub use axiolid_route::*;
}

/// Geometry generation: discrete sweeps plus focused certified trimmed arrangements.
///
/// Broad profile/path generators still return explicit meshes. The certified affine
/// surface-pair constructor returns an analytic `ExactBRep` arrangement with an
/// explicit residual certificate and never substitutes a mesh fallback (ADR 0029).
#[cfg(feature = "generate")]
pub mod generate {
    pub use axiolid_construct::*;
}

/// Frame-neutral sampled layered-field values and configuration.
#[cfg(feature = "field")]
pub mod field {
    pub use axiolid_field::*;
}

/// Sampling, morphology, clearance, and optional navigation over layered fields.
#[cfg(feature = "field-ops")]
pub mod field_ops {
    pub use axiolid_field_ops::*;
}

#[cfg(feature = "heal")]
pub mod heal {
    pub use axiolid_heal::*;
}

#[cfg(feature = "contracts")]
pub mod contracts {
    pub use axiolid_contracts::*;
}

#[cfg(feature = "mesh-contracts")]
pub mod mesh_contracts {
    pub use axiolid_mesh_contracts::*;
}

#[cfg(feature = "mesh-boolean")]
pub mod mesh_boolean {
    pub use axiolid_mesh_boolean_contract::*;
}

#[cfg(feature = "mesh-section")]
pub mod mesh_section {
    pub use axiolid_mesh_section_contract::*;
}

#[cfg(feature = "graph-compile")]
pub mod graph_compile {
    pub use axiolid_mesh_compile_contract::*;
}

#[cfg(any(feature = "dispatch-mesh-boolean", feature = "dispatch-mesh-section"))]
pub mod dispatch {
    pub use axiolid_dispatch::*;
}

#[cfg(feature = "cpu")]
pub mod cpu {
    pub use axiolid_backend_cpu::*;
}

#[cfg(feature = "gpu")]
pub mod gpu {
    pub use axiolid_backend_gpu::*;
}
