#![forbid(unsafe_code)]

//! Feature-gated facade for Axiolid geometry.
//!
//! Nothing is compiled unless it is asked for: `default = []`, so a consumer
//! names the capability they need and pays for that alone (kernel#9). Building
//! with no feature at all raises a `compile_error!` listing the options rather
//! than leaving an empty crate whose every call site fails as "not found".
//! `standard` reproduces the pre-0.4 default (`mesh + cpu + integration`) for
//! consumers who want the old behaviour in one line.
//!
//! # Exact geometry is the primary currency
//!
//! No entry point in this facade converts exact geometry into a mesh
//! unless the caller asked for one and supplied a tolerance. Tessellation
//! is a requested OUTPUT, never a fallback:
//!
//! - `generate::GenerationRequest::ExactBRep` carries no tolerance field,
//!   so a mesh cannot be produced for it. A request it cannot satisfy
//!   exactly is refused rather than approximated.
//! - `generate::TessellationRequest` and `TessellationOptions` have no
//!   `Default`. There is no global chord error to fall back to, because
//!   acceptable error depends on source units and downstream use.
//! - A tessellated result carries the tolerance it was built to, so an
//!   approximation cannot later be mistaken for an exact value.
//!
//! `tests/exact_primary.rs` holds this to account (kernel#36).

// No capability feature is enabled, so this crate would compile to an empty
// facade and every call site would fail with a confusing "not found in this
// scope" error instead of naming the real cause (kernel#9). Fail loudly here.
//
// The guard lists ONE feature per capability family rather than every
// feature: each family entry is implied by its own bundle, so a consumer who
// enabled `brep` (which implies `surfaces` -> `curves` -> `linear`) already
// satisfies it. Listing leaves instead would make this fire spuriously.
#[cfg(not(any(
    feature = "mesh",
    feature = "linear",
    feature = "predicates",
    feature = "profiles",
    feature = "curves",
    feature = "primitives",
    feature = "overlay",
    feature = "field",
    feature = "pointcloud",
    feature = "contracts",
    feature = "cpu",
)))]
compile_error!(
    r#"the `axiolid` facade has no capability features enabled, so it exports nothing.

As of 0.4 `default = []`: this crate is pay-for-what-you-use and compiles only
the geometry you name. Add at least one feature in Cargo.toml.

MIGRATING FROM <=0.3 (default was `mesh + cpu + integration`):
    axiolid = { version = "0.4", features = ["standard"] }

PICK BY TASK -- narrowest feature that does the job:
    triangle meshes, no operations .. "mesh"
    line/plane/ray intersection ..... "linear-intersection"
    volume, area, centroid .......... "measure"
    ray casting against a mesh ...... "ray-mesh"
    boolean union/difference ........ "application"
    planar section / contours ....... "application"
    NURBS curves and surfaces ....... "nurbs"
    exact B-rep topology ............ "brep"
    tessellate exact -> mesh ........ "tessellation"
    point clouds to mesh ............ "pointcloud-provider"
    2D fields, morphology ........... "field-ops"
    mesh repair / healing ........... "heal"

BUNDLES -- use when you want breadth over a minimal build:
    "standard" ... the pre-0.4 default: mesh + cpu + integration
    "discrete" ... mesh stack: booleans, sections, measure, spatial, generate
    "parametric" . exact stack: curves, surfaces, topology, NURBS
    "advanced" ... discrete + parametric + heal
    "full" ....... everything, including parallel/simd/gpu adapters

IMPORTANT -- features are ADDITIVE and imply their prerequisites, so name the
capability you need, not its dependencies. "brep" already pulls surfaces,
curves and linear; adding them by hand only widens your build.

Operations that dispatch to a provider (booleans, sections) need BOTH a
contract and a registered provider. "application" is the supported
combination; "mesh-boolean" alone gives you the trait with nothing behind it.

Full matrix: `cargo metadata` on this crate, or docs/architecture/closure-profiles.md
for the exact crate set each profile compiles."#
);

/// Versioned description of the compiled downstream surface.
#[cfg(feature = "integration")]
pub mod integration;

/// Caller-held broad-phase index for repeated ray casts.
// Gated on BOTH: the cache needs the ray crate for the narrow phase and
// the spatial crate for the broad phase. Enabling ray-mesh alone must
// still build, so the cache is absent there and the facade falls back.
#[cfg(all(feature = "ray-mesh", feature = "spatial"))]
pub mod ray_index;

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
