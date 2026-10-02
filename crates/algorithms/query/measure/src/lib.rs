#![forbid(unsafe_code)]

//! Metric-property contracts.
//!
//! Algorithms are generic over representation and return structured failures;
//! an open shell must not silently report a plausible volume.

#[cfg(feature = "exact")]
pub mod exact;
#[cfg(feature = "exact")]
pub mod exact_distance;
#[cfg(feature = "exact")]
mod exact_domain;
#[cfg(feature = "exact")]
mod exact_face;
#[cfg(feature = "exact")]
pub mod exact_hausdorff;
pub mod frechet;
pub mod measure;
pub mod mesh;
pub mod mesh_hausdorff;
pub mod mesh_measure;
pub mod mesh_proximity;
pub mod polyline_hausdorff;
pub mod properties;
pub mod proximity;
pub mod winding;

#[cfg(feature = "exact")]
pub use exact::{exact_properties, ExactMeasureError};
#[cfg(feature = "exact")]
pub use exact_distance::{
    boundary_clearance, boundary_distance, plan_boundary_clearance, plan_boundary_distance,
    plan_overlap, Clearance, DistanceBounds, PlanOverlap,
};
#[cfg(feature = "exact")]
pub use exact_domain::FaceDomain;
#[cfg(feature = "exact")]
pub use exact_hausdorff::{
    boundary_hausdorff_distance, one_sided_boundary_hausdorff,
    one_sided_boundary_hausdorff_with_budget, BoundaryHausdorff,
};
pub use frechet::{
    discrete_frechet_distance, discrete_frechet_distance_2d, frechet_at_most, frechet_at_most_2d,
    frechet_distance, frechet_distance_2d, FrechetError,
};
pub use measure::Measure;
pub use mesh::{
    second_moments, surface_properties, volume_properties, MeshMeasureError, SurfaceProperties,
    VolumeProperties,
};
pub use mesh_hausdorff::{
    hausdorff_distance, one_sided_hausdorff, HausdorffBounds, HausdorffError, MeshHausdorff,
};
pub use mesh_measure::MeshMeasure;
pub use mesh_proximity::{
    mesh_distance, proximity_components, MeshDistance, MeshProximityError, ProximityComponent,
};
pub use polyline_hausdorff::{
    frechet_decide_certified, frechet_decide_certified_2d, one_sided_polyline_hausdorff_distance,
    one_sided_polyline_hausdorff_distance_2d, polyline_hausdorff_distance,
    polyline_hausdorff_distance_2d, FrechetDecision,
};
pub use properties::MassProperties;
pub use proximity::{
    closest_point_on_triangle, closest_points_on_segments, closest_points_on_triangles,
    ClosestPoints3, ProximityError,
};
pub use winding::{WindingError, WindingMesh, WindingNumber};
