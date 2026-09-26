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
pub mod measure;
pub mod mesh;
pub mod mesh_measure;
pub mod mesh_proximity;
pub mod properties;
pub mod proximity;
pub mod winding;

#[cfg(feature = "exact")]
pub use exact::{exact_properties, ExactMeasureError};
#[cfg(feature = "exact")]
pub use exact_distance::{boundary_clearance, boundary_distance, Clearance, DistanceBounds};
pub use measure::Measure;
pub use mesh::{
    second_moments, surface_properties, volume_properties, MeshMeasureError, SurfaceProperties,
    VolumeProperties,
};
pub use mesh_measure::MeshMeasure;
pub use mesh_proximity::{
    mesh_distance, proximity_components, MeshDistance, MeshProximityError, ProximityComponent,
};
pub use properties::MassProperties;
pub use proximity::{
    closest_point_on_triangle, closest_points_on_segments, closest_points_on_triangles,
    ClosestPoints3, ProximityError,
};
pub use winding::{WindingError, WindingMesh, WindingNumber};
