#![forbid(unsafe_code)]

//! Mesh queries: clearance, near-coplanar faces, thin sheets, containment, ray casting,
//! and genus (#89).
//!
//! Containment is the exact ray-parity test #77 built for its boolean and
//! kept private. It is promoted here rather than written a second time, so
//! the boolean and the public query cannot disagree.

pub mod cast;
pub mod clearance;
pub mod containment;
pub mod coplanar;
pub mod genus;
pub mod overlap;
pub mod planes;
pub mod sheet;
pub mod sight;
pub mod topology;
pub mod winding;

pub use cast::{ray_cast, Hit};
pub use clearance::{min_gap, surface_gap};
pub use containment::{contains, winding_number};
pub use coplanar::{near_coplanar, NearCoplanar};
pub use genus::{genus, GenusError};
pub use overlap::{
    difference_volume, enclosed_volume, intersection_volume, Operand, OverlapError, VolumeInterval,
};
pub use planes::{detect_planes, DetectedPlane, PlaneError, PlaneTolerance};
pub use sheet::{thin_sheet, ThinSheet};
pub use sight::{line_of_sight, line_of_sight_within, Sight, SightError, MAX_SIGHT_CELLS};
pub use topology::{
    topology, ComponentTopology, EdgeLoop, MeshTopology, SurfaceKind, TopologyError,
};
pub use winding::{certify_unit_winding, WindingError};
