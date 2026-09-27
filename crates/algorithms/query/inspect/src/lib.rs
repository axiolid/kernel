#![forbid(unsafe_code)]

//! Mesh queries: clearance, containment, ray casting, and genus (#89).
//!
//! Containment is the exact ray-parity test #77 built for its boolean and
//! kept private. It is promoted here rather than written a second time, so
//! the boolean and the public query cannot disagree.

pub mod cast;
pub mod clearance;
pub mod containment;
pub mod genus;
pub mod overlap;
pub mod sight;

pub use cast::{ray_cast, Hit};
pub use clearance::min_gap;
pub use containment::{contains, winding_number};
pub use genus::{genus, GenusError};
pub use overlap::{
    difference_volume, enclosed_volume, intersection_volume, Operand, OverlapError, VolumeInterval,
};
pub use sight::{line_of_sight, line_of_sight_within, Sight, SightError, MAX_SIGHT_CELLS};
