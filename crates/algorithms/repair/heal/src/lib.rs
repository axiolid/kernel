#![forbid(unsafe_code)]

//! Explicit geometry diagnosis and opt-in repair.
//!
//! Healing never runs inside another algorithm. Callers diagnose first, choose
//! a narrow [`RepairPlan`], and retain the resulting [`RepairReport`] for audit.

mod carry;
pub mod diagnose;
pub mod diagnosis;
pub mod intersect;
pub mod mesh;
pub mod repair;
pub mod traits;

pub use diagnose::diagnose;
pub use diagnosis::{Defect, DefectKind, Diagnosis};
pub use intersect::{self_intersections, self_intersections_brute_force, IntersectingPair};
pub use repair::{RepairAction, RepairPlan, RepairReport};
pub use traits::{Diagnose, Repair};
