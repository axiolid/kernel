#![forbid(unsafe_code)]

//! The mesh boolean provider behind
//! [`axiolid_mesh_boolean_contract::MeshBoolean`].
//!
//! The algorithm itself lives in the private `csg` module, absorbed from
//! `boolmesh` 0.1.9 rather than depended upon (ADR 0047, superseding 0014).
//! Absorption
//! makes the hot paths and the known defects reachable: profiling put over
//! 99% of a boolean's runtime inside upstream's single `compute_boolean`
//! entry point, and a depth-2 Menger sponge panics inside it.
//!
//! The absorbed files keep Saki Komikado's copyright headers and MPL-2.0
//! notice. Axiolid is MPL-2.0 itself, so this adds no licence obligation the
//! project had not already accepted -- but the headers must survive
//! refactoring.

mod box_detect;
mod cellular;
mod convert;
mod csg;
mod grouping;
mod provider;

pub use provider::BoolmeshBoolean;
