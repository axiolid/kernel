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
// Lint posture for the absorbed algorithm (ADR 0047).
//
// `csg` is a faithful port of boolmesh 0.1.9. These lints fire on
// upstream's style rather than on defects: dense numeric kernels take
// many parameters, and several helpers are currently unused because the
// `compose` module that called them was deliberately not absorbed.
// Scoped to this module so the rest of the crate keeps the workspace's
// normal lint level, and revisited as the code is reworked -- making that
// rework possible is why it was absorbed.
#[allow(clippy::too_many_arguments)]
#[allow(clippy::cast_abs_to_unsigned)]
#[allow(clippy::collapsible_else_if)]
#[allow(clippy::manual_swap)]
#[allow(clippy::if_same_then_else)]
#[allow(clippy::mut_range_bound)]
#[allow(clippy::field_reassign_with_default)]
#[allow(clippy::mut_range_bound)]
#[allow(clippy::needless_range_loop)]
#[allow(dead_code)]
mod csg;
mod grouping;
mod provider;

pub use provider::BoolmeshBoolean;
