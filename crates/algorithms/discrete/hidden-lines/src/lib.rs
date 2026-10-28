#![forbid(unsafe_code)]
//! Exact orthographic hidden-line drawings of triangle meshes (ADR 0090).
//!
//! [`ScalarHiddenLines`] implements
//! [`MeshHiddenLines`](axiolid_mesh_hidden_lines_contract::MeshHiddenLines).
//! Visibility is decided exactly in the view's affine frame: the image of a
//! vertex is an exact polynomial of the binary64 inputs and the frame, every
//! visibility question along an edge is an affine function of the edge
//! parameter, and each sign is decided by the interval filter and, when that
//! cannot decide, by dyadic big-integer arithmetic (`axiolid-exact`). Nothing
//! is rounded before a decision; only the endpoints of the output pieces are,
//! each the exact point rounded once to the nearest binary64.
//!
//! Module map: `view` (the affine map), `recipe` (the sign questions),
//! `param` (parameters along an edge and their interval algebra), `edges`
//! (topology and edge kinds), `broad` (candidate triangles from a BVH), `occlusion` (one edge against the triangles), `caps` (the cross-sections
//! of solids a near plane cuts),
//! `output` (correctly rounded endpoints), `provider` (validation and the
//! provider itself).

mod broad;
mod caps;
mod decide;
mod edges;
mod occlusion;
mod output;
mod param;
mod provider;
mod recipe;
mod view;

pub use provider::ScalarHiddenLines;
#[doc(hidden)]
pub use provider::ScanHiddenLines;
