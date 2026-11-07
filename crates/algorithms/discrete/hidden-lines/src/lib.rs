#![forbid(unsafe_code)]
//! Exact orthographic and perspective hidden-line drawings of triangle meshes
//! (ADR 0090, 0107).
//!
//! [`ScalarHiddenLines`] implements
//! [`MeshHiddenLines`](axiolid_mesh_hidden_lines_contract::MeshHiddenLines).
//! Visibility is decided exactly in the view's affine frame: the image of a
//! vertex is an exact polynomial of the binary64 inputs and the frame, every
//! visibility question along an edge is an affine function of the edge
//! parameter, and each sign is decided by the interval filter and, when that
//! cannot decide, by dyadic big-integer arithmetic (`axiolid-exact`). Nothing
//! is rounded before a decision; only the endpoints of the output pieces are,
//! each the exact point rounded once to the nearest binary64. With
//! `EdgeSelection::group_boundary` it also draws the edges between face
//! groups, and `FaceGroups::Planar` derives the groups exactly (ADR 0106).
//!
//! A perspective view (ADR 0107) is a pinhole at the frame origin: the image
//! of `(U, V, W)` is `(U / -W, V / -W)`. Every sign question is then an exact
//! determinant of view-space vectors from the eye (`recipe`), affine in the
//! edge parameter as before, so nothing divides before it is decided and a
//! triangle behind the eye needs no case; the endpoints are the exact
//! rationals rounded once, and a near-plane cap is decided by ray parity along
//! the eye ray (`caps`). The visible-face provider below is orthographic only.
//!
//! [`ScalarVisibleFaces`] implements
//! [`MeshVisibleFaces`](axiolid_mesh_hidden_lines_contract::MeshVisibleFaces)
//! (ADR 0094) on the same view: the regions of the image where a face group
//! (or the cap of a solid cut by the near plane) is frontmost. It cuts the
//! image with an exact arrangement (`axiolid-exact-arrangement`) of the
//! visible edge pieces, clip segments, cross-sections of interpenetrating
//! triangles and the crop box, labels each side of each piece exactly, and
//! snap-rounds the regions to the binary64 grid (ADR 0100): where no edge
//! passes through the rounding cell of another vertex, each region vertex is
//! the exact point rounded once and equals the hidden-line endpoint of the
//! same point.
//!
//! Module map: `view` (the affine map), `recipe` (the sign questions),
//! `param` (parameters along an edge and their interval algebra), `edges`
//! (topology and edge kinds), `groups` (face groups, planar ones derived), `broad` (candidate triangles from a BVH), `occlusion` (one edge against the triangles), `caps` (the cross-sections
//! of solids a near plane cuts),
//! `output` (correctly rounded endpoints), `prepare` (request validation and the
//! scene shared with the visible-face operation), `provider` (the hidden-line
//! provider itself), `faces` (the visible-face provider; its module docs give
//! the algorithm).

mod broad;
mod caps;
mod decide;
mod edges;
mod faces;
mod groups;
mod occlusion;
mod output;
mod param;
mod prepare;
mod provider;
mod recipe;
mod view;

pub use faces::ScalarVisibleFaces;
#[doc(hidden)]
pub use faces::{ScanVisibleFaces, UnprunedVisibleFaces};
pub use provider::ScalarHiddenLines;
#[doc(hidden)]
pub use provider::ScanHiddenLines;
