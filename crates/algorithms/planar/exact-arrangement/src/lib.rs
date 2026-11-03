// SPDX-License-Identifier: MPL-2.0
#![forbid(unsafe_code)]
#![warn(missing_docs)]

//! Exact arrangement of straight segments on rational points, with face
//! tracing (ADR 0093, ADR 0094).
//!
//! [`ExactPoints`] stores points in homogeneous coordinates over exact
//! numbers and decides orientation, order and containment exactly. [`arrange`]
//! cuts segments at every crossing, T-junction and collinear overlap into
//! pieces and joins equal points. [`trace`] closes labelled half-edges into
//! rings and nests holes in outer boundaries. No tolerance enters any
//! decision; a caller rounds the points it returns once each, which makes a
//! vertex shared by two rings bit-identical in both.
//!
//! The crate knows no mesh, plane or face. The callers give segments and say
//! which side of a piece carries which label.

mod arrangement;
mod points;
mod ratio;
mod ring;
mod trace;

pub use arrangement::{
    arrange, boxes_overlap, union_box, Arrangement, Crossing, PairBudget, Piece, PieceId, Segment,
    Vertex, VertexId, CROSSING_SCRATCH_BYTES, SPLIT_SCRATCH_BYTES,
};
pub use points::{ExactPoints, Location, PointId, RingBox};
pub use ratio::round_ratio;
pub use ring::{canonical_rotation, canonical_start, point_order, sequence_order};
pub use trace::{trace, HalfEdge, TracedRegion, TracedRing};
