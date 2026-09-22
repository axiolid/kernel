// SPDX-License-Identifier: MPL-2.0

//! Building an arrangement from polygon boundaries.

use axiolid_core::Point2;

use crate::entity::{Face, HalfEdge, Vertex};
use crate::id::{FaceId, HalfEdgeId, VertexId};
use crate::Arrangement;

/// Why an arrangement could not be built.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
#[non_exhaustive]
pub enum BuildError {
    /// Fewer than three distinct corners were supplied.
    DegenerateBoundary,
    /// The boundary encloses no area.
    ZeroArea,
}

impl core::fmt::Display for BuildError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        match self {
            Self::DegenerateBoundary => write!(f, "a face needs at least three corners"),
            Self::ZeroArea => write!(f, "the boundary encloses no area"),
        }
    }
}

impl core::error::Error for BuildError {}

impl Arrangement {
    /// Build a single bounded face from a closed polygon.
    ///
    /// Corners are taken counter-clockwise; a clockwise ring is reversed so
    /// the bounded face is the one inside it. Silently accepting either
    /// winding is the right call here because "which side is inside" is
    /// unambiguous for a simple closed ring, and forcing callers to normalise
    /// first only moves the same code outward.
    ///
    /// # Errors
    ///
    /// [`BuildError::DegenerateBoundary`] or [`BuildError::ZeroArea`].
    pub fn from_polygon(corners: &[Point2]) -> Result<Self, BuildError> {
        if corners.len() < 3 {
            return Err(BuildError::DegenerateBoundary);
        }
        let mut ring: Vec<Point2> = corners.to_vec();
        if signed_area(&ring) < 0.0 {
            ring.reverse();
        }
        if signed_area(&ring) <= 0.0 {
            return Err(BuildError::ZeroArea);
        }

        let mut arrangement = Self::new();
        let n = ring.len();
        for &corner in &ring {
            arrangement.vertices.push(Vertex {
                position: corner,
                outgoing: None,
            });
        }

        // Half-edge `2i` runs corner i -> i+1 with the bounded face on its
        // left; `2i + 1` is its twin, bordering the unbounded face.
        let inner = FaceId::from_index(1);
        for i in 0..n {
            let next_i = (i + 1) % n;
            let prev_i = (i + n - 1) % n;
            arrangement.halfedges.push(HalfEdge {
                origin: VertexId::from_index(i),
                twin: HalfEdgeId::from_index(2 * i + 1),
                next: HalfEdgeId::from_index(2 * next_i),
                prev: HalfEdgeId::from_index(2 * prev_i),
                face: inner,
            });
            arrangement.halfedges.push(HalfEdge {
                origin: VertexId::from_index(next_i),
                twin: HalfEdgeId::from_index(2 * i),
                // The outer boundary runs the opposite way around.
                next: HalfEdgeId::from_index(2 * prev_i + 1),
                prev: HalfEdgeId::from_index(2 * next_i + 1),
                face: FaceId::OUTER,
            });
            arrangement.vertices[i].outgoing = Some(HalfEdgeId::from_index(2 * i));
        }

        arrangement.faces.push(Face {
            boundary: Some(HalfEdgeId::from_index(0)),
        });
        arrangement.faces[0].boundary = Some(HalfEdgeId::from_index(1));
        Ok(arrangement)
    }
}

/// Shoelace area of a ring, about its first corner.
fn signed_area(ring: &[Point2]) -> f64 {
    let base = ring[0];
    let mut twice = 0.0;
    for window in ring[1..].windows(2) {
        let a = window[0] - base;
        let b = window[1] - base;
        twice += a.x * b.y - a.y * b.x;
    }
    twice / 2.0
}
