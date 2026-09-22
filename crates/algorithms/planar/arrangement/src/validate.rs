// SPDX-License-Identifier: MPL-2.0

//! Structural audit of the arrangement.
//!
//! Every edit is supposed to preserve these invariants. Having them as a
//! checkable report rather than scattered `debug_assert!`s means a test can
//! prove an edit left the structure sound, and a caller integrating a new
//! edit can find out where it went wrong instead of getting a hang.

use crate::id::{FaceId, HalfEdgeId};
use crate::Arrangement;

/// What an audit found.
///
/// Counts rather than a bare bool: "the structure is broken" is not
/// actionable, "three half-edges have a twin that does not point back" is.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Default)]
pub struct ArrangementHealth {
    /// Half-edges whose twin does not name them back.
    pub broken_twins: usize,
    /// Half-edges where `next`'s `prev` is not the half-edge itself.
    pub broken_links: usize,
    /// Half-edges whose `next` lies on a different face.
    pub face_mismatches: usize,
    /// Boundary walks that did not close within the arena bound.
    pub unclosed_faces: usize,
    /// Bounded faces whose signed area is not strictly positive.
    pub inverted_faces: usize,
}

impl ArrangementHealth {
    /// Whether every invariant holds.
    #[must_use]
    pub fn is_sound(&self) -> bool {
        *self == Self::default()
    }
}

impl Arrangement {
    /// Check every structural invariant.
    #[must_use]
    pub fn audit(&self) -> ArrangementHealth {
        let mut health = ArrangementHealth::default();

        for (index, he) in self.halfedges.iter().enumerate() {
            let id = HalfEdgeId::from_index(index);
            if self.halfedges[he.twin.index()].twin != id {
                health.broken_twins += 1;
            }
            if self.halfedges[he.next.index()].prev != id {
                health.broken_links += 1;
            }
            if self.halfedges[he.next.index()].face != he.face {
                health.face_mismatches += 1;
            }
        }

        for index in 0..self.faces.len() {
            let face = FaceId::from_index(index);
            let Some(start) = self.faces[index].boundary else {
                continue;
            };
            // Walk the boundary and confirm it returns to the start.
            let mut steps = 0usize;
            let mut current = self.halfedges[start.index()].next;
            while current != start {
                steps += 1;
                if steps > self.halfedges.len() {
                    health.unclosed_faces += 1;
                    break;
                }
                current = self.halfedges[current.index()].next;
            }
            if face != FaceId::OUTER && self.face_area(face) <= 0.0 {
                health.inverted_faces += 1;
            }
        }

        health
    }
}
