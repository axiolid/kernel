//! Which face a point is in.
//!
//! Every decision is an exact orientation of the point against binary64
//! ring edges (`orient_doubles`), so the answer is exact for the rings as
//! returned. A point within the snap displacement of a boundary may be
//! answered as either neighbour (ADR 0100). The scan is linear in the
//! number of ring edges, with a box test per ring; no performance claim.
//!
//! The answer is the first that holds:
//!
//! 1. On a ring edge of a face: scanning faces, outer ring then holes, edges
//!    in output order.
//! 2. Inside a face: inside its outer ring and in none of its holes. A point
//!    on a dangling barrier or bridge is inside, since the piece bounds
//!    nothing there.
//! 3. In or on a closed gap, likewise.
//! 4. In or on an obstacle: the closed outer ring, not strictly inside a
//!    hole; the lowest index.
//! 5. Otherwise unbounded.

use axiolid_core::Point2;
use axiolid_guarantees::Sign;

use super::{FaceLocation, FreeFace, FreeFaces, FreeFacesError};
use crate::exact_arc::orient_doubles;

/// The box of a ring contains the point, edges included.
fn in_box(ring: &[Point2], at: Point2) -> bool {
    let (mut low, mut high) = (ring[0], ring[0]);
    for point in ring {
        low = Point2::new(low.x.min(point.x), low.y.min(point.y));
        high = Point2::new(high.x.max(point.x), high.y.max(point.y));
    }
    low.x <= at.x && at.x <= high.x && low.y <= at.y && at.y <= high.y
}

/// The first edge of `ring` the point lies on.
fn on_edge(ring: &[Point2], at: Point2) -> Option<usize> {
    if !in_box(ring, at) {
        return None;
    }
    (0..ring.len()).find(|&index| {
        let (a, b) = (ring[index], ring[(index + 1) % ring.len()]);
        orient_doubles(a, b, at) == Sign::Zero
            && a.x.min(b.x) <= at.x
            && at.x <= a.x.max(b.x)
            && a.y.min(b.y) <= at.y
            && at.y <= a.y.max(b.y)
    })
}

/// Whether the winding number of the closed `ring` about a point that is
/// not on it is not zero.
fn encloses(ring: &[Point2], at: Point2) -> bool {
    if !in_box(ring, at) {
        return false;
    }
    let mut winding = 0i32;
    for index in 0..ring.len() {
        let (a, b) = (ring[index], ring[(index + 1) % ring.len()]);
        let side = orient_doubles(a, b, at);
        if a.y <= at.y {
            if b.y > at.y && side == Sign::Positive {
                winding += 1;
            }
        } else if b.y <= at.y && side == Sign::Negative {
            winding -= 1;
        }
    }
    winding != 0
}

/// Whether `at`, on no ring edge of `face`, is inside it.
fn inside(face: &FreeFace, at: Point2) -> bool {
    encloses(&face.outer.points, at) && !face.holes.iter().any(|hole| encloses(&hole.points, at))
}

impl FreeFaces {
    /// Where `point` lies: in a face, on a ring edge, in a closed gap, in an
    /// obstacle, or in none of them.
    ///
    /// # Errors
    ///
    /// [`FreeFacesError::NonFinitePoint`] for a point with a non-finite
    /// coordinate.
    pub fn locate(&self, point: Point2) -> Result<FaceLocation, FreeFacesError> {
        if !point.is_finite() {
            return Err(FreeFacesError::NonFinitePoint);
        }
        for (face, found) in self.faces.iter().enumerate() {
            let rings = std::iter::once(&found.outer).chain(&found.holes);
            for (ring, found) in rings.enumerate() {
                if let Some(edge) = on_edge(&found.points, point) {
                    return Ok(FaceLocation::OnBoundary { face, ring, edge });
                }
            }
        }
        if let Some(face) = self.faces.iter().position(|face| inside(face, point)) {
            return Ok(FaceLocation::Inside(face));
        }
        let in_gap = |gap: &FreeFace| {
            std::iter::once(&gap.outer)
                .chain(&gap.holes)
                .any(|ring| on_edge(&ring.points, point).is_some())
                || inside(gap, point)
        };
        if let Some(gap) = self.closed_gaps.iter().position(in_gap) {
            return Ok(FaceLocation::InClosedGap(gap));
        }
        let in_obstacle = |polygon: &crate::Polygon| {
            let rings = std::iter::once(&polygon.outer).chain(&polygon.holes);
            rings
                .clone()
                .any(|ring| on_edge(&ring.points, point).is_some())
                || (encloses(&polygon.outer.points, point)
                    && !polygon
                        .holes
                        .iter()
                        .any(|hole| encloses(&hole.points, point)))
        };
        if let Some(obstacle) = self.obstacles.iter().position(in_obstacle) {
            return Ok(FaceLocation::InObstacle(obstacle));
        }
        Ok(FaceLocation::Unbounded)
    }
}
