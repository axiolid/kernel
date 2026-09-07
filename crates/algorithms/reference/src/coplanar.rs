//! Overlap between two coplanar triangles.
//!
//! # Why this exists
//!
//! Two coplanar faces meet in an AREA, not a curve, so the intersection-curve
//! module cannot describe them and used to refuse the whole operation. Real
//! building models hit this constantly: flush walls, stacked slabs, a column
//! sharing a face with the floor it stands on.
//!
//! # What it computes
//!
//! The overlap polygon of the two triangles, exactly. Both are convex, so the
//! overlap is a convex polygon of at most six vertices, obtained by clipping
//! one against each edge of the other.
//!
//! Two results matter to the caller and are deliberately distinguished:
//!
//! - **Empty** -- the faces share a plane but no area. They contribute
//!   nothing, and refusing the operation over them would be wrong. Measured:
//!   two walls in one plane five metres apart produce sixteen coplanar face
//!   pairs and zero shared area.
//! - **Non-empty** -- the polygon's boundary is where the surfaces stop
//!   coinciding, and those edges constrain the retriangulation of both faces.
//!
//! # Exactness
//!
//! Every inside/outside decision is an `orient2d` sign. Vertex positions of
//! the clipped polygon are computed only after the crossing they represent
//! has been proven to exist, the same discipline the intersection curve uses.

use axiolid_contracts::Sign;
use axiolid_core::{Point2, Point3};

use crate::orient2d;

/// The dominant axis of `normal`, dropped when projecting to 2D.
///
/// Dropping the largest component keeps the projection non-degenerate: the
/// triangle cannot collapse to a line in the remaining two coordinates.
#[must_use]
pub fn dominant_axis(normal: Point3) -> usize {
    let (x, y, z) = (normal.x.abs(), normal.y.abs(), normal.z.abs());
    if x >= y && x >= z {
        0
    } else if y >= z {
        1
    } else {
        2
    }
}

/// Drop `axis` from `point`, giving the 2D projection.
#[must_use]
pub fn project(point: Point3, axis: usize) -> Point2 {
    match axis {
        0 => Point2::new(point.y, point.z),
        1 => Point2::new(point.x, point.z),
        _ => Point2::new(point.x, point.y),
    }
}

/// Lift a 2D point back onto the plane of `reference`, along `axis`.
///
/// The dropped coordinate is recovered from the plane equation rather than
/// carried along, so the lifted point lies on the plane by construction
/// instead of by accumulated arithmetic.
fn lift(flat: Point2, axis: usize, reference: [Point3; 3]) -> Point3 {
    let normal = (reference[1] - reference[0]).cross(reference[2] - reference[0]);
    let origin = reference[0];
    // Solve normal . (p - origin) = 0 for the dropped coordinate.
    match axis {
        0 => {
            let (y, z) = (flat.x, flat.y);
            let x = origin.x - (normal.y * (y - origin.y) + normal.z * (z - origin.z)) / normal.x;
            Point3::new(x, y, z)
        }
        1 => {
            let (x, z) = (flat.x, flat.y);
            let y = origin.y - (normal.x * (x - origin.x) + normal.z * (z - origin.z)) / normal.y;
            Point3::new(x, y, z)
        }
        _ => {
            let (x, y) = (flat.x, flat.y);
            let z = origin.z - (normal.x * (x - origin.x) + normal.y * (y - origin.y)) / normal.z;
            Point3::new(x, y, z)
        }
    }
}

/// The overlap of two coplanar triangles, as a 3D polygon.
///
/// `subject` and `clip` must be coplanar; the caller establishes that with
/// [`triangle_triangle_relation`](crate::triangle_triangle_relation).
///
/// Returns the overlap polygon's vertices in order, or an empty vector when
/// the triangles share a plane but no area. An empty result is a normal
/// answer, not a failure: coplanar faces that do not overlap are common and
/// contribute nothing to a boolean.
///
/// A result with fewer than three vertices is returned empty. Such a polygon
/// has no area -- the triangles meet along an edge or at a point -- and
/// reporting it as an overlap would create degenerate faces downstream.
#[must_use]
pub fn coplanar_overlap(subject: [Point3; 3], clip: [Point3; 3]) -> Vec<Point3> {
    let normal = (clip[1] - clip[0]).cross(clip[2] - clip[0]);
    let axis = dominant_axis(normal);

    let flat_clip = clip.map(|p| project(p, axis));
    // Clip against consistently-wound edges, so 'inside' is one fixed sign
    // rather than depending on how the caller happened to wind the triangle.
    let clip_ring = match orientation(flat_clip) {
        Sign::Negative => [flat_clip[0], flat_clip[2], flat_clip[1]],
        _ => flat_clip,
    };

    let mut polygon: Vec<Point2> = subject.map(|p| project(p, axis)).to_vec();

    for corner in 0..3 {
        let edge_start = clip_ring[corner];
        let edge_end = clip_ring[(corner + 1) % 3];
        polygon = clip_to_halfplane(&polygon, edge_start, edge_end);
        if polygon.is_empty() {
            return Vec::new();
        }
    }

    // Collapse points repeated by clipping through a shared vertex, then
    // reject anything that is not a genuine area.
    polygon.dedup();
    if polygon.len() > 1 && polygon[0] == polygon[polygon.len() - 1] {
        polygon.pop();
    }
    if polygon.len() < 3 {
        return Vec::new();
    }

    polygon.into_iter().map(|p| lift(p, axis, clip)).collect()
}

/// Keep the part of `polygon` on the inside of the directed line a->b.
///
/// Inside means a non-negative `orient2d` sign, so points exactly ON the line
/// are kept. That choice matters: a subject edge lying along a clip edge is a
/// real part of the overlap boundary, and dropping it would open a gap.
fn clip_to_halfplane(polygon: &[Point2], a: Point2, b: Point2) -> Vec<Point2> {
    let mut out = Vec::with_capacity(polygon.len() + 1);

    for index in 0..polygon.len() {
        let current = polygon[index];
        let next = polygon[(index + 1) % polygon.len()];
        let current_side = side_of(a, b, current);
        let next_side = side_of(a, b, next);

        if current_side != Sign::Negative {
            out.push(current);
        }
        // A strict sign change crosses the line, so the crossing point joins
        // the output. Equality or a zero endpoint does not cross: the shared
        // point was already emitted above.
        if (current_side == Sign::Positive && next_side == Sign::Negative)
            || (current_side == Sign::Negative && next_side == Sign::Positive)
        {
            out.push(line_crossing(current, next, a, b));
        }
    }
    out
}

/// Exact side of the directed line a->b that `point` lies on.
fn side_of(a: Point2, b: Point2, point: Point2) -> Sign {
    orient2d(a, b, point)
        .sign()
        .expect("certified predicates are total")
}

/// Winding of a triangle, as an exact sign.
fn orientation([a, b, c]: [Point2; 3]) -> Sign {
    orient2d(a, b, c)
        .sign()
        .expect("certified predicates are total")
}

/// Where segment `start`->`end` crosses the line a->b.
///
/// The caller has already proven the crossing exists with exact signs; this
/// computes only its position, so rounding can move the point slightly but
/// cannot invent or remove it.
fn line_crossing(start: Point2, end: Point2, a: Point2, b: Point2) -> Point2 {
    let edge = b - a;
    let start_height = edge.x * (start.y - a.y) - edge.y * (start.x - a.x);
    let end_height = edge.x * (end.y - a.y) - edge.y * (end.x - a.x);
    let span = start_height - end_height;
    if span == 0.0 {
        // Proven to cross, so this is unreachable; returning an endpoint
        // keeps the polygon finite rather than emitting a NaN.
        return start;
    }
    let t = start_height / span;
    Point2::new(
        start.x + (end.x - start.x) * t,
        start.y + (end.y - start.y) * t,
    )
}
