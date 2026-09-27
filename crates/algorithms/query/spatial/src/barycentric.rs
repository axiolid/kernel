//! Barycentric and mean-value coordinates, for interpolating values given at
//! the corners of a triangle, a tetrahedron or a polygon (#143).
//!
//! Every function returns one weight per corner. The weights sum to one and
//! reproduce the point: `sum(w_i * corner_i) == point`, up to rounding. A
//! point on a corner gets weight exactly one there and zero elsewhere, so
//! corner values are reproduced exactly.
//!
//! Shapes thinner than the linear tolerance are refused rather than given
//! weights that are mostly rounding noise; a triangle's or tetrahedron's
//! thinness is its least altitude. Weights outside `[0, 1]` are not an
//! error: they place the point outside the triangle or tetrahedron, which is
//! how callers extrapolate or test containment.

use axiolid_core::{Point2, Point3, Polygon2, Scalar, Tolerance, Triangle2, Triangle3};

/// Why coordinates could not be computed.
#[derive(Debug, Clone, Copy, PartialEq)]
#[non_exhaustive]
pub enum BarycentricError {
    /// A corner or the query point is not finite.
    NonFinite,
    /// The triangle, tetrahedron or polygon is thinner than the linear
    /// tolerance: its least altitude (for a polygon, twice its area over its
    /// perimeter) is `thickness`.
    Degenerate {
        /// The shape's thickness, in length units.
        thickness: Scalar,
    },
    /// A polygon has fewer than three vertices.
    TooFewVertices {
        /// The number of vertices given.
        count: usize,
    },
    /// Polygon edge `index`, from vertex `index` to the next, is no longer
    /// than the linear tolerance.
    ShortEdge {
        /// The edge's index.
        index: usize,
    },
    /// Polygon edges `first` and `second` cross, touch or fold back onto
    /// each other within the linear tolerance, so the polygon is not simple.
    SelfIntersecting {
        /// The lower edge index.
        first: usize,
        /// The higher edge index.
        second: usize,
    },
    /// The mean-value weights cancel at this point, which can happen only
    /// outside a non-convex polygon: the coordinates are undefined there.
    Undefined,
}

impl core::fmt::Display for BarycentricError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        match self {
            Self::NonFinite => f.write_str("a corner or the query point is not finite"),
            Self::Degenerate { thickness } => write!(
                f,
                "shape is {thickness} thick, not thicker than the linear tolerance"
            ),
            Self::TooFewVertices { count } => {
                write!(f, "polygon has {count} vertices, need at least 3")
            }
            Self::ShortEdge { index } => {
                write!(f, "polygon edge {index} is not longer than the tolerance")
            }
            Self::SelfIntersecting { first, second } => {
                write!(f, "polygon edges {first} and {second} meet")
            }
            Self::Undefined => {
                f.write_str("mean-value weights cancel here (outside a non-convex polygon)")
            }
        }
    }
}

impl core::error::Error for BarycentricError {}

/// Twice the signed area of `x, y, z`, positive counter-clockwise.
fn area2(x: Point2, y: Point2, z: Point2) -> Scalar {
    (y - x).perp_dot(z - x)
}

/// Six times the signed volume of `x, y, z, w`.
fn volume6(x: Point3, y: Point3, z: Point3, w: Point3) -> Scalar {
    (y - x).dot((z - x).cross(w - x))
}

/// Barycentric coordinates of `point` in a 2D triangle, as weights of
/// `a`, `b` and `c`.
///
/// Each weight is the signed area of the sub-triangle opposite its corner
/// over the triangle's area, so the result does not depend on the
/// triangle's winding.
///
/// # Errors
///
/// [`BarycentricError::NonFinite`], or [`BarycentricError::Degenerate`] when
/// the triangle's least altitude is not above `tolerance.linear()`.
pub fn triangle_barycentric2(
    triangle: &Triangle2,
    point: Point2,
    tolerance: Tolerance,
) -> Result<[Scalar; 3], BarycentricError> {
    let Triangle2 { a, b, c } = *triangle;
    if ![a, b, c, point].iter().all(|p| p.is_finite()) {
        return Err(BarycentricError::NonFinite);
    }
    let total = area2(a, b, c);
    let longest = (b - a).length().max((c - b).length()).max((a - c).length());
    thick_enough(total.abs(), longest, tolerance)?;
    Ok([
        area2(point, b, c) / total,
        area2(a, point, c) / total,
        area2(a, b, point) / total,
    ])
}

/// Barycentric coordinates of `point` in a 3D triangle, as weights of
/// `a`, `b` and `c`.
///
/// A point off the triangle's plane gets the coordinates of its orthogonal
/// projection onto that plane: the weights reproduce the projection, not the
/// point.
///
/// # Errors
///
/// As [`triangle_barycentric2`].
pub fn triangle_barycentric3(
    triangle: &Triangle3,
    point: Point3,
    tolerance: Tolerance,
) -> Result<[Scalar; 3], BarycentricError> {
    let Triangle3 { a, b, c } = *triangle;
    if ![a, b, c, point].iter().all(|p| p.is_finite()) {
        return Err(BarycentricError::NonFinite);
    }
    let normal = (b - a).cross(c - a);
    let longest = (b - a).length().max((c - b).length()).max((a - c).length());
    thick_enough(normal.length(), longest, tolerance)?;
    let total = normal.dot(normal);
    Ok([
        normal.dot((b - point).cross(c - point)) / total,
        normal.dot((c - point).cross(a - point)) / total,
        normal.dot((a - point).cross(b - point)) / total,
    ])
}

/// Barycentric coordinates of `point` in a tetrahedron, as weights of its
/// four corners in order.
///
/// Each weight is the signed volume of the sub-tetrahedron opposite its
/// corner over the tetrahedron's volume, so the corner order's handedness
/// does not matter.
///
/// # Errors
///
/// [`BarycentricError::NonFinite`], or [`BarycentricError::Degenerate`] when
/// the tetrahedron's least altitude is not above `tolerance.linear()`.
pub fn tetrahedron_barycentric(
    corners: [Point3; 4],
    point: Point3,
    tolerance: Tolerance,
) -> Result<[Scalar; 4], BarycentricError> {
    let [a, b, c, d] = corners;
    if ![a, b, c, d, point].iter().all(|p| p.is_finite()) {
        return Err(BarycentricError::NonFinite);
    }
    let total = volume6(a, b, c, d);
    // The least altitude is the volume over the largest face.
    let largest_face = [(b, c, d), (a, c, d), (a, b, d), (a, b, c)]
        .iter()
        .map(|&(x, y, z)| (y - x).cross(z - x).length())
        .fold(0.0, Scalar::max);
    thick_enough(total.abs(), largest_face, tolerance)?;
    Ok([
        volume6(point, b, c, d) / total,
        volume6(a, point, c, d) / total,
        volume6(a, b, point, d) / total,
        volume6(a, b, c, point) / total,
    ])
}

/// Refuse a shape whose `measure / base` (an altitude) is not above the
/// linear tolerance.
fn thick_enough(
    measure: Scalar,
    base: Scalar,
    tolerance: Tolerance,
) -> Result<(), BarycentricError> {
    let thickness = if base > 0.0 { measure / base } else { 0.0 };
    if thickness > tolerance.linear() {
        Ok(())
    } else {
        Err(BarycentricError::Degenerate { thickness })
    }
}

/// Mean-value coordinates of `point` in a simple polygon, one weight per
/// vertex (Floater 2003, in the form for arbitrary polygons of Hormann and
/// Floater 2006).
///
/// Inside a simple polygon, convex or not, the weights are positive-sum,
/// smooth, and reproduce the point. On the boundary, within the linear
/// tolerance, they are the boundary's own linear interpolation: one at a
/// vertex, or the two endpoint weights of an edge. Outside the polygon they
/// are defined wherever the raw weights do not cancel, which is everywhere
/// outside a convex polygon. The polygon may wind either way.
///
/// The polygon is checked for simplicity pairwise, `O(n^2)` in its vertex
/// count, since mean-value coordinates of a crossing polygon mean nothing.
///
/// # Errors
///
/// [`BarycentricError::NonFinite`], [`BarycentricError::TooFewVertices`],
/// [`BarycentricError::ShortEdge`], [`BarycentricError::SelfIntersecting`],
/// [`BarycentricError::Degenerate`] for a polygon with no area to speak of,
/// and [`BarycentricError::Undefined`] where the weights cancel.
pub fn mean_value_coordinates2(
    polygon: &Polygon2,
    point: Point2,
    tolerance: Tolerance,
) -> Result<Vec<Scalar>, BarycentricError> {
    let v = &polygon.vertices;
    let n = v.len();
    if !point.is_finite() || !v.iter().all(|p| p.is_finite()) {
        return Err(BarycentricError::NonFinite);
    }
    if n < 3 {
        return Err(BarycentricError::TooFewVertices { count: n });
    }
    check_simple(v, tolerance)?;
    let perimeter: Scalar = (0..n).map(|i| (v[(i + 1) % n] - v[i]).length()).sum();
    thick_enough(2.0 * polygon.signed_area().abs(), perimeter, tolerance)?;

    let linear = tolerance.linear();
    let mut weights = vec![0.0; n];
    // On a vertex: that vertex alone.
    let s: Vec<Point2> = v.iter().map(|&q| q - point).collect();
    let r: Vec<Scalar> = s.iter().map(|q| q.length()).collect();
    if let Some(i) = (0..n)
        .filter(|&i| r[i] <= linear)
        .min_by(|&i, &j| r[i].total_cmp(&r[j]))
    {
        weights[i] = 1.0;
        return Ok(weights);
    }
    // On an edge: its linear interpolation.
    for i in 0..n {
        let j = (i + 1) % n;
        let edge = v[j] - v[i];
        let length = edge.length();
        let t = (point - v[i]).dot(edge) / (length * length);
        if (0.0..=1.0).contains(&t) && s[i].perp_dot(s[j]).abs() / length <= linear {
            weights[i] = 1.0 - t;
            weights[j] = t;
            return Ok(weights);
        }
    }
    // tan(alpha_i / 2) for the signed angle alpha_i the edge from vertex i
    // to vertex i + 1 subtends at the point: A / (r_i r_j + D), which stays
    // finite because the point is on no edge.
    let half_tangent: Vec<Scalar> = (0..n)
        .map(|i| {
            let j = (i + 1) % n;
            s[i].perp_dot(s[j]) / (r[i] * r[j] + s[i].dot(s[j]))
        })
        .collect();
    let mut sum = 0.0;
    let mut magnitude = 0.0;
    for i in 0..n {
        let w = (half_tangent[(i + n - 1) % n] + half_tangent[i]) / r[i];
        weights[i] = w;
        sum += w;
        magnitude += w.abs();
    }
    if !sum.is_finite() || sum.abs() <= Scalar::EPSILON * n as Scalar * magnitude {
        return Err(BarycentricError::Undefined);
    }
    for w in &mut weights {
        *w /= sum;
    }
    Ok(weights)
}

/// Refuse short edges, and non-adjacent edges that come within the linear
/// tolerance of each other.
///
/// Adjacent edges need no test of their own. Two that fold back onto each
/// other leave the far end of the shorter on the longer, and that end also
/// belongs to an edge not adjacent to the longer, which this test sees; in a
/// triangle, where every pair is adjacent, a fold leaves no area and the
/// thickness test refuses it.
fn check_simple(v: &[Point2], tolerance: Tolerance) -> Result<(), BarycentricError> {
    let n = v.len();
    let linear = tolerance.linear();
    let edge = |i: usize| (v[i], v[(i + 1) % n]);
    for i in 0..n {
        let (p, q) = edge(i);
        if (q - p).length() <= linear {
            return Err(BarycentricError::ShortEdge { index: i });
        }
    }
    for i in 0..n {
        // Skip the edge itself, its successor, and (from edge 0) the last
        // edge, which precedes it.
        let last = if i == 0 { n - 1 } else { n };
        for j in i + 2..last {
            let (p, q) = edge(i);
            let (r, s) = edge(j);
            if segment_distance(p, q, r, s) <= linear {
                return Err(BarycentricError::SelfIntersecting {
                    first: i,
                    second: j,
                });
            }
        }
    }
    Ok(())
}

/// Distance between segments `pq` and `rs`: zero when they cross.
fn segment_distance(p: Point2, q: Point2, r: Point2, s: Point2) -> Scalar {
    let o1 = area2(p, q, r);
    let o2 = area2(p, q, s);
    let o3 = area2(r, s, p);
    let o4 = area2(r, s, q);
    if o1 * o2 < 0.0 && o3 * o4 < 0.0 {
        return 0.0;
    }
    point_segment(r, p, q)
        .min(point_segment(s, p, q))
        .min(point_segment(p, r, s))
        .min(point_segment(q, r, s))
}

fn point_segment(x: Point2, p: Point2, q: Point2) -> Scalar {
    let d = q - p;
    let t = ((x - p).dot(d) / d.dot(d)).clamp(0.0, 1.0);
    (p + d * t - x).length()
}
