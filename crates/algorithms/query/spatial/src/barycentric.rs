//! Barycentric and mean-value coordinates, for interpolating values given at
//! the corners of a triangle, a tetrahedron, a polygon or a closed triangle
//! mesh (#143).
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
//!
//! Four polygon schemes are offered, trading generality for positivity:
//! [`mean_value_coordinates2`] answers for any simple polygon, convex or
//! not; [`wachspress_coordinates2`] and [`discrete_harmonic_coordinates2`]
//! require a strictly convex polygon (refused by name otherwise), and the
//! latter can still produce an infinite or cancelled weight at a pole, which
//! is refused the same way mean-value cancellation is.
//! [`mean_value_coordinates3`] extends the 2D mean-value scheme to a closed
//! triangle mesh in space.

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
    /// Wachspress and discrete harmonic coordinates require a strictly
    /// convex polygon; vertex `index` turns the other way.
    NotConvex {
        /// The reflex (or ambiguous) vertex.
        index: usize,
    },
    /// A face of a triangle mesh names vertex `index`, but the mesh has only
    /// `len` vertices.
    VertexIndex {
        /// The face's index.
        face: usize,
        /// The out-of-range vertex index it named.
        index: usize,
        /// The number of vertices given.
        len: usize,
    },
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
            Self::NotConvex { index } => {
                write!(f, "polygon is not convex at vertex {index}")
            }
            Self::VertexIndex { face, index, len } => {
                write!(f, "face {face} names vertex {index}, but only {len} exist")
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

/// `point` on a vertex or an edge of `v`, answered by that boundary's own
/// linear interpolation: shared by [`wachspress_coordinates2`] and
/// [`discrete_harmonic_coordinates2`], both of which reduce to ordinary
/// barycentric interpolation on the boundary, same as
/// [`mean_value_coordinates2`].
fn vertex_or_edge_weights(
    v: &[Point2],
    point: Point2,
    tolerance: Tolerance,
) -> Option<Vec<Scalar>> {
    let n = v.len();
    let linear = tolerance.linear();
    let mut weights = vec![0.0; n];
    let r: Vec<Scalar> = v.iter().map(|&q| (q - point).length()).collect();
    if let Some(i) = (0..n)
        .filter(|&i| r[i] <= linear)
        .min_by(|&i, &j| r[i].total_cmp(&r[j]))
    {
        weights[i] = 1.0;
        return Some(weights);
    }
    for i in 0..n {
        let j = (i + 1) % n;
        let edge = v[j] - v[i];
        let length = edge.length();
        let t = (point - v[i]).dot(edge) / (length * length);
        let s_i = v[i] - point;
        let s_j = v[j] - point;
        if (0.0..=1.0).contains(&t) && s_i.perp_dot(s_j).abs() / length <= linear {
            weights[i] = 1.0 - t;
            weights[j] = t;
            return Some(weights);
        }
    }
    None
}

/// Refuse a polygon with a reflex vertex: Wachspress and discrete harmonic
/// coordinates are defined here only for a strictly convex polygon (as
/// CGAL's `Wachspress_coordinates_2` requires), any winding.
///
/// Each turn `area2(prev, cur, next)` is compared to the polygon's own
/// orientation, with a threshold in the same area units as the linear
/// tolerance times the perimeter, so a numerically flat corner (near-180
/// degrees) is accepted rather than flagged on rounding noise.
fn check_convex(v: &[Point2], tolerance: Tolerance) -> Result<(), BarycentricError> {
    let n = v.len();
    let signed_area: Scalar = (0..n)
        .map(|i| v[i].perp_dot(v[(i + 1) % n]))
        .sum::<Scalar>()
        * 0.5;
    let orientation = if signed_area >= 0.0 { 1.0 } else { -1.0 };
    let perimeter: Scalar = (0..n).map(|i| (v[(i + 1) % n] - v[i]).length()).sum();
    let threshold = tolerance.linear() * perimeter;
    for i in 0..n {
        let prev = v[(i + n - 1) % n];
        let next = v[(i + 1) % n];
        let turn = area2(prev, v[i], next);
        if turn * orientation < -threshold {
            return Err(BarycentricError::NotConvex { index: i });
        }
    }
    Ok(())
}

/// Wachspress coordinates of `point` in a strictly convex polygon, one
/// weight per vertex (Wachspress 1975, in the determinant form of Meyer,
/// Lee, Barr and Desbrun 2002).
///
/// For vertex `i`, `w_i = C_i / (A_{i-1} * A_i)` where `C_i` is twice the
/// area of the corner triangle `(v_{i-1}, v_i, v_{i+1})` and `A_i` is twice
/// the area of `(v_i, v_{i+1}, point)`; the result is `w_i` normalized to sum
/// to one. Inside the polygon the weights are positive, smooth, and
/// reproduce the point; on the boundary, within the linear tolerance, they
/// are the boundary's own linear interpolation. Outside the polygon they are
/// the formula's rational extension wherever it does not divide by zero. The
/// polygon may wind either way.
///
/// # Errors
///
/// [`BarycentricError::NonFinite`], [`BarycentricError::TooFewVertices`],
/// [`BarycentricError::ShortEdge`], [`BarycentricError::SelfIntersecting`],
/// [`BarycentricError::Degenerate`] for a polygon with no area to speak of,
/// [`BarycentricError::NotConvex`] for a polygon that is not strictly
/// convex, and [`BarycentricError::Undefined`] where a denominator vanishes
/// (on the line through a non-adjacent edge, extended past its endpoints).
pub fn wachspress_coordinates2(
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
    check_convex(v, tolerance)?;
    if let Some(w) = vertex_or_edge_weights(v, point, tolerance) {
        return Ok(w);
    }

    let a: Vec<Scalar> = (0..n).map(|i| area2(v[i], v[(i + 1) % n], point)).collect();
    let mut weights = vec![0.0; n];
    let mut sum = 0.0;
    let mut magnitude = 0.0;
    for i in 0..n {
        let prev = (i + n - 1) % n;
        let c = area2(v[prev], v[i], v[(i + 1) % n]);
        let w = c / (a[prev] * a[i]);
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

/// Discrete harmonic coordinates of `point` in a strictly convex polygon,
/// one weight per vertex (Pinkall and Polthier's cotangent weights, applied
/// to the triangle fan `(point, v_i, v_{i+1})`; Floater's survey calls this
/// scheme "discrete harmonic").
///
/// For vertex `i`, sharing a fan edge with the two triangles on either side
/// of it, `w_i = (cot(gamma) + cot(beta)) / r_i^2`, where `beta` and `gamma`
/// are the unsigned angles at `v_i` in those two triangles and `r_i =
/// |point - v_i|`; the result is `w_i` normalized to sum to one. On the
/// boundary, within the linear tolerance, the weights are the boundary's own
/// linear interpolation. Unlike Wachspress coordinates, discrete harmonic
/// weights are not guaranteed positive inside a convex polygon, and have
/// poles where a fan angle is a multiple of a straight angle (point,
/// vertex and a neighbor collinear): both refused as
/// [`BarycentricError::Undefined`], since the caller cannot act on an
/// infinite or cancelled weight either way. The polygon may wind either
/// way.
///
/// # Errors
///
/// As [`wachspress_coordinates2`].
pub fn discrete_harmonic_coordinates2(
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
    check_convex(v, tolerance)?;
    if let Some(w) = vertex_or_edge_weights(v, point, tolerance) {
        return Ok(w);
    }

    let mut weights = vec![0.0; n];
    let mut sum = 0.0;
    let mut magnitude = 0.0;
    for i in 0..n {
        let u = point - v[i];
        let prev = v[(i + n - 1) % n] - v[i];
        let next = v[(i + 1) % n] - v[i];
        let cot_prev = u.dot(prev) / u.perp_dot(prev).abs();
        let cot_next = u.dot(next) / u.perp_dot(next).abs();
        let r2 = u.length_squared();
        let w = (cot_prev + cot_next) / r2;
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

/// Mean-value coordinates of `point` in a closed triangle mesh, one weight
/// per vertex (Floater, Kos and Reimers 2005).
///
/// `faces` lists each triangle as three indices into `vertices`. The mesh's
/// closedness (every edge shared by exactly two faces, consistently
/// wound) is the caller's responsibility, same as `mean_value_coordinates2`
/// leaves polygon repair to its caller: this function does not audit
/// topology, only computes weights from the faces given. On a vertex, that
/// vertex alone gets weight one. When `point` lies exactly on a face's
/// plane, inside that triangle, the weights are that triangle's own 2D
/// barycentric coordinates (all other vertices zero) -- the formula's own
/// closed form for the coplanar case, not a fallback. Elsewhere on a face's
/// plane but outside its triangle, that face contributes nothing (its
/// spherical angles fold flat without enclosing `point`); a mesh with no
/// other face to carry the point returns [`BarycentricError::Undefined`].
///
/// # Errors
///
/// [`BarycentricError::NonFinite`] for a non-finite vertex or point,
/// [`BarycentricError::TooFewVertices`] for fewer than three vertices or no
/// faces, [`BarycentricError::VertexIndex`] for a face naming a vertex index
/// past the end of `vertices`, and [`BarycentricError::Undefined`] where the
/// weights do not sum to a nonzero, finite total.
pub fn mean_value_coordinates3(
    vertices: &[Point3],
    faces: &[[usize; 3]],
    point: Point3,
    tolerance: Tolerance,
) -> Result<Vec<Scalar>, BarycentricError> {
    let n = vertices.len();
    if !point.is_finite() || !vertices.iter().all(|p| p.is_finite()) {
        return Err(BarycentricError::NonFinite);
    }
    if n < 3 || faces.is_empty() {
        return Err(BarycentricError::TooFewVertices { count: n });
    }
    for (f, face) in faces.iter().enumerate() {
        for &index in face {
            if index >= n {
                return Err(BarycentricError::VertexIndex {
                    face: f,
                    index,
                    len: n,
                });
            }
        }
    }

    let linear = tolerance.linear();
    let d: Vec<Scalar> = vertices.iter().map(|&q| (q - point).length()).collect();
    if let Some(i) = (0..n)
        .filter(|&i| d[i] <= linear)
        .min_by(|&i, &j| d[i].total_cmp(&d[j]))
    {
        let mut weights = vec![0.0; n];
        weights[i] = 1.0;
        return Ok(weights);
    }
    let u: Vec<Point3> = vertices
        .iter()
        .zip(&d)
        .map(|(&q, &r)| (q - point) / r)
        .collect();

    let mut weights = vec![0.0; n];
    for &[fi, fj, fk] in faces {
        let (ui, uj, uk) = (u[fi], u[fj], u[fk]);
        let l = [(uj - uk).length(), (uk - ui).length(), (ui - uj).length()];
        let theta = [
            2.0 * (l[0] / 2.0).clamp(-1.0, 1.0).asin(),
            2.0 * (l[1] / 2.0).clamp(-1.0, 1.0).asin(),
            2.0 * (l[2] / 2.0).clamp(-1.0, 1.0).asin(),
        ];
        let h = (theta[0] + theta[1] + theta[2]) / 2.0;
        if core::f64::consts::PI - h < linear {
            // `point` is coplanar with, and inside, this triangle: its own
            // 2D barycentric coordinates, exact at corners.
            let tri = Triangle3 {
                a: vertices[fi],
                b: vertices[fj],
                c: vertices[fk],
            };
            let w = triangle_barycentric3(&tri, point, tolerance)?;
            let mut weights = vec![0.0; n];
            weights[fi] = w[0];
            weights[fj] = w[1];
            weights[fk] = w[2];
            return Ok(weights);
        }
        let sin_theta = [theta[0].sin(), theta[1].sin(), theta[2].sin()];
        if sin_theta.iter().any(|s| s.abs() <= Scalar::EPSILON) {
            // A fan angle at 0 or pi: this face's spherical triangle is
            // degenerate at `point`, so it cannot enclose it. Skip.
            continue;
        }
        let c = [
            2.0 * h.sin() * (h - theta[0]).sin() / (sin_theta[1] * sin_theta[2]) - 1.0,
            2.0 * h.sin() * (h - theta[1]).sin() / (sin_theta[2] * sin_theta[0]) - 1.0,
            2.0 * h.sin() * (h - theta[2]).sin() / (sin_theta[0] * sin_theta[1]) - 1.0,
        ];
        let sign = ui.dot(uj.cross(uk)).signum();
        let mut degenerate = false;
        let s: Vec<Scalar> = c
            .iter()
            .map(|&ci| {
                let value = sign * (1.0 - ci * ci).max(0.0).sqrt();
                if value.abs() <= Scalar::EPSILON {
                    degenerate = true;
                }
                value
            })
            .collect();
        if degenerate {
            // `point` is on this face's plane but outside its triangle: no
            // contribution from this face.
            continue;
        }
        weights[fi] +=
            (theta[0] - c[1] * theta[2] - c[2] * theta[1]) / (d[fi] * sin_theta[1] * s[2]);
        weights[fj] +=
            (theta[1] - c[2] * theta[0] - c[0] * theta[2]) / (d[fj] * sin_theta[2] * s[0]);
        weights[fk] +=
            (theta[2] - c[0] * theta[1] - c[1] * theta[0]) / (d[fk] * sin_theta[0] * s[1]);
    }
    let sum: Scalar = weights.iter().sum();
    let magnitude: Scalar = weights.iter().map(|w| w.abs()).sum();
    if !sum.is_finite() || sum.abs() <= Scalar::EPSILON * n as Scalar * magnitude {
        return Err(BarycentricError::Undefined);
    }
    for w in &mut weights {
        *w /= sum;
    }
    Ok(weights)
}
