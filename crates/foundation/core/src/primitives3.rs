//! Bounded three-dimensional primitives.
//!
//! These are plain data with public fields, matching the rest of [`crate`]:
//! adapters may construct degenerate or dirty values, and algorithms validate
//! against the operation's own tolerance rather than a global epsilon hidden
//! here. Nothing in this module refuses input.
//!
//! There is deliberately no read-only/mutable type pair. Languages without
//! ownership express that distinction by shipping `Triangle3` beside a
//! `MTriangle3`; Rust expresses it as `&T` versus `&mut T` on one type, checked
//! by the compiler. Mirroring every type would double the surface and prove
//! nothing the borrow checker does not already prove.

use crate::{Point3, Scalar, Vec3};

/// A three-dimensional triangle.
///
/// Corner order defines the winding, and therefore the sign of [`normal`]. A
/// degenerate triangle -- collinear or coincident corners -- is representable
/// on purpose: whether that is an error depends on the caller's tolerance, so
/// the decision belongs to the algorithm consuming it.
///
/// [`normal`]: Triangle3::normal
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Triangle3 {
    /// First corner.
    pub a: Point3,
    /// Second corner.
    pub b: Point3,
    /// Third corner.
    pub c: Point3,
}

impl Triangle3 {
    /// Construct a triangle without validating it.
    pub const fn new(a: Point3, b: Point3, c: Point3) -> Self {
        Self { a, b, c }
    }

    /// Unnormalized right-handed normal.
    ///
    /// Left unnormalized because its magnitude is twice the triangle area, so
    /// callers that need area get it without a second cross product, and
    /// callers that need a direction normalize explicitly against their own
    /// tolerance. A degenerate triangle yields a zero vector rather than a
    /// NaN-bearing unit vector.
    pub fn normal(&self) -> Vec3 {
        (self.b - self.a).cross(self.c - self.a)
    }

    /// Triangle area, always non-negative.
    pub fn area(&self) -> Scalar {
        self.normal().length() * 0.5
    }

    /// Centroid of the three corners.
    pub fn centroid(&self) -> Point3 {
        (self.a + self.b + self.c) / 3.0
    }
}

/// A planar rectangle in three-dimensional space.
///
/// Stored as a corner and two edge vectors rather than four corners, so the
/// parallelogram property holds by construction and cannot drift as corners
/// are edited independently. Whether `x` and `y` are perpendicular is not
/// enforced: an importer may produce a sheared quad, and refusing it here
/// would lose data the caller may still want to inspect or repair.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Rectangle3 {
    /// Corner the edge vectors start from.
    pub origin: Point3,
    /// First edge vector, from `origin`.
    pub x: Vec3,
    /// Second edge vector, from `origin`.
    pub y: Vec3,
}

impl Rectangle3 {
    /// Construct from a corner and two edge vectors.
    pub const fn new(origin: Point3, x: Vec3, y: Vec3) -> Self {
        Self { origin, x, y }
    }

    /// The four corners, counter-clockwise about [`normal`].
    ///
    /// [`normal`]: Rectangle3::normal
    pub fn corners(&self) -> [Point3; 4] {
        [
            self.origin,
            self.origin + self.x,
            self.origin + self.x + self.y,
            self.origin + self.y,
        ]
    }

    /// Unnormalized normal of the plane the rectangle lies in.
    pub fn normal(&self) -> Vec3 {
        self.x.cross(self.y)
    }

    /// Area of the rectangle, which is the parallelogram area when sheared.
    pub fn area(&self) -> Scalar {
        self.normal().length()
    }

    /// Centre point.
    pub fn center(&self) -> Point3 {
        self.origin + (self.x + self.y) * 0.5
    }
}

/// An oriented box: the generalization of a rectangle to three dimensions.
///
/// Distinct from [`Aabb`], which is axis-aligned and exists for broad-phase
/// rejection. This one carries its own orientation, so it can bound a rotated
/// object tightly where an axis-aligned box would not.
///
/// The three edge vectors are not required to be mutually perpendicular, for
/// the same reason [`Rectangle3`] does not require it: a dirty import stays
/// representable and inspectable.
///
/// [`Aabb`]: crate::Aabb
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Box3 {
    /// Corner the edge vectors start from.
    pub origin: Point3,
    /// First edge vector, from `origin`.
    pub x: Vec3,
    /// Second edge vector, from `origin`.
    pub y: Vec3,
    /// Third edge vector, from `origin`.
    pub z: Vec3,
}

impl Box3 {
    /// Construct from a corner and three edge vectors.
    pub const fn new(origin: Point3, x: Vec3, y: Vec3, z: Vec3) -> Self {
        Self { origin, x, y, z }
    }

    /// The eight corners.
    ///
    /// Ordered so bit 0 selects `x`, bit 1 selects `y`, and bit 2 selects `z`:
    /// index 0 is `origin` and index 7 is the far corner. That ordering lets a
    /// caller index a corner by axis mask instead of memorising a winding.
    pub fn corners(&self) -> [Point3; 8] {
        let mut out = [self.origin; 8];
        for (index, corner) in out.iter_mut().enumerate() {
            let mut point = self.origin;
            if index & 1 != 0 {
                point += self.x;
            }
            if index & 2 != 0 {
                point += self.y;
            }
            if index & 4 != 0 {
                point += self.z;
            }
            *corner = point;
        }
        out
    }

    /// Signed volume. Negative when the edge vectors are left-handed.
    ///
    /// Signed rather than absolute so a caller can detect an inverted frame,
    /// which is a common symptom of a mirrored or badly converted import.
    pub fn signed_volume(&self) -> Scalar {
        self.x.cross(self.y).dot(self.z)
    }

    /// Centre point.
    pub fn center(&self) -> Point3 {
        self.origin + (self.x + self.y + self.z) * 0.5
    }
}

/// A simple polygon in three-dimensional space, without holes.
///
/// The vertices are assumed coplanar and non-self-intersecting, and neither is
/// checked here: both are tolerance-dependent judgements. The closing edge is
/// implicit, so the last vertex joins the first and repeating it creates a
/// zero-length edge rather than a closed ring.
#[derive(Debug, Clone, PartialEq)]
pub struct Polygon3 {
    /// Boundary vertices in order. The closing edge is implicit.
    pub vertices: Vec<Point3>,
}

impl Polygon3 {
    /// Construct from boundary vertices in order.
    pub const fn new(vertices: Vec<Point3>) -> Self {
        Self { vertices }
    }

    /// Number of boundary vertices, which equals the number of edges.
    pub fn len(&self) -> usize {
        self.vertices.len()
    }

    /// Whether the polygon has no vertices.
    pub fn is_empty(&self) -> bool {
        self.vertices.is_empty()
    }

    /// Twice the vector area, summed by the shoelace formula in three
    /// dimensions.
    ///
    /// The direction is the polygon's normal and the magnitude is twice its
    /// area, so this is the 3D analogue of the signed 2D shoelace sum. It is
    /// exact for a planar polygon and degrades gracefully for a slightly
    /// non-planar one, which is what imported data usually is.
    ///
    /// Fewer than three vertices enclose nothing and yield a zero vector.
    pub fn vector_area2(&self) -> Vec3 {
        if self.vertices.len() < 3 {
            return Vec3::ZERO;
        }
        // Summed about the first vertex rather than the origin: the origin can
        // be arbitrarily far from the data in a georeferenced model, and the
        // cancellation that follows costs precision for no benefit.
        let base = self.vertices[0];
        let mut total = Vec3::ZERO;
        for pair in self.vertices[1..].windows(2) {
            total += (pair[0] - base).cross(pair[1] - base);
        }
        total
    }

    /// Unnormalized normal implied by the vertex winding.
    pub fn normal(&self) -> Vec3 {
        self.vector_area2()
    }

    /// Polygon area, always non-negative.
    pub fn area(&self) -> Scalar {
        self.vector_area2().length() * 0.5
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn triangle_normal_follows_winding_and_area_is_half_its_length() {
        let triangle = Triangle3::new(
            Point3::ZERO,
            Point3::new(2.0, 0.0, 0.0),
            Point3::new(0.0, 3.0, 0.0),
        );
        assert_eq!(triangle.normal(), Vec3::new(0.0, 0.0, 6.0));
        assert_eq!(triangle.area(), 3.0);

        // Reversing the winding flips the normal but not the area.
        let reversed = Triangle3::new(triangle.a, triangle.c, triangle.b);
        assert_eq!(reversed.normal(), -triangle.normal());
        assert_eq!(reversed.area(), triangle.area());
    }

    #[test]
    fn degenerate_triangle_has_zero_normal_rather_than_nan() {
        let collinear = Triangle3::new(
            Point3::ZERO,
            Point3::new(1.0, 1.0, 1.0),
            Point3::new(2.0, 2.0, 2.0),
        );
        assert_eq!(collinear.normal(), Vec3::ZERO);
        assert_eq!(collinear.area(), 0.0);
    }

    #[test]
    fn rectangle_corners_close_the_loop_and_area_matches_the_edges() {
        let rectangle = Rectangle3::new(
            Point3::new(1.0, 0.0, 0.0),
            Vec3::new(2.0, 0.0, 0.0),
            Vec3::new(0.0, 4.0, 0.0),
        );
        let corners = rectangle.corners();
        assert_eq!(corners[0], rectangle.origin);
        // Opposite corners share a midpoint when the quad really is planar.
        assert_eq!((corners[0] + corners[2]) * 0.5, rectangle.center());
        assert_eq!((corners[1] + corners[3]) * 0.5, rectangle.center());
        assert_eq!(rectangle.area(), 8.0);
    }

    #[test]
    fn box_corner_index_selects_edges_by_bit() {
        let unit = Box3::new(Point3::ZERO, Vec3::X, Vec3::Y, Vec3::Z);
        let corners = unit.corners();
        assert_eq!(corners[0], Point3::ZERO);
        assert_eq!(corners[1], Vec3::X);
        assert_eq!(corners[2], Vec3::Y);
        assert_eq!(corners[4], Vec3::Z);
        assert_eq!(corners[7], Vec3::ONE);
        assert_eq!(unit.center(), Vec3::splat(0.5));
    }

    #[test]
    fn box_volume_is_signed_so_a_mirrored_frame_is_detectable() {
        let right_handed = Box3::new(Point3::ZERO, Vec3::X, Vec3::Y, Vec3::Z);
        assert_eq!(right_handed.signed_volume(), 1.0);

        let mirrored = Box3::new(Point3::ZERO, Vec3::Y, Vec3::X, Vec3::Z);
        assert_eq!(mirrored.signed_volume(), -1.0);
    }

    #[test]
    fn polygon_area_is_winding_independent_and_normal_is_not() {
        let square = Polygon3::new(vec![
            Point3::ZERO,
            Point3::new(2.0, 0.0, 0.0),
            Point3::new(2.0, 2.0, 0.0),
            Point3::new(0.0, 2.0, 0.0),
        ]);
        assert_eq!(square.area(), 4.0);
        assert_eq!(square.normal(), Vec3::new(0.0, 0.0, 8.0));

        let mut reversed = square.vertices.clone();
        reversed.reverse();
        let reversed = Polygon3::new(reversed);
        assert_eq!(reversed.area(), square.area());
        assert_eq!(reversed.normal(), -square.normal());
    }

    #[test]
    fn polygon_with_fewer_than_three_vertices_encloses_nothing() {
        assert_eq!(Polygon3::new(Vec::new()).vector_area2(), Vec3::ZERO);
        assert_eq!(
            Polygon3::new(vec![Point3::ZERO, Vec3::X]).vector_area2(),
            Vec3::ZERO
        );
    }

    #[test]
    fn polygon_area_is_translation_invariant_far_from_the_origin() {
        // Georeferenced models sit millions of units from the origin, which is
        // where a shoelace sum about the origin loses precision.
        let offset = Vec3::splat(6_000_000.0);
        let local = Polygon3::new(vec![
            Point3::ZERO,
            Point3::new(1.0, 0.0, 0.0),
            Point3::new(1.0, 1.0, 0.0),
            Point3::new(0.0, 1.0, 0.0),
        ]);
        let far = Polygon3::new(local.vertices.iter().map(|v| *v + offset).collect());
        assert_eq!(far.area(), local.area());
    }
}
