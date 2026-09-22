//! Bounded two-dimensional primitives.
//!
//! The XY-plane counterparts of [`crate::primitives3`], following the same
//! rules: plain data with public fields, no validation in constructors, and no
//! hidden epsilon. A degenerate value stays representable so the algorithm
//! consuming it can judge it against its own [`Tolerance`].
//!
//! [`Tolerance`]: crate::Tolerance

use crate::{Point2, Scalar, Vec2};

/// Axis-aligned two-dimensional bounding box.
///
/// The XY-plane counterpart of [`Aabb`], with the same empty-absorbs-first
/// behaviour so a fold over points needs no special case for the first one.
///
/// [`Aabb`]: crate::Aabb
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Aabb2 {
    /// Minimum corner.
    pub min: Point2,
    /// Maximum corner.
    pub max: Point2,
}

impl Aabb2 {
    /// An empty box that absorbs the first extended point.
    ///
    /// Inverted on purpose: `min` starts at positive infinity and `max` at
    /// negative infinity, so the first [`extend`] overwrites both rather than
    /// leaving a spurious box around the origin.
    ///
    /// [`extend`]: Aabb2::extend
    pub const fn empty() -> Self {
        Self {
            min: Vec2::splat(Scalar::INFINITY),
            max: Vec2::splat(Scalar::NEG_INFINITY),
        }
    }

    /// Construct a non-empty box containing exactly one point.
    #[inline]
    pub const fn from_point(point: Point2) -> Self {
        Self {
            min: point,
            max: point,
        }
    }

    /// Grow to include a point.
    #[inline]
    pub fn extend(&mut self, point: Point2) {
        self.min = self.min.min(point);
        self.max = self.max.max(point);
    }

    /// Grow to include every point in another box.
    #[inline]
    pub fn union(&mut self, other: &Self) {
        self.min = self.min.min(other.min);
        self.max = self.max.max(other.max);
    }

    /// Whether both corners are finite coordinates.
    #[inline]
    pub fn is_finite(&self) -> bool {
        self.min.is_finite() && self.max.is_finite()
    }

    /// Whether the boxes overlap. Touching counts as overlap.
    #[inline]
    pub fn intersects(&self, other: &Self) -> bool {
        self.min.x <= other.max.x
            && self.max.x >= other.min.x
            && self.min.y <= other.max.y
            && self.max.y >= other.min.y
    }

    /// Whether the box contains a point. The boundary counts as inside.
    #[inline]
    pub fn contains(&self, point: Point2) -> bool {
        point.x >= self.min.x
            && point.x <= self.max.x
            && point.y >= self.min.y
            && point.y <= self.max.y
    }

    /// Whether no point has been added.
    #[inline]
    pub fn is_empty(&self) -> bool {
        self.min.x > self.max.x
    }

    /// Diagonal vector, or zero for an empty box.
    pub fn diagonal(&self) -> Vec2 {
        if self.is_empty() {
            Vec2::ZERO
        } else {
            self.max - self.min
        }
    }

    /// Centre point, or zero for an empty box.
    #[inline]
    pub fn center(&self) -> Point2 {
        if self.is_empty() {
            Vec2::ZERO
        } else {
            (self.min + self.max) * 0.5
        }
    }

    /// Enclosed area, or zero for an empty box.
    pub fn area(&self) -> Scalar {
        let span = self.diagonal();
        span.x * span.y
    }
}

impl Default for Aabb2 {
    fn default() -> Self {
        Self::empty()
    }
}

/// A two-dimensional triangle.
///
/// Corner order sets the winding, which [`signed_area`] reports: positive is
/// counter-clockwise. Degenerate triangles are representable.
///
/// [`signed_area`]: Triangle2::signed_area
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Triangle2 {
    /// First corner.
    pub a: Point2,
    /// Second corner.
    pub b: Point2,
    /// Third corner.
    pub c: Point2,
}

impl Triangle2 {
    /// Construct a triangle without validating it.
    pub const fn new(a: Point2, b: Point2, c: Point2) -> Self {
        Self { a, b, c }
    }

    /// Twice the signed area: positive counter-clockwise, negative clockwise.
    ///
    /// Returned doubled and unscaled because this is the same quantity as the
    /// orientation determinant, so a caller testing which side of `ab` the
    /// point `c` lies on can use it directly without a multiply that would
    /// only add rounding.
    pub fn signed_area2(&self) -> Scalar {
        let ab = self.b - self.a;
        let ac = self.c - self.a;
        ab.perp_dot(ac)
    }

    /// Signed area: positive counter-clockwise, negative clockwise.
    pub fn signed_area(&self) -> Scalar {
        self.signed_area2() * 0.5
    }

    /// Area, always non-negative.
    pub fn area(&self) -> Scalar {
        self.signed_area().abs()
    }

    /// Centroid of the three corners.
    pub fn centroid(&self) -> Point2 {
        (self.a + self.b + self.c) / 3.0
    }
}

/// A rectangle in the XY-plane, possibly rotated.
///
/// Stored as a corner and two edge vectors so the parallelogram property holds
/// by construction, matching [`Rectangle3`]. Perpendicularity is not enforced:
/// a sheared import stays representable and inspectable.
///
/// [`Rectangle3`]: crate::Rectangle3
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Rectangle2 {
    /// Corner the edge vectors start from.
    pub origin: Point2,
    /// First edge vector, from `origin`.
    pub x: Vec2,
    /// Second edge vector, from `origin`.
    pub y: Vec2,
}

impl Rectangle2 {
    /// Construct from a corner and two edge vectors.
    pub const fn new(origin: Point2, x: Vec2, y: Vec2) -> Self {
        Self { origin, x, y }
    }

    /// An axis-aligned rectangle spanning `bounds`.
    ///
    /// Returns `None` for an empty box, because an inverted corner pair has no
    /// meaningful edge vectors and silently producing a zero-sized rectangle
    /// at the origin would hide the emptiness from the caller.
    pub fn from_aabb(bounds: &Aabb2) -> Option<Self> {
        if bounds.is_empty() {
            return None;
        }
        let span = bounds.diagonal();
        Some(Self {
            origin: bounds.min,
            x: Vec2::new(span.x, 0.0),
            y: Vec2::new(0.0, span.y),
        })
    }

    /// The four corners, in order around the boundary.
    pub fn corners(&self) -> [Point2; 4] {
        [
            self.origin,
            self.origin + self.x,
            self.origin + self.x + self.y,
            self.origin + self.y,
        ]
    }

    /// Signed area: positive when the edge vectors are counter-clockwise.
    pub fn signed_area(&self) -> Scalar {
        self.x.perp_dot(self.y)
    }

    /// Area, always non-negative. Equals the parallelogram area when sheared.
    pub fn area(&self) -> Scalar {
        self.signed_area().abs()
    }

    /// Centre point.
    pub fn center(&self) -> Point2 {
        self.origin + (self.x + self.y) * 0.5
    }
}

/// A simple polygon in the XY-plane, without holes.
///
/// The closing edge is implicit, so the last vertex joins the first. Repeating
/// the first vertex at the end creates a zero-length edge rather than closing
/// the ring, which is a common source of degenerate imported data.
///
/// Simplicity -- no self-intersection -- is assumed and not checked: that is a
/// tolerance-dependent judgement, and the types here hold data rather than
/// enforce policy. A polygon set with holes is a different concept; see the
/// planar overlay crate's region type for that.
#[derive(Debug, Clone, PartialEq)]
pub struct Polygon2 {
    /// Boundary vertices in order. The closing edge is implicit.
    pub vertices: Vec<Point2>,
}

impl Polygon2 {
    /// Construct from boundary vertices in order.
    pub const fn new(vertices: Vec<Point2>) -> Self {
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

    /// Twice the signed area by the shoelace formula.
    ///
    /// Summed about the first vertex rather than the origin. The two agree
    /// mathematically, but a polygon in georeferenced coordinates sits far
    /// from the origin, where the terms are large and nearly cancel; rebasing
    /// keeps the terms the size of the polygon instead of the size of the
    /// coordinate system.
    ///
    /// Fewer than three vertices enclose nothing and give zero.
    pub fn signed_area2(&self) -> Scalar {
        if self.vertices.len() < 3 {
            return 0.0;
        }
        let base = self.vertices[0];
        let mut total = 0.0;
        for pair in self.vertices[1..].windows(2) {
            total += (pair[0] - base).perp_dot(pair[1] - base);
        }
        total
    }

    /// Signed area: positive counter-clockwise, negative clockwise.
    pub fn signed_area(&self) -> Scalar {
        self.signed_area2() * 0.5
    }

    /// Area, always non-negative.
    pub fn area(&self) -> Scalar {
        self.signed_area().abs()
    }

    /// Whether the vertex order is counter-clockwise.
    ///
    /// A degenerate polygon has no winding; this reports `false` for it rather
    /// than inventing one.
    pub fn is_counter_clockwise(&self) -> bool {
        self.signed_area2() > 0.0
    }

    /// Axis-aligned bounds of the boundary vertices.
    pub fn bounds(&self) -> Aabb2 {
        let mut bounds = Aabb2::empty();
        for vertex in &self.vertices {
            bounds.extend(*vertex);
        }
        bounds
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn empty_bounds_absorb_the_first_point() {
        let mut bounds = Aabb2::default();
        assert!(bounds.is_empty());
        bounds.extend(Point2::new(3.0, -1.0));
        assert_eq!(bounds.min, bounds.max);
        assert!(!bounds.is_empty());
        assert_eq!(bounds.area(), 0.0);
    }

    #[test]
    fn bounds_intersect_on_touch_and_contain_their_boundary() {
        let mut left = Aabb2::from_point(Point2::ZERO);
        left.extend(Point2::new(1.0, 1.0));
        let mut right = Aabb2::from_point(Point2::new(1.0, 0.0));
        right.extend(Point2::new(2.0, 1.0));
        assert!(left.intersects(&right), "touching boxes overlap");
        assert!(left.contains(Point2::new(1.0, 1.0)), "boundary is inside");

        let mut away = Aabb2::from_point(Point2::new(5.0, 5.0));
        away.extend(Point2::new(6.0, 6.0));
        assert!(!left.intersects(&away));
    }

    #[test]
    fn triangle_signed_area_carries_winding_but_area_does_not() {
        let ccw = Triangle2::new(Point2::ZERO, Point2::new(4.0, 0.0), Point2::new(0.0, 2.0));
        assert_eq!(ccw.signed_area(), 4.0);
        assert!(ccw.signed_area2() > 0.0);

        let cw = Triangle2::new(ccw.a, ccw.c, ccw.b);
        assert_eq!(cw.signed_area(), -4.0);
        assert_eq!(cw.area(), ccw.area());
    }

    #[test]
    fn collinear_triangle_has_zero_area() {
        let degenerate = Triangle2::new(Point2::ZERO, Point2::new(1.0, 1.0), Point2::new(3.0, 3.0));
        assert_eq!(degenerate.signed_area2(), 0.0);
        assert_eq!(degenerate.area(), 0.0);
    }

    #[test]
    fn rectangle_from_bounds_matches_the_box_it_came_from() {
        let mut bounds = Aabb2::from_point(Point2::new(1.0, 2.0));
        bounds.extend(Point2::new(4.0, 6.0));
        let rectangle = Rectangle2::from_aabb(&bounds).expect("non-empty bounds");
        assert_eq!(rectangle.area(), bounds.area());
        assert_eq!(rectangle.center(), bounds.center());
        assert_eq!(rectangle.corners()[2], bounds.max);
    }

    #[test]
    fn rectangle_from_empty_bounds_is_refused_rather_than_zero_sized() {
        assert!(Rectangle2::from_aabb(&Aabb2::empty()).is_none());
    }

    #[test]
    fn rotated_rectangle_keeps_its_area() {
        // A 45-degree rotated unit square: axis-aligned bounds would overstate
        // the area, which is why Rectangle2 is not stored as an Aabb2.
        let diagonal = Vec2::new(1.0, 1.0);
        let rectangle = Rectangle2::new(Point2::ZERO, diagonal, Vec2::new(-1.0, 1.0));
        assert_eq!(rectangle.area(), 2.0);
    }

    #[test]
    fn polygon_winding_flips_with_vertex_order() {
        let square = Polygon2::new(vec![
            Point2::ZERO,
            Point2::new(2.0, 0.0),
            Point2::new(2.0, 2.0),
            Point2::new(0.0, 2.0),
        ]);
        assert_eq!(square.area(), 4.0);
        assert!(square.is_counter_clockwise());

        let mut reversed = square.vertices.clone();
        reversed.reverse();
        let reversed = Polygon2::new(reversed);
        assert_eq!(reversed.area(), square.area());
        assert!(!reversed.is_counter_clockwise());
        assert_eq!(reversed.signed_area(), -square.signed_area());
    }

    #[test]
    fn polygon_with_fewer_than_three_vertices_encloses_nothing() {
        assert_eq!(Polygon2::new(Vec::new()).signed_area2(), 0.0);
        assert_eq!(
            Polygon2::new(vec![Point2::ZERO, Point2::new(1.0, 0.0)]).signed_area2(),
            0.0
        );
        assert!(!Polygon2::new(Vec::new()).is_counter_clockwise());
    }

    #[test]
    fn polygon_area_is_translation_invariant_far_from_the_origin() {
        // Site coordinates routinely sit millions of units out. Summing the
        // shoelace terms about the origin there loses significant digits.
        let offset = Vec2::splat(6_000_000.0);
        let local = Polygon2::new(vec![
            Point2::ZERO,
            Point2::new(1.0, 0.0),
            Point2::new(1.0, 1.0),
            Point2::new(0.0, 1.0),
        ]);
        let far = Polygon2::new(local.vertices.iter().map(|v| *v + offset).collect());
        assert_eq!(far.area(), local.area());
    }

    #[test]
    fn polygon_bounds_cover_every_vertex() {
        let polygon = Polygon2::new(vec![
            Point2::new(-1.0, 4.0),
            Point2::new(3.0, -2.0),
            Point2::new(0.0, 0.0),
        ]);
        let bounds = polygon.bounds();
        assert_eq!(bounds.min, Point2::new(-1.0, -2.0));
        assert_eq!(bounds.max, Point2::new(3.0, 4.0));
        for vertex in &polygon.vertices {
            assert!(bounds.contains(*vertex));
        }
    }
}
