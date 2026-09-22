//! Coordinate, direction, transform, and analytic support types.

use crate::Scalar;

/// Double-precision two-dimensional vector.
pub type Vec2 = glam::DVec2;
/// Double-precision three-dimensional vector.
pub type Vec3 = glam::DVec3;
/// A semantic alias used when a value is a two-dimensional position.
pub type Point2 = Vec2;
/// A semantic alias used when a value is a three-dimensional position.
pub type Point3 = Vec3;
/// Double-precision 3x3 matrix.
pub type Mat3 = glam::DMat3;
/// Double-precision 2D affine transform.
pub type Transform2 = glam::DAffine2;
/// Double-precision affine transform.
///
/// Affine transforms cannot represent perspective: the implicit bottom row is
/// always `[0, 0, 0, 1]`. Use [`Mat4`] when a projection is needed.
pub type Transform3 = glam::DAffine3;
/// Double-precision general 4x4 matrix, including projective transforms.
///
/// Distinct from [`Transform3`], which is affine and therefore cannot express
/// perspective. This type carries a real fourth row, so it can represent a
/// projection whose `w` varies per point -- which is exactly what makes a
/// homogeneous divide meaningful.
///
/// Prefer [`Transform3`] for the ordinary placement, scaling, and rotation of
/// geometry: it is smaller, composes faster, and its inverse is always
/// well-defined. Reach for `Mat4` only when the transform genuinely projects.
///
/// Note that the previous release aliased `Mat4` to [`Transform3`], so it was
/// affine despite the name. Code that wants the old meaning should name
/// [`Transform3`] explicitly; the mismatch in API surface makes that a compile
/// error rather than a silent change in behaviour.
pub type Mat4 = glam::DMat4;

/// Right-handed 2D local frame. Algorithms validate orthonormality explicitly.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Frame2 {
    /// Local origin.
    pub origin: Point2,
    /// Local x axis.
    pub x: Vec2,
    /// Local y axis.
    pub y: Vec2,
}

/// Right-handed 3D local frame. Dirty imported frames remain representable.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Frame3 {
    /// Local origin.
    pub origin: Point3,
    /// Local x axis.
    pub x: Vec3,
    /// Local y axis.
    pub y: Vec3,
    /// Local z axis.
    pub z: Vec3,
}

/// A finite parameter interval. The endpoint order carries orientation.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Interval {
    /// Start parameter.
    pub start: Scalar,
    /// End parameter.
    pub end: Scalar,
}

impl Interval {
    /// Unit parameter interval.
    pub const UNIT: Self = Self {
        start: 0.0,
        end: 1.0,
    };

    /// Construct an oriented interval without sorting its endpoints.
    pub const fn new(start: Scalar, end: Scalar) -> Self {
        Self { start, end }
    }

    /// Absolute parameter span.
    pub fn length(self) -> Scalar {
        (self.end - self.start).abs()
    }
}

/// A plane represented by an origin and unit-normal candidate.
///
/// Adapters may construct dirty input. Algorithms validate normalization using
/// the operation's tolerance instead of hiding a global epsilon here.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Plane3 {
    /// Point on the plane.
    pub origin: Point3,
    /// Expected outward normal.
    pub normal: Vec3,
}

/// A parametric three-dimensional ray.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Ray3 {
    /// Ray start.
    pub origin: Point3,
    /// Ray direction. It need not be normalized at the storage boundary.
    pub direction: Vec3,
}

#[cfg(test)]
mod tests {
    use super::*;

    /// The reason `Mat4` stopped being an alias for [`Transform3`].
    ///
    /// A perspective projection makes `w` vary per point, so recovering a
    /// cartesian coordinate requires dividing by it. An affine transform
    /// cannot express that at all, which is why the old alias could not close
    /// this gap no matter how it was called.
    #[test]
    fn projective_matrix_divides_by_w_and_affine_cannot() {
        // Looking down -Z with a 90 degree vertical field of view, so a point
        // at depth d has its height scaled by 1/d.
        let projection = Mat4::perspective_rh(std::f64::consts::FRAC_PI_2, 1.0, 1.0, 100.0);

        let near = projection.project_point3(Point3::new(0.0, 1.0, -2.0));
        let far = projection.project_point3(Point3::new(0.0, 1.0, -4.0));

        // Same world height, twice the depth: the projected height halves.
        // That ratio is the homogeneous divide doing its job.
        assert!(
            (near.y / far.y - 2.0).abs() < 1e-12,
            "near {near:?} far {far:?}"
        );

        // The fourth row is what carries it. An affine transform's implicit
        // bottom row is [0, 0, 0, 1], so w is constant and no such ratio can
        // arise: identical input keeps its height at both depths.
        let affine = Transform3::IDENTITY;
        let near_affine = affine.transform_point3(Point3::new(0.0, 1.0, -2.0));
        let far_affine = affine.transform_point3(Point3::new(0.0, 1.0, -4.0));
        assert_eq!(near_affine.y, far_affine.y);
    }

    /// `Mat4` and `Transform3` are now genuinely different types.
    ///
    /// Asserted so a future "simplification" back to an alias fails here
    /// rather than silently removing projection support again.
    #[test]
    fn a_projective_matrix_round_trips_through_its_affine_subset() {
        let affine = Transform3::from_translation(Vec3::new(1.0, 2.0, 3.0));
        let promoted = Mat4::from(affine);
        let point = Point3::new(0.5, -0.5, 2.0);
        // Promoting an affine transform must not change what it does.
        assert_eq!(
            promoted.project_point3(point),
            affine.transform_point3(point)
        );
    }
}
