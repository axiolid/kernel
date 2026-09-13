//! Elementary analytic surfaces.

use axiolid_core::{Frame3, Scalar};

/// Infinite plane.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Plane {
    /// Local frame; `z` is the normal.
    pub frame: Frame3,
}

/// Infinite circular cylinder.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Cylinder {
    /// Local frame; `z` is the axis.
    pub frame: Frame3,
    /// Radius.
    pub radius: Scalar,
}

/// Infinite elliptical cylinder.
///
/// Distinct from [`Cylinder`] rather than a special case of it: a circular
/// cylinder's outward normal is its radial direction, and for an ellipse that
/// is only true at the four axis points -- elsewhere the two disagree by up to
/// 53 degrees for a 3:1 ellipse. Anything that assumes radial normals is wrong
/// here, so the type is separate and forces the question.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct EllipticalCylinder {
    /// Local frame; `z` is the axis, `x` and `y` the semi-axis directions.
    pub frame: Frame3,
    /// Semi-axis along local `x`.
    pub semi_axis_x: Scalar,
    /// Semi-axis along local `y`.
    pub semi_axis_y: Scalar,
}

/// Infinite right circular cone.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Cone {
    /// Local frame; `z` is the axis.
    pub frame: Frame3,
    /// Radius at the local origin plane.
    pub radius: Scalar,
    /// Semi-angle in radians.
    pub semi_angle: Scalar,
}

/// Sphere.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Sphere {
    /// Local frame.
    pub frame: Frame3,
    /// Radius.
    pub radius: Scalar,
}

/// Torus.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Torus {
    /// Local frame; `z` is the revolution axis.
    pub frame: Frame3,
    /// Radius from frame origin to tube center.
    pub major_radius: Scalar,
    /// Tube radius.
    pub minor_radius: Scalar,
}
