#![forbid(unsafe_code)]

//! Exact, format-neutral surface representations and evaluation contracts.
//!
//! Bounded, swept, offset, and curve-on-surface relationships are nodes in
//! `axiolid-model`; this crate contains atomic surface data only.

pub mod elementary;
pub mod evaluate;
pub mod spline;

pub use elementary::{Cone, Cylinder, EllipticalCylinder, Plane, Sphere, Torus};
pub use evaluate::SurfaceEvaluator;
pub use spline::BSplineSurface;

/// Atomic three-dimensional surface values.
#[non_exhaustive]
#[derive(Debug, Clone, PartialEq)]
pub enum Surface {
    /// Plane.
    Plane(Plane),
    /// Circular cylinder.
    Cylinder(Cylinder),
    /// Elliptical cylinder.
    EllipticalCylinder(EllipticalCylinder),
    /// Right circular cone.
    Cone(Cone),
    /// Sphere.
    Sphere(Sphere),
    /// Torus.
    Torus(Torus),
    /// Polynomial or rational B-spline surface.
    BSpline(BSplineSurface),
}
