//! Bounded analytic solid primitives.

use axiolid_core::Scalar;

/// Exact primitive solid; tessellation is a separate operation.
#[non_exhaustive]
#[derive(Debug, Clone, Copy, PartialEq)]
#[cfg_attr(
    feature = "serde",
    derive(serde::Serialize, serde::Deserialize),
    serde(deny_unknown_fields)
)]
pub enum Primitive {
    /// Axis-aligned box in local coordinates.
    Block {
        /// Extent along x.
        x: Scalar,
        /// Extent along y.
        y: Scalar,
        /// Extent along z.
        z: Scalar,
    },
    /// Sphere centered at the local origin.
    Sphere {
        /// Radius.
        radius: Scalar,
    },
    /// Cylinder along local +z.
    Cylinder {
        /// Radius.
        radius: Scalar,
        /// Height.
        height: Scalar,
    },
    /// Cone along local +z.
    Cone {
        /// Base radius.
        radius: Scalar,
        /// Height.
        height: Scalar,
    },
    /// Rectangular pyramid along local +z.
    Pyramid {
        /// Base extent along x.
        x: Scalar,
        /// Base extent along y.
        y: Scalar,
        /// Height.
        height: Scalar,
    },
    /// Ring torus centred on the local origin, about local +z.
    ///
    /// The tube's centre circle has radius `major_radius` in the plane
    /// z = 0. Only a ring torus, `0 < minor_radius < major_radius`, is a
    /// valid solid: at `minor == major` (a horn torus) the tube pinches the
    /// axis to a point, and past it (a spindle torus) the tube overlaps
    /// itself, so neither bounds a two-manifold solid.
    Torus {
        /// Radius of the tube's centre circle.
        major_radius: Scalar,
        /// Radius of the tube.
        minor_radius: Scalar,
    },
    /// Wedge: a box whose top face is narrowed or shifted, as OCCT's
    /// `BRepPrimAPI_MakeWedge(dx, dy, dz, xmin, zmin, xmax, zmax)` with
    /// OCCT's y (the height) read as local +z here.
    ///
    /// The base is the rectangle `[0, x] x [0, y]` at z = 0 and the top is
    /// `[top_x_min, top_x_max] x [top_y_min, top_y_max]` at z = `height`.
    /// Every face is planar and the solid is the convex hull of the two
    /// rectangles. The top may collapse to a segment (`top_x_min ==
    /// top_x_max` or `top_y_min == top_y_max`, a wedge proper) or to a point
    /// (a pyramid with its apex anywhere); OCCT's `ltx` form is `top_x_min =
    /// 0`, `top_x_max = ltx`, `top_y_min = 0`, `top_y_max = y`.
    Wedge {
        /// Base extent along x.
        x: Scalar,
        /// Base extent along y.
        y: Scalar,
        /// Height along z.
        height: Scalar,
        /// Top face's least x.
        top_x_min: Scalar,
        /// Top face's greatest x.
        top_x_max: Scalar,
        /// Top face's least y.
        top_y_min: Scalar,
        /// Top face's greatest y.
        top_y_max: Scalar,
    },
}
