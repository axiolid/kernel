//! Faces of a triangle mesh selected by their normal (#417).
//!
//! Takeoff quantities need the area of the downward faces (a floor), or of
//! the faces normal to a wall with their extent along it. A face is selected
//! when the angle between its normal and a cone's axis is within the cone's
//! half-angle plus the tolerance's angular slack. The angle is
//! `atan2(|N x axis|, N . axis)` of the triangle's raw doubled-area vector
//! `N`, which stays accurate at angles `acos` of a normalised dot product
//! cannot resolve.
//!
//! The mesh must be surface-usable with no degenerate triangle, exactly as
//! for [`crate::surface_properties`], so a normal is always defined. Areas
//! are summed in ascending triangle order about the mesh's first vertex, as
//! `surface_properties` does, so selecting every face reproduces its area bit
//! for bit.

use core::fmt;

use axiolid_core::{Interval, Tolerance, Vec3};
use axiolid_mesh::{audit_mesh, MeshHealth, TriangleMeshView};

use crate::mesh::triangle_area;

/// A direction and an angular reach about it.
///
/// The axis is stored normalised, so the stored value does not depend on the
/// length the caller passed.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct NormalCone {
    axis: Vec3,
    half_angle: f64,
}

impl NormalCone {
    /// A cone about `axis` reaching `half_angle` radians.
    ///
    /// # Errors
    /// [`FaceSelectionError::InvalidAxis`] when `axis` is not finite or has
    /// no direction; [`FaceSelectionError::InvalidHalfAngle`] when
    /// `half_angle` is not finite or outside `[0, pi]`.
    pub fn new(axis: Vec3, half_angle: f64) -> Result<Self, FaceSelectionError> {
        let axis = unit(axis).ok_or(FaceSelectionError::InvalidAxis)?;
        if !(0.0..=core::f64::consts::PI).contains(&half_angle) {
            return Err(FaceSelectionError::InvalidHalfAngle);
        }
        Ok(Self { axis, half_angle })
    }

    /// The unit axis.
    pub fn axis(self) -> Vec3 {
        self.axis
    }

    /// The half-angle in radians.
    pub fn half_angle(self) -> f64 {
        self.half_angle
    }
}

/// Triangles chosen by [`select_by_normal`].
#[non_exhaustive]
#[derive(Debug, Clone, PartialEq)]
pub struct FaceSelection {
    /// Selected triangle indices, ascending.
    pub triangles: Vec<usize>,
    /// Total area of the selected triangles; `0.0` when none is selected.
    pub area: f64,
}

/// Why a selection or an extent was refused.
#[non_exhaustive]
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum FaceSelectionError {
    /// The cone axis is not finite or has no direction.
    InvalidAxis,
    /// The half-angle is not finite or outside `[0, pi]`.
    InvalidHalfAngle,
    /// The direction to measure along is not finite or has no direction.
    InvalidDirection,
    /// The mesh is not surface-usable or has a degenerate triangle.
    MeshNotSurfaceUsable(MeshHealth),
    /// A listed triangle does not exist.
    TriangleOutOfRange {
        /// The listed index.
        triangle: usize,
    },
    /// A corner of a triangle does not address a position.
    IndexOutOfRange {
        /// The triangle holding the corner.
        triangle: usize,
    },
    /// A triangle's normal or area overflows. (A total of finite areas cannot:
    /// a triangle's area is below `1e154`.)
    NonFiniteNormal {
        /// The triangle concerned.
        triangle: usize,
    },
    /// A listed triangle has a non-finite position.
    NonFinitePosition {
        /// The triangle concerned.
        triangle: usize,
    },
}

impl fmt::Display for FaceSelectionError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        f.write_str(match self {
            Self::InvalidAxis => "cone axis is not a finite direction",
            Self::InvalidHalfAngle => "cone half-angle is not finite within [0, pi]",
            Self::InvalidDirection => "direction is not a finite direction",
            Self::MeshNotSurfaceUsable(_) => "mesh is not suitable for face selection",
            Self::TriangleOutOfRange { .. } => "triangle index is out of range",
            Self::IndexOutOfRange { .. } => "triangle corner does not address a position",
            Self::NonFiniteNormal { .. } => "triangle normal or area is not finite",
            Self::NonFinitePosition { .. } => "triangle has a non-finite position",
        })
    }
}

impl std::error::Error for FaceSelectionError {}

/// `direction` as a unit vector. It is scaled by its largest component first,
/// so squaring neither underflows (`1e-170`) nor overflows (`1e200`).
fn unit(direction: Vec3) -> Option<Vec3> {
    let largest = direction
        .x
        .abs()
        .max(direction.y.abs())
        .max(direction.z.abs());
    if !largest.is_finite() || largest == 0.0 {
        return None;
    }
    let scaled = direction / largest;
    let unit = scaled / scaled.length();
    unit.is_finite().then_some(unit)
}

/// Select the triangles whose normal lies within `cone`.
///
/// `tolerance.linear()` decides degenerate triangles as
/// [`crate::surface_properties`] does (any degenerate triangle refuses the
/// mesh); `tolerance.angular()` is slack added to the half-angle. A triangle
/// is selected when `atan2(|N x axis|, N . axis) <= half_angle + angular`.
/// A floor tilted by `1e-12` is therefore selected at half-angle zero with
/// [`Tolerance::METRE`] but not with a zero angular tolerance, and an exactly
/// level floor is selected either way. A half-angle of pi selects every
/// triangle.
///
/// # Errors
/// [`FaceSelectionError::MeshNotSurfaceUsable`] when the audit finds the mesh
/// unusable or degenerate, [`FaceSelectionError::NonFiniteNormal`] when a
/// triangle's normal or area is not finite.
pub fn select_by_normal<M: TriangleMeshView + ?Sized>(
    mesh: &M,
    cone: NormalCone,
    tolerance: Tolerance,
) -> Result<FaceSelection, FaceSelectionError> {
    let health = audit_mesh(mesh, tolerance);
    if !health.is_surface_usable() || health.degenerate_triangles != 0 {
        return Err(FaceSelectionError::MeshNotSurfaceUsable(health));
    }
    let reach = cone.half_angle + tolerance.angular();
    let base = mesh.position(0);
    let (mut triangles, mut area) = (Vec::new(), 0.0f64);
    for i in 0..mesh.triangle_count() {
        let [a, b, c] = mesh.triangle(i).map(|j| mesh.position(j as usize) - base);
        let normal = (b - a).cross(c - a);
        // Decided before the angle: an overflowed normal or area would make
        // the angle wrong, not merely non-finite.
        let weight = triangle_area([a, b, c]);
        if !normal.is_finite() || !weight.is_finite() {
            return Err(FaceSelectionError::NonFiniteNormal { triangle: i });
        }
        let angle = normal
            .cross(cone.axis)
            .length()
            .atan2(normal.dot(cone.axis));
        if angle <= reach {
            triangles.push(i);
            area += weight;
        }
    }
    Ok(FaceSelection { triangles, area })
}

/// Extent of the listed triangles along `along`: the least and greatest of
/// `p . u` over every corner, `u = along / |along|`, as `start` and `end`.
///
/// Duplicate and unordered triangle indices are fine. An empty list gives
/// `None`. The triangles need not form a surface; only their corners are read.
///
/// # Errors
/// [`FaceSelectionError::InvalidDirection`] for an unusable `along`;
/// [`FaceSelectionError::TriangleOutOfRange`],
/// [`FaceSelectionError::IndexOutOfRange`] or
/// [`FaceSelectionError::NonFinitePosition`] naming the first offending
/// listed triangle.
pub fn extent_along<M: TriangleMeshView + ?Sized>(
    mesh: &M,
    triangles: &[usize],
    along: Vec3,
) -> Result<Option<Interval>, FaceSelectionError> {
    let direction = unit(along).ok_or(FaceSelectionError::InvalidDirection)?;
    let mut extent: Option<(f64, f64)> = None;
    for &triangle in triangles {
        if triangle >= mesh.triangle_count() {
            return Err(FaceSelectionError::TriangleOutOfRange { triangle });
        }
        for corner in mesh.triangle(triangle) {
            let index = usize::try_from(corner)
                .ok()
                .filter(|&index| index < mesh.position_count())
                .ok_or(FaceSelectionError::IndexOutOfRange { triangle })?;
            let position = mesh.position(index);
            let height = position.dot(direction);
            if !position.is_finite() || !height.is_finite() {
                return Err(FaceSelectionError::NonFinitePosition { triangle });
            }
            extent = Some(match extent {
                Some((low, high)) => (low.min(height), high.max(height)),
                None => (height, height),
            });
        }
    }
    Ok(extent.map(|(low, high)| Interval::new(low, high)))
}
