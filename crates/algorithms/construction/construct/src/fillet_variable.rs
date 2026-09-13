//! Variable-radius (tapered) fillet blends (ADR 0051).
//!
//! # Why this is not a cone
//!
//! A constant-radius fillet on a vertical edge sweeps a cylinder. Tapering
//! the radius linearly with height does NOT give a right circular cone: the
//! blend centre moves along the bisector as the radius grows, so the axis
//! and the rulings disagree. The surface is a cone in the projective sense
//! -- every ruling passes through one apex -- but an OBLIQUE one, and
//! `Surface::Cone` is right-circular only. Measured obliqueness is far above
//! numerical noise at every interior angle and taper tested (ADR 0051).
//!
//! Each horizontal section is still an exact circular arc, so the surface is
//! a linear loft between two rational quadratic arcs. That is representable
//! exactly as a rational B-spline, and that is what this module builds.

use axiolid_contracts::{GeomError, GeomResult};
use axiolid_core::{Point2, Point3, Scalar};
use axiolid_curve::KnotSpec;
use axiolid_surface::BSplineSurface;

use crate::feature::BlendCorner;

/// A fillet whose radius changes linearly along the extrusion.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct TaperedFillet {
    /// Corner index in the profile ring.
    pub corner: usize,
    /// Radius at the bottom cap.
    pub bottom_radius: Scalar,
    /// Radius at the top cap.
    pub top_radius: Scalar,
}

/// Control points and weight of one rational quadratic arc.
///
/// A circular arc of sweep `angle` is exact as a rational quadratic with
/// control points start, shoulder, end and middle weight `cos(angle/2)`.
/// The shoulder is the intersection of the two end tangents, which sits at
/// `radius / cos(angle/2)` from the centre along the arc's bisector.
fn arc_control_points(
    centre: Point2,
    start: Point2,
    sweep: Scalar,
    height: Scalar,
) -> Option<([Point3; 3], Scalar)> {
    let half = sweep.abs() / 2.0;
    let weight = half.cos();
    // A half-turn or more puts the tangent intersection at infinity, so a
    // single rational quadratic cannot span it.
    if weight <= Scalar::EPSILON {
        return None;
    }
    let radial = start - centre;
    let radius = radial.length();
    if radius <= Scalar::EPSILON {
        return None;
    }

    let (sin, cos) = sweep.sin_cos();
    // Rotate the start radius by the full sweep to reach the end point,
    // which keeps both ends exactly on the circle by construction.
    let end = Point2::new(
        centre.x + radial.x * cos - radial.y * sin,
        centre.y + radial.x * sin + radial.y * cos,
    );
    // Shoulder lies on the bisector of the two radii.
    let mid_angle = sweep / 2.0;
    let (mid_sin, mid_cos) = mid_angle.sin_cos();
    let bisector = Point2::new(
        radial.x * mid_cos - radial.y * mid_sin,
        radial.x * mid_sin + radial.y * mid_cos,
    );
    let shoulder_distance = radius / weight;
    let shoulder = Point2::new(
        centre.x + bisector.x / radius * shoulder_distance,
        centre.y + bisector.y / radius * shoulder_distance,
    );

    Some((
        [
            Point3::new(start.x, start.y, height),
            Point3::new(shoulder.x, shoulder.y, height),
            Point3::new(end.x, end.y, height),
        ],
        weight,
    ))
}

/// The tapered blend surface: a linear loft between two rational arcs.
///
/// `u` runs along the arc (degree 2, rational), `v` along the extrusion
/// (degree 1, exact for a linear taper). Both bounding arcs are exact
/// circles, so every horizontal section of the loft is an exact circle too.
pub fn tapered_blend_surface(
    bottom: &BlendCorner,
    top: &BlendCorner,
    height: Scalar,
) -> GeomResult<BSplineSurface> {
    let (bottom_points, bottom_weight) =
        arc_control_points(bottom.centre, bottom.start, bottom.sweep, 0.0).ok_or_else(|| {
            GeomError::Degenerate("tapered fillet arc spans half a turn or more".to_owned())
        })?;
    let (top_points, top_weight) = arc_control_points(top.centre, top.start, top.sweep, height)
        .ok_or_else(|| {
            GeomError::Degenerate("tapered fillet arc spans half a turn or more".to_owned())
        })?;

    Ok(BSplineSurface {
        u_degree: 2,
        v_degree: 1,
        control_points: vec![
            vec![bottom_points[0], top_points[0]],
            vec![bottom_points[1], top_points[1]],
            vec![bottom_points[2], top_points[2]],
        ],
        u_knots: vec![0.0, 1.0],
        u_multiplicities: vec![3, 3],
        v_knots: vec![0.0, 1.0],
        v_multiplicities: vec![2, 2],
        // The arc weight is constant along v, so the same pair repeats.
        weights: Some(vec![
            vec![1.0, 1.0],
            vec![bottom_weight, top_weight],
            vec![1.0, 1.0],
        ]),
        u_closed: false,
        v_closed: false,
        knot_spec: KnotSpec::Unspecified,
        self_intersect: Some(false),
    })
}
