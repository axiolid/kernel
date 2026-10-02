//! Exact swept disks along a straight segment and along a circular arc
//! (#223).
//!
//! # What a swept disk is here
//!
//! A disk of radius `r` (optionally with a concentric bore of radius
//! `r_i`) moved along a directrix, held perpendicular to it. Along a
//! straight segment that is a circular cylinder capped by two discs: an
//! exact extrusion of the disk. Along an arc of radius `R` about an axis it
//! is a torus wedge capped by two discs, or a whole torus for a full
//! circle: an exact revolution of the disk about the arc's axis. Both are
//! built in a local frame by the exact extrusion and revolution, then
//! placed with [`ExactBRep::transformed`], so the walls are genuine
//! cylinders and tori and the caps genuine planes -- never a tessellation.
//!
//! # What is refused
//!
//! A disk whose radius reaches the arc's axis (`r >= R`) would sweep
//! through itself; a bore not strictly inside the disk is not a ring; an
//! arc beyond a full turn sweeps through itself. All are refused by name
//! or as invalid input. A directrix with corners, or a curved directrix
//! other than an arc, is the compiler's to refuse: these functions take a
//! single segment or arc.

use std::f64::consts::TAU;

use axiolid_brep::{ExactBRep, TransformError};
use axiolid_contracts::{GeomError, GeomResult, Operation};
use axiolid_core::{Interval, Mat3, Point3, Scalar, Tolerance, Transform3, Vec3};
use axiolid_curve::Circle3;
use axiolid_profile::{CircleProfile, Profile};

use crate::extrude_exact::extrude_profile_exact;
use crate::revolve_exact::revolve_profile_exact;
use crate::BACKEND_ID;

fn unsupported(input: &'static str) -> GeomError {
    GeomError::UnsupportedInput {
        backend: BACKEND_ID,
        operation: Operation::Sweep,
        input,
    }
}

/// The disk as a profile in its local plane, centred at the origin.
fn disk(radius: Scalar, inner_radius: Option<Scalar>) -> GeomResult<Profile> {
    if !radius.is_finite() || radius <= 0.0 {
        return Err(GeomError::InvalidInput(format!(
            "swept disk radius must be positive and finite, got {radius}"
        )));
    }
    match inner_radius {
        // A bore of radius zero is no bore.
        None | Some(0.0) => Ok(Profile::Circle(CircleProfile {
            radius,
            thickness: None,
        })),
        Some(inner) if inner.is_finite() && inner > 0.0 && inner < radius => {
            Ok(Profile::Circle(CircleProfile {
                radius,
                thickness: Some(radius - inner),
            }))
        }
        Some(inner) => Err(GeomError::InvalidInput(format!(
            "swept disk inner radius must lie in [0, {radius}), got {inner}"
        ))),
    }
}

/// Place a B-rep built in a local frame; the frames here are orthonormal by
/// construction, so a refusal names the support that could not move.
pub(crate) fn place(solid: &ExactBRep, transform: &Transform3) -> GeomResult<ExactBRep> {
    solid.transformed(transform).map_err(|error| match error {
        TransformError::Unsupported(what) => unsupported(what),
        other => GeomError::InvalidInput(format!("swept disk placement: {other}")),
    })
}

/// A unit vector perpendicular to the unit `d`, chosen from the coordinate
/// axis least aligned with it so the choice is deterministic and well
/// conditioned.
fn perpendicular(d: Vec3) -> Vec3 {
    let a = d.abs();
    let helper = if a.x <= a.y && a.x <= a.z {
        Vec3::X
    } else if a.y <= a.z {
        Vec3::Y
    } else {
        Vec3::Z
    };
    (helper - d * helper.dot(d)).normalize()
}

/// A disk of `radius`, with an optional bore of `inner_radius`, swept along
/// the straight segment from `start` to `end`: an exact capped cylinder
/// (#223).
///
/// # Errors
///
/// Invalid input for a non-positive radius, a bore not strictly inside the
/// disk, and a segment shorter than the tolerance or not finite.
pub fn swept_disk_along_line_exact(
    start: Point3,
    end: Point3,
    radius: Scalar,
    inner_radius: Option<Scalar>,
    tolerance: Tolerance,
) -> GeomResult<ExactBRep> {
    let profile = disk(radius, inner_radius)?;
    let along = end - start;
    let length = along.length();
    if !length.is_finite() || length <= tolerance.linear() {
        return Err(GeomError::InvalidInput(format!(
            "swept disk segment must be finite and longer than the tolerance, got {length}"
        )));
    }
    let z = along / length;
    let x = perpendicular(z);
    let y = z.cross(x);
    // A hollow disk extrudes through the composite path, which builds the
    // bore as an inner wall; the dedicated circle path takes a full disk.
    let profile = match profile {
        hollow @ Profile::Circle(CircleProfile {
            thickness: Some(_), ..
        }) => Profile::Composite(vec![hollow]),
        full => full,
    };
    let local = extrude_profile_exact(&profile, Vec3::Z, length, tolerance)?;
    place(
        &local,
        &Transform3::from_mat3_translation(Mat3::from_cols(x, y, z), start),
    )
}

/// A disk of `radius`, with an optional bore of `inner_radius`, swept along
/// `arc` over the angle span `span` (radians, start to end, in `arc`'s own
/// parameter): an exact torus wedge capped by two discs, or a whole torus
/// for a span of a full turn (#223).
///
/// The disk is perpendicular to the arc at every point; the arc's frame
/// must be orthonormal.
///
/// # Errors
///
/// Invalid input for a non-positive radius, a bore not strictly inside the
/// disk, a zero or non-finite span, and a non-orthonormal arc frame;
/// refused by name for a span beyond a full turn and a disk reaching the
/// arc's axis.
pub fn swept_disk_along_arc_exact(
    arc: &Circle3,
    span: Interval,
    radius: Scalar,
    inner_radius: Option<Scalar>,
    tolerance: Tolerance,
) -> GeomResult<ExactBRep> {
    let profile = disk(radius, inner_radius)?;
    let bend = arc.radius;
    if !bend.is_finite() || bend <= 0.0 {
        return Err(GeomError::InvalidInput(format!(
            "swept disk arc radius must be positive and finite, got {bend}"
        )));
    }
    let angle = span.end - span.start;
    if !angle.is_finite() || angle == 0.0 {
        return Err(GeomError::InvalidInput(format!(
            "swept disk arc span must be finite and non-zero, got {angle}"
        )));
    }
    if angle.abs() > TAU + tolerance.linear() {
        return Err(unsupported("swept disk along an arc beyond a full turn"));
    }
    if radius >= bend - tolerance.linear() {
        return Err(unsupported("swept disk reaching its arc's axis"));
    }
    let frame = &arc.frame;
    let axes = [frame.x, frame.y, frame.z];
    let orthonormal = axes.iter().all(|a| a.is_finite())
        && frame.origin.is_finite()
        && (0..3).all(|i| {
            (i..3).all(|j| {
                let want = if i == j { 1.0 } else { 0.0 };
                (axes[i].dot(axes[j]) - want).abs() <= 1e-12
            })
        })
        && (frame.x.cross(frame.y) - frame.z).length() <= 1e-12;
    if !orthonormal {
        return Err(GeomError::InvalidInput(
            "swept disk arc frame must be a right-handed orthonormal frame".to_owned(),
        ));
    }

    // Local frame: the disk in the plane z = 0 centred at the origin, the
    // axis the line x = -R along +y. A turn by `phi` about +y carries the
    // centre to `A + R (cos phi, 0, -sin phi)`.
    let local = revolve_profile_exact(
        &profile,
        Point3::new(-bend, 0.0, 0.0),
        Vec3::Y,
        angle,
        tolerance,
    )?;
    // World: local x to the radial direction at the start, local y to the
    // arc's normal, local z against the arc's tangent there, so a turn by
    // `phi` about the normal follows the arc from `span.start`.
    let (sin, cos) = span.start.sin_cos();
    let radial = frame.x * cos + frame.y * sin;
    let tangent = frame.z.cross(radial);
    let rotation = Mat3::from_cols(radial, frame.z, -tangent);
    let translation = frame.origin + radial * bend;
    place(
        &local,
        &Transform3::from_mat3_translation(rotation, translation),
    )
}
