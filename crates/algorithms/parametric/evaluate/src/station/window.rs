//! Seam-snapping windows (#294, ADR 0082 amendment).
//!
//! An authoring rule may put a station ON a seam of its basis when it lies
//! within the model's precision of one (IFC 4.3 ADD2 8.9.3.48.3: the
//! previous segment's tangent governs a station "within precision limits"
//! of a tangent discontinuity). That precision, up to a millimetre, is far
//! wider than [`ARC_LENGTH_TOLERANCE`], within which every station is on a
//! seam already ([`super::seam`]). A station may therefore carry a WINDOW,
//! a non-negative distance in its basis's station measure:
//! [`snapped_distance2`], [`snapped_distance3`] and
//! [`CompositeBasis::snapped_distance`] return the distance it is read at.
//!
//! # The rule
//!
//! The seams considered are those whose frame may jump (not
//! [`StationSeam::smooth`](super::StationSeam::smooth)): a composite's
//! joints, a polyline's vertices, a B-spline's corner knots, an elevated
//! curve's grade breaks, a banked curve's cant and pivot seams. Each is
//! located by a [`DistanceInterval`] that contains its true distance: a
//! point where the distance is read from stored data, a certified interval
//! no wider than [`WINDOW_CERTIFICATION`] times the window where it is a
//! quadrature (after an offset piece, an ellipse or a B-spline, and a
//! B-spline's corner knot; see [`super::certify`]).
//!
//! - A seam whose interval lies inside `[d - w, d + w]` is within the
//!   window; one whose interval lies wholly outside it is not.
//! - Exactly one seam within the window: the station is read at that
//!   seam's own distance, where the kernel reads the seam, so the seam
//!   rule takes over and its `SeamSide` names the piece. The station's
//!   point is the seam's: the station has moved onto it.
//! - None: the distance is read as it is.
//! - Two or more: refused by a typed [`GeomError::UnsupportedInput`]
//!   naming [`SEAM_WINDOW_AMBIGUOUS`], since the side would not say which
//!   seam.
//! - A seam whose interval straddles the window's edge, so that inside and
//!   outside cannot be told: refused by name, [`SEAM_WINDOW_STRADDLED`].
//!
//! A window of `0` reads the distance as it is, the arc-length tolerance
//! alone; a negative or non-finite window is refused by name.
//!
//! The same test of an interval against a window, [`DistanceInterval::meets`],
//! picks the joints [`CompositeBasis::joint_intervals_within`] (#305) and the
//! curve-evaluation contract's `path_joint_distances_within` return.

use axiolid_contracts::{BackendId, GeomError, GeomResult, Operation};
use axiolid_core::Scalar;
use axiolid_curve::{Curve2, Curve3};

use super::certify::{atomic_seams, CertifiedSeam, DistanceInterval};
use super::composite::StationCurve;
use super::{station_seams2, station_seams3};
use crate::arc_parameter::ARC_LENGTH_TOLERANCE;

#[cfg(doc)]
use super::CompositeBasis;

/// How finely a seam's distance is certified for a window: the certified
/// interval is at most this fraction of the window wide.
pub const WINDOW_CERTIFICATION: Scalar = 1e-3;

/// The input a window holding two seams is refused for.
pub const SEAM_WINDOW_AMBIGUOUS: &str =
    "a station whose seam-snapping window holds two seams of its basis: its seam side does not \
     say which one it lies on";

/// The input a seam straddling a window's edge is refused for.
pub const SEAM_WINDOW_STRADDLED: &str =
    "a station whose seam-snapping window has a seam on its edge: the seam's certified position \
     straddles it, so whether the station lies within the window cannot be told";

fn refused(input: &'static str) -> GeomError {
    GeomError::UnsupportedInput {
        backend: BackendId::new("axiolid-evaluate"),
        operation: Operation::CurveEvaluation,
        input,
    }
}

/// The window, checked: finite and not negative.
///
/// # Errors
///
/// Any other, by name.
pub(crate) fn checked_window(window: Scalar) -> GeomResult<Scalar> {
    if window.is_finite() && window >= 0.0 {
        Ok(window)
    } else {
        Err(GeomError::InvalidInput(format!(
            "station: the seam-snapping window {window} is negative or not finite"
        )))
    }
}

/// The tolerance a window certifies seams to.
pub(crate) fn certification(window: Scalar) -> Scalar {
    window * WINDOW_CERTIFICATION
}

/// The distance a station at `distance` with `window` is read at, among
/// `seams`; see the [module documentation](self#the-rule).
///
/// # Errors
///
/// Two seams within the window, and a seam straddling its edge, by name.
pub(crate) fn snap(
    seams: &[CertifiedSeam],
    distance: Scalar,
    window: Scalar,
) -> GeomResult<Scalar> {
    let (lo, hi) = (distance - window, distance + window);
    let mut inside: Option<Scalar> = None;
    let mut ambiguous = false;
    let mut straddled = false;
    for seam in seams {
        let interval = seam.interval;
        if !interval.meets(&DistanceInterval::new(lo, hi)) {
            continue;
        }
        if interval.lower >= lo && interval.upper <= hi {
            match inside {
                Some(first)
                    if (first - seam.distance).abs()
                        > ARC_LENGTH_TOLERANCE * first.abs().max(1.0) =>
                {
                    ambiguous = true;
                }
                _ => inside = Some(seam.distance),
            }
        } else {
            straddled = true;
        }
    }
    if ambiguous {
        return Err(refused(SEAM_WINDOW_AMBIGUOUS));
    }
    if straddled {
        return Err(refused(SEAM_WINDOW_STRADDLED));
    }
    Ok(inside.unwrap_or(distance))
}

/// The distance a station at `distance` along a 2D curve, with the
/// seam-snapping `window`, is read at; see the
/// [module documentation](self#the-rule).
///
/// # Errors
///
/// A negative or non-finite window, two seams within it, a seam straddling
/// its edge, and a seam whose distance cannot be certified, by name.
pub fn snapped_distance2(curve: &Curve2, distance: Scalar, window: Scalar) -> GeomResult<Scalar> {
    if checked_window(window)? == 0.0 || !distance.is_finite() {
        return Ok(distance);
    }
    let seams = atomic_seams(
        StationCurve::Two(curve),
        station_seams2(curve)?,
        certification(window),
    )?;
    snap(&seams, distance, window)
}

/// [`snapped_distance2`] along a 3D curve.
///
/// # Errors
///
/// As [`snapped_distance2`].
pub fn snapped_distance3(curve: &Curve3, distance: Scalar, window: Scalar) -> GeomResult<Scalar> {
    if checked_window(window)? == 0.0 || !distance.is_finite() {
        return Ok(distance);
    }
    let seams = atomic_seams(
        StationCurve::Three(curve),
        station_seams3(curve)?,
        certification(window),
    )?;
    snap(&seams, distance, window)
}
