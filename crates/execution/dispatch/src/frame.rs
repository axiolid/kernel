//! View-frame validation, shared by the registries that take a frame.

use axiolid_contracts::{GeomError, GeomResult};
use axiolid_core::{Frame3, SpaceFrame, Tolerance};

/// Check that `frame` is a finite right-handed orthonormal frame.
///
/// `what` names the frame in the refusal ("section frame", "hidden-line view
/// frame"). Structural frame validity cannot be relaxed into zero/scaled axes
/// by a coarse geometric tolerance, so the caller's angular tolerance is
/// capped here before it reaches the shared validator. Tighter caller
/// tolerances still apply; looser ones do not.
pub(crate) fn validate_frame(frame: Frame3, angular_tolerance: f64, what: &str) -> GeomResult<()> {
    let limit = angular_tolerance.clamp(128.0 * f64::EPSILON, 1.0e-6);
    let tolerance = Tolerance::new(limit, limit)
        .map_err(|_| GeomError::InvalidInput(format!("{what} tolerance is invalid")))?;
    SpaceFrame::new(frame.origin, frame.x, frame.y, frame.z, tolerance)
        .map(|_| ())
        .map_err(|error| GeomError::InvalidInput(format!("{what} is invalid: {error}")))
}
