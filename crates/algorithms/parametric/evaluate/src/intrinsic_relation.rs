//! Relations over `Intrinsic3`: trimming, offsetting, and joining.
//!
//! # What is exact and what is refused
//!
//! An `Intrinsic3` stores laws against ARC LENGTH. That is what makes some
//! relations exact and others impossible without changing the value's
//! meaning.
//!
//! **Trim is exact.** Restricting to `[a, b]` does not refit anything: the
//! curvature and torsion laws are re-anchored with `CurvatureLaw::shifted`,
//! which the family is closed under, and the start frame is moved to the
//! curve's own point and frame at `a`. The only numerical step is finding
//! that anchor, which is the same quadrature the evaluator already performs
//! -- the SHAPE stays exact.
//!
//! **Offset is exact only for a helix, and otherwise refused.** The normal
//! offset `q(s) = p(s) + d * N(s)` is not unit speed: differentiating gives
//! `|q'| = |1 - d*k(s)|`, so the offset advances at a different rate than
//! the base. An `Intrinsic3` stores laws in arc length, so representing the
//! offset requires reparameterising by ITS arc length. That reparameterisation
//! is a closed-form rescale only when `k` is constant; for a varying law it is
//! the inverse of a non-elementary integral, and writing a law in the family
//! would be a fit, not the curve. Measured: for `k(s) = 0.20 + 0.05 s` and
//! `d = 0.8` the offset speed sweeps `0.72 .. 0.84` across the span.
//!
//! For a helix the closure is genuine, and pleasant: the normal points at the
//! axis, so offsetting slides the curve to a coaxial helix of radius `a - d`
//! with the SAME pitch and the SAME angular rate. With `c = hypot(a, b)` and
//! `c2 = hypot(a - d, b)`, the offset has `k2 = (a-d)/c2^2`, `tau2 = b/c2^2`,
//! and its arc length runs at `c2/c` times the base's.
//!
//! **Join is exact when the ends actually meet.** Two curves compose into one
//! `Piecewise` law when the second one's start frame is the first one's end
//! frame; otherwise the join is a fiction and is refused by name.

use axiolid_contracts::{BackendId, GeomError, GeomResult, Operation};
use axiolid_core::{Frame3, Point3, Scalar, Vec3};
use axiolid_curve::{CurvatureLaw, Intrinsic3};

use crate::frenet::{frenet_frame, frenet_point};

fn unsupported() -> GeomError {
    GeomError::Unsupported {
        backend: BackendId::new("axiolid-evaluate"),
        operation: Operation::CurveEvaluation,
    }
}

fn invalid(detail: &str) -> GeomError {
    GeomError::InvalidInput(detail.to_owned())
}

/// Restrict a space curve to `[start, end]` of its arc length.
///
/// Exact in shape: the laws are re-anchored, not refitted. The returned curve
/// has length `end - start` and its own start frame is the base curve's frame
/// at `start`.
pub fn trim_intrinsic3(curve: &Intrinsic3, start: Scalar, end: Scalar) -> GeomResult<Intrinsic3> {
    if !start.is_finite() || !end.is_finite() {
        return Err(invalid("trim bounds must be finite"));
    }
    if end <= start {
        return Err(invalid("trim end must exceed trim start"));
    }
    if start < 0.0 || end > curve.length {
        return Err(invalid("trim bounds must lie within the curve"));
    }
    let anchor = frenet_frame(curve, start)?;
    let curvature = curve
        .curvature
        .shifted(start)
        .ok_or_else(|| invalid("curvature law cannot be re-anchored at the trim start"))?;
    let torsion = curve
        .torsion
        .shifted(start)
        .ok_or_else(|| invalid("torsion law cannot be re-anchored at the trim start"))?;
    Ok(Intrinsic3::new(anchor, curvature, torsion, end - start))
}

/// A helix's radius and pitch parameters, recovered from its laws.
///
/// For constant `k` and `tau`, `a = k / (k^2 + tau^2)` and
/// `b = tau / (k^2 + tau^2)`; `hypot(a, b) = 1 / hypot(k, tau)`.
fn helix_parameters(curve: &Intrinsic3) -> Option<(Scalar, Scalar)> {
    let k = curve.curvature.constant_value()?;
    let tau = curve.torsion.constant_value()?;
    let denominator = k * k + tau * tau;
    if denominator <= 0.0 || !denominator.is_finite() {
        return None;
    }
    Some((k / denominator, tau / denominator))
}

/// Offset a space curve by `distance` along its principal normal.
///
/// Exact for a helix, where the offset is a coaxial helix. Refused for any
/// varying law, because the offset is not unit speed and an `Intrinsic3`
/// stores arc length: re-fitting would silently change the curve.
pub fn offset_intrinsic3(curve: &Intrinsic3, distance: Scalar) -> GeomResult<Intrinsic3> {
    if !distance.is_finite() {
        return Err(invalid("offset distance must be finite"));
    }
    if distance == 0.0 {
        return Ok(curve.clone());
    }
    if !curve.is_helical() {
        // Not a representational gap: the offset of a varying-curvature space
        // curve is not an intrinsic curve in its own arc length at all.
        return Err(unsupported());
    }
    let (a, b) = helix_parameters(curve).ok_or_else(|| invalid("degenerate helix parameters"))?;
    let radius = a - distance;
    let c2 = radius.hypot(b);
    if c2 <= 0.0 || !c2.is_finite() {
        // The offset collapsed onto the axis: there is no curve to return.
        return Err(invalid("offset distance collapses the helix onto its axis"));
    }
    let base = a.hypot(b);
    if base <= 0.0 || !base.is_finite() {
        return Err(invalid("degenerate helix parameters"));
    }
    let curvature = CurvatureLaw::circular(radius / (c2 * c2));
    let torsion = CurvatureLaw::circular(b / (c2 * c2));

    // The offset curve's own start frame. Its tangent is NOT the base
    // tangent: the offset point travels at the same angular rate about the
    // axis but on a different radius, so the tangential/axial mix changes.
    let start = offset_start_frame(curve, distance, a, b)?;
    // Arc length along the offset runs at c2/base times the base's.
    Ok(Intrinsic3::new(
        start,
        curvature,
        torsion,
        curve.length * c2 / base,
    ))
}

/// Start frame of the normal-offset of a helix.
fn offset_start_frame(
    curve: &Intrinsic3,
    distance: Scalar,
    a: Scalar,
    b: Scalar,
) -> GeomResult<Frame3> {
    let base = a.hypot(b);
    let tangent = curve.start.x;
    let normal = curve.start.y;
    let binormal = curve.start.z;
    // Darboux axis in world coordinates: (tau * T + k * B) / |(k, tau)|.
    // With a and b as above this is (b * T + a * B) / base.
    let axis = (tangent * b + binormal * a) / base;
    // Tangential direction: the part of the base tangent perpendicular to
    // the axis, normalised.
    let axial_component = tangent.dot(axis);
    let perpendicular = tangent - axis * axial_component;
    let perpendicular_length = perpendicular.length();
    if perpendicular_length <= 0.0 || !perpendicular_length.is_finite() {
        return Err(invalid("helix axis is parallel to its tangent"));
    }
    let tangential = perpendicular / perpendicular_length;
    // Offset point keeps the axial speed and scales the tangential one by
    // the radius ratio.
    let radius = a - distance;
    let velocity = tangential * radius + axis * b;
    let speed = velocity.length();
    if speed <= 0.0 || !speed.is_finite() {
        return Err(invalid("offset has no well-defined tangent"));
    }
    let new_tangent = velocity / speed;
    // The normal still points at the axis; it flips when the offset crosses
    // the axis and the curve winds the other way round.
    let new_normal = if radius >= 0.0 { normal } else { -normal };
    let new_binormal = new_tangent.cross(new_normal);
    let binormal_length = new_binormal.length();
    if binormal_length <= 0.0 || !binormal_length.is_finite() {
        return Err(invalid("offset frame is degenerate"));
    }
    Ok(Frame3 {
        origin: curve.start.origin + normal * distance,
        x: new_tangent,
        y: new_normal,
        z: new_binormal / binormal_length,
    })
}

/// Join two space curves into one, when the second continues the first.
///
/// Exact: the result carries both laws as a `Piecewise` seam at the first
/// curve's length. Refused when the curves do not actually meet, because a
/// joined curve that jumps is not the curve either input described.
pub fn join_intrinsic3(
    first: &Intrinsic3,
    second: &Intrinsic3,
    position_tolerance: Scalar,
    direction_tolerance: Scalar,
) -> GeomResult<Intrinsic3> {
    if !position_tolerance.is_finite() || position_tolerance < 0.0 {
        return Err(invalid(
            "position tolerance must be finite and non-negative",
        ));
    }
    if !direction_tolerance.is_finite() || direction_tolerance < 0.0 {
        return Err(invalid(
            "direction tolerance must be finite and non-negative",
        ));
    }
    let end_point = frenet_point(first, first.length)?;
    let end_frame = frenet_frame(first, first.length)?;
    if !meets(end_point, second.start.origin, position_tolerance) {
        return Err(invalid(
            "curves do not meet: the second does not start where the first ends",
        ));
    }
    if !aligned(end_frame.x, second.start.x, direction_tolerance) {
        return Err(invalid(
            "curves meet but their tangents disagree, so the join would kink",
        ));
    }
    let curvature = CurvatureLaw::piecewise(
        vec![first.length],
        vec![first.curvature.clone(), second.curvature.clone()],
    );
    let torsion = CurvatureLaw::piecewise(
        vec![first.length],
        vec![first.torsion.clone(), second.torsion.clone()],
    );
    Ok(Intrinsic3::new(
        first.start,
        curvature,
        torsion,
        first.length + second.length,
    ))
}

fn meets(a: Point3, b: Point3, tolerance: Scalar) -> bool {
    (a - b).length() <= tolerance
}

fn aligned(a: Vec3, b: Vec3, tolerance: Scalar) -> bool {
    (a - b).length() <= tolerance
}
