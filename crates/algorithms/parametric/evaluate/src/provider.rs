//! The reference curve-evaluation provider.
//!
//! Implements [`CurveEvaluator`] over the families where a DISTANCE can
//! be honoured exactly, and refuses by name where it cannot. See
//! `docs/adr/0063-curve-evaluation-contract.md`.

use axiolid_contracts::{
    Backend, BackendDescriptor, BackendId, Determinism, ExecutionTarget, GeomError, GeomResult,
};
use axiolid_core::{Frame3, Point3, Scalar, Vec3};
use axiolid_curve::Curve3;
use axiolid_curve_evaluate_contract::{CurveEvaluator, CurveMeasure, DistanceConvention};

use crate::arc_length::{elevated_point, elevated_tangent};
use crate::frenet::{frenet_point, frenet_tangent};
use crate::polyline_length::polyline_parameter;

/// Reference curve evaluator.
///
/// `up` is the reference direction the oriented frame is built against.
/// It is carried by the provider rather than passed per call so that one
/// consumer cannot silently place two objects against different
/// conventions on the same alignment.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct ReferenceCurveEvaluator {
    up: Vec3,
}

impl ReferenceCurveEvaluator {
    /// Backend identity.
    pub const ID: BackendId = BackendId::new("axiolid-evaluate");

    /// Evaluator whose reference up is global `+Z`.
    ///
    /// The right default for civil and building work, where `z` is up by
    /// construction.
    #[must_use]
    pub const fn new() -> Self {
        Self { up: Vec3::Z }
    }

    /// Evaluator with an explicit reference up direction.
    ///
    /// Returns `None` for a non-finite or zero vector, which cannot
    /// define a roll convention.
    #[must_use]
    pub fn with_up(up: Vec3) -> Option<Self> {
        let length = up.length();
        if !length.is_finite() || length <= 0.0 {
            return None;
        }
        Some(Self { up: up / length })
    }

    /// The reference up direction this evaluator builds frames against.
    #[must_use]
    pub const fn up(&self) -> Vec3 {
        self.up
    }
}

impl Default for ReferenceCurveEvaluator {
    fn default() -> Self {
        Self::new()
    }
}

impl Backend for ReferenceCurveEvaluator {
    fn descriptor(&self) -> BackendDescriptor {
        BackendDescriptor::new(Self::ID, ExecutionTarget::PortableCpu)
    }
}

fn unsupported() -> GeomError {
    GeomError::Unsupported {
        backend: ReferenceCurveEvaluator::ID,
        operation: axiolid_contracts::Operation::CurveEvaluation,
    }
}

fn invalid(detail: &str) -> GeomError {
    GeomError::InvalidInput(detail.into())
}

/// The carried number, rejected here if it is not finite.
///
/// A NaN parameter would otherwise reach the evaluators and be refused
/// as a `curve parameter`, naming an internal concept rather than what
/// the caller actually passed.
fn finite_value(at: CurveMeasure) -> GeomResult<Scalar> {
    let value = at.value();
    if !value.is_finite() {
        return Err(invalid("curve measure must be finite"));
    }
    Ok(value)
}
/// Which distance this provider can honour for `curve`.
///
/// The families divide by whether distance is RECOVERABLE in closed
/// form, not by whether they can be evaluated at all:
///
/// - `Intrinsic` is already parameterised by arc length, so distance is
///   the native parameter -- nothing to convert.
/// - `Line` and `Circle` have an elementary arc length, so a distance
///   maps to a parameter exactly (`d / |direction|`, `d / radius`).
/// - `Elevated` is authored against PLAN distance, and that is what it
///   is reported as. It is deliberately NOT called arc length: on a 5%
///   grade the true 3D length exceeds the plan distance by 0.125 m per
///   100 m, and quietly conflating the two would misplace an object by
///   that much.
/// - `Polyline` has an exact arc length: a finite sum of segment
///   lengths, located by a running sum and one linear interpolation. Its
///   only transcendental is the same per-segment `sqrt` that `Line`
///   already reports as `ArcLength3d`, so refusing the sequence while
///   accepting each element was not defensible (kernel#107). Seam,
///   degenerate-segment and closed-wrap behaviour is pinned in
///   `polyline_length`.
/// - `Ellipse` and `BSpline` have no closed-form arc length (the ellipse
///   needs an elliptic integral, the spline needs quadrature plus
///   numeric inversion), so distance is refused rather than approximated
///   behind an exact-looking signature.
fn convention_for(curve: &Curve3) -> DistanceConvention {
    match curve {
        Curve3::Intrinsic(_) | Curve3::Line(_) | Curve3::Circle(_) | Curve3::Polyline(_) => {
            DistanceConvention::ArcLength3d
        }
        Curve3::Elevated(_) => DistanceConvention::PlanDistance,
        _ => DistanceConvention::Unsupported,
    }
}

/// Convert a distance to the curve's native parameter, exactly.
fn parameter_for(curve: &Curve3, distance: Scalar) -> GeomResult<Scalar> {
    if !distance.is_finite() {
        return Err(invalid("distance along a curve must be finite"));
    }
    match curve {
        // Already arc length.
        Curve3::Intrinsic(_) | Curve3::Elevated(_) => Ok(distance),
        // The parameter advances |direction| per unit, so a caller-facing
        // distance must be divided by it. Import adapters may preserve a
        // non-unit direction, so this is not a no-op in practice.
        Curve3::Line(l) => {
            let speed = l.direction.length();
            if !speed.is_finite() || speed <= 0.0 {
                return Err(invalid("line has no direction, so no distance along it"));
            }
            Ok(distance / speed)
        }
        // Angle = arc / radius.
        Curve3::Circle(c) => {
            if !c.radius.is_finite() || c.radius <= 0.0 {
                return Err(invalid("circle has no positive radius"));
            }
            Ok(distance / c.radius)
        }
        // Running sum over segment lengths; exact, no iteration.
        Curve3::Polyline(p) => polyline_parameter(p, distance),
        _ => Err(unsupported()),
    }
}

impl CurveEvaluator for ReferenceCurveEvaluator {
    fn distance_convention(&self, curve: &Curve3) -> DistanceConvention {
        convention_for(curve)
    }

    /// Bitwise: every path is deterministic floating-point arithmetic
    /// with no hashing, threading or iteration-order dependence.
    fn determinism(&self) -> Determinism {
        Determinism::Bitwise
    }

    fn point_at(&self, curve: &Curve3, at: CurveMeasure) -> GeomResult<Point3> {
        // A native parameter needs no conversion and no convention, so it
        // works for EVERY family -- including the ones whose arc length has
        // no closed form and whose distance route is refused.
        let CurveMeasure::Distance(distance) = at else {
            return crate::curve::evaluate3(curve, finite_value(at)?);
        };
        match curve {
            // Native arc-length families: straight through, no conversion.
            Curve3::Intrinsic(i) => frenet_point(i, distance),
            Curve3::Elevated(e) => elevated_point(e, distance),
            Curve3::Line(_) | Curve3::Circle(_) | Curve3::Polyline(_) => {
                crate::curve::evaluate3(curve, parameter_for(curve, distance)?)
            }
            _ => Err(unsupported()),
        }
    }

    fn tangent_at(&self, curve: &Curve3, at: CurveMeasure) -> GeomResult<Vec3> {
        let raw = match at {
            CurveMeasure::Parameter(_) => crate::curve::derivative3(curve, finite_value(at)?)?,
            CurveMeasure::Distance(distance) => match curve {
                Curve3::Intrinsic(i) => frenet_tangent(i, distance)?,
                Curve3::Elevated(e) => elevated_tangent(e, distance)?,
                Curve3::Line(_) | Curve3::Circle(_) | Curve3::Polyline(_) => {
                    crate::curve::derivative3(curve, parameter_for(curve, distance)?)?
                }
                _ => return Err(unsupported()),
            },
            // `CurveMeasure` is #[non_exhaustive]. An unknown method of
            // measurement is refused by name rather than guessed at.
            _ => return Err(unsupported()),
        };
        // The contract promises a UNIT tangent. The intrinsic and elevated
        // evaluators already return one; a line's derivative is its raw
        // direction and a circle's scales with radius, so normalising here
        // is what makes the families interchangeable to a caller.
        let length = raw.length();
        if !length.is_finite() || length <= 0.0 {
            return Err(invalid("curve has no tangent direction there"));
        }
        Ok(raw / length)
    }

    fn frame_at(&self, curve: &Curve3, at: CurveMeasure) -> GeomResult<Frame3> {
        let origin = self.point_at(curve, at)?;
        let tangent = self.tangent_at(curve, at)?;
        // Reference-up construction. `right` is perpendicular to both the
        // tangent and the reference direction; `up` is then recovered from
        // those two so the triad is exactly orthonormal even when the
        // tangent is not perpendicular to the reference direction.
        let right = tangent.cross(self.up);
        let magnitude = right.length();
        // Parallel tangent and reference: roll is genuinely undetermined.
        // Refuse rather than invent one, because a silently arbitrary roll
        // rotates whatever is placed by it.
        if !magnitude.is_finite() || magnitude <= 1e-12 {
            return Err(invalid(
                "tangent is parallel to the reference up direction, so roll is undefined",
            ));
        }
        let right = right / magnitude;
        let up = right.cross(tangent);
        Ok(Frame3 {
            origin,
            x: tangent,
            y: up,
            z: right,
        })
    }
}
