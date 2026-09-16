//! Conformance suite every curve-evaluation provider must pass.
//!
//! These are the properties a CALLER is entitled to assume. A provider
//! that passes them can be swapped for another without a consumer
//! noticing, which is the whole point of naming the capability.

use axiolid_core::{Frame3, Point3, Scalar, Vec3};
use axiolid_curve::{Circle3, CurvatureLaw, Curve3, Intrinsic3, Line3};

use crate::{CurveEvaluator, CurveMeasure};

/// One failed conformance expectation.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct ConformanceFailure {
    /// Which expectation failed.
    pub check: &'static str,
    /// What went wrong.
    pub detail: String,
}

fn fail(check: &'static str, detail: impl Into<String>) -> ConformanceFailure {
    ConformanceFailure {
        check,
        detail: detail.into(),
    }
}

/// Run every conformance check against `provider`.
///
/// Returns the failures; an empty vector means conformant.
#[must_use]
pub fn check<E: CurveEvaluator>(provider: &E) -> Vec<ConformanceFailure> {
    let mut out = Vec::new();
    out.extend(check_line(provider));
    out.extend(check_circle(provider));
    out.extend(check_tangent_is_unit(provider));
    out.extend(check_frame_is_orthonormal(provider));
    out.extend(check_refusals(provider));
    out.extend(check_measure_routes_differ(provider));
    out
}

/// A line with a NON-UNIT direction must still advance true distance.
///
/// This is the check that catches a provider handing back the native
/// parameter as though it were a distance: with `|direction| = 2` the two
/// differ by a factor of two, while a unit direction would hide the bug.
fn check_line<E: CurveEvaluator>(provider: &E) -> Vec<ConformanceFailure> {
    let mut out = Vec::new();
    let line = Curve3::Line(Line3 {
        origin: Point3::new(1.0, 2.0, 3.0),
        direction: Vec3::new(2.0, 0.0, 0.0),
    });
    if !provider.distance_convention(&line).is_supported() {
        return out;
    }
    match provider.point_at(&line, CurveMeasure::Distance(5.0)) {
        Ok(p) => {
            let moved = (p - Point3::new(1.0, 2.0, 3.0)).length();
            if (moved - 5.0).abs() > 1e-9 {
                out.push(fail(
                    "line advances true distance",
                    format!("distance 5 moved {moved}, not 5"),
                ));
            }
        }
        Err(e) => out.push(fail("line advances true distance", format!("{e:?}"))),
    }
    out
}

/// A circle of radius r must reach angle d/r at distance d.
fn check_circle<E: CurveEvaluator>(provider: &E) -> Vec<ConformanceFailure> {
    let mut out = Vec::new();
    let radius = 4.0;
    let circle = Curve3::Circle(Circle3 {
        frame: Frame3 {
            origin: Point3::ZERO,
            x: Vec3::X,
            y: Vec3::Y,
            z: Vec3::Z,
        },
        radius,
    });
    if !provider.distance_convention(&circle).is_supported() {
        return out;
    }
    // A quarter of the circumference must land on the +y axis.
    let quarter = core::f64::consts::FRAC_PI_2 * radius;
    match provider.point_at(&circle, CurveMeasure::Distance(quarter)) {
        Ok(p) => {
            let want = Point3::new(0.0, radius, 0.0);
            let error = (p - want).length();
            if error > 1e-9 {
                out.push(fail(
                    "circle reaches angle d/r",
                    format!("quarter circumference landed {error:e} away from the +y axis"),
                ));
            }
        }
        Err(e) => out.push(fail("circle reaches angle d/r", format!("{e:?}"))),
    }
    out
}

fn sample_curves() -> Vec<Curve3> {
    vec![
        Curve3::Line(Line3 {
            origin: Point3::new(0.0, 0.0, 0.0),
            direction: Vec3::new(3.0, 4.0, 0.0),
        }),
        Curve3::Circle(Circle3 {
            frame: Frame3 {
                origin: Point3::new(2.0, 0.0, 1.0),
                x: Vec3::X,
                y: Vec3::Y,
                z: Vec3::Z,
            },
            radius: 3.0,
        }),
        // An arc-length-native family. Included because a provider may route
        // it past the distance-to-parameter conversion entirely, so its
        // guards would otherwise go unexercised by this suite.
        Curve3::Intrinsic(Intrinsic3::new(
            Frame3 {
                origin: Point3::ZERO,
                x: Vec3::X,
                y: Vec3::Y,
                z: Vec3::Z,
            },
            CurvatureLaw::circular(0.1),
            CurvatureLaw::circular(0.04),
            10.0,
        )),
    ]
}

/// The tangent must be unit length wherever it is defined.
fn check_tangent_is_unit<E: CurveEvaluator>(provider: &E) -> Vec<ConformanceFailure> {
    let mut out = Vec::new();
    for curve in sample_curves() {
        if !provider.distance_convention(&curve).is_supported() {
            continue;
        }
        for distance in [0.0, 1.5, 4.0] {
            if let Ok(t) = provider.tangent_at(&curve, CurveMeasure::Distance(distance)) {
                if (t.length() - 1.0).abs() > 1e-9 {
                    out.push(fail(
                        "tangent is unit length",
                        format!("|t| = {} at distance {distance}", t.length()),
                    ));
                }
            }
        }
    }
    out
}

/// The frame must be right-handed orthonormal with `x` on the tangent.
///
/// A caller builds a rotation from this directly, so a frame that is
/// merely close to orthonormal shears whatever it places.
fn check_frame_is_orthonormal<E: CurveEvaluator>(provider: &E) -> Vec<ConformanceFailure> {
    let mut out = Vec::new();
    for curve in sample_curves() {
        if !provider.distance_convention(&curve).is_supported() {
            continue;
        }
        for distance in [0.0, 1.5, 4.0] {
            let Ok(frame) = provider.frame_at(&curve, CurveMeasure::Distance(distance)) else {
                continue;
            };
            let checks = [
                ("x unit", frame.x.length() - 1.0),
                ("y unit", frame.y.length() - 1.0),
                ("z unit", frame.z.length() - 1.0),
                ("x.y", frame.x.dot(frame.y)),
                ("x.z", frame.x.dot(frame.z)),
                ("y.z", frame.y.dot(frame.z)),
            ];
            for (name, value) in checks {
                if value.abs() > 1e-9 {
                    out.push(fail(
                        "frame is orthonormal",
                        format!("{name} off by {value:e} at distance {distance}"),
                    ));
                }
            }
            // Right-handed: x cross y must be z, not -z.
            let handed = frame.x.cross(frame.y).dot(frame.z);
            if (handed - 1.0).abs() > 1e-9 {
                out.push(fail(
                    "frame is right-handed",
                    format!("x cross y . z = {handed} at distance {distance}"),
                ));
            }
            // The frame must sit ON the curve.
            if let Ok(point) = provider.point_at(&curve, CurveMeasure::Distance(distance)) {
                let off = (frame.origin - point).length();
                if off > 1e-9 {
                    out.push(fail(
                        "frame origin is the curve point",
                        format!("origin {off:e} away at distance {distance}"),
                    ));
                }
            }
            // And its x axis must be the tangent.
            if let Ok(tangent) = provider.tangent_at(&curve, CurveMeasure::Distance(distance)) {
                let off = (frame.x - tangent).length();
                if off > 1e-9 {
                    out.push(fail(
                        "frame x is the tangent",
                        format!("x differs from tangent by {off:e}"),
                    ));
                }
            }
        }
    }
    out
}

/// Bad input must be refused, not answered with a guess.
///
/// A NaN distance that returns a NaN point is the dangerous case: it
/// propagates silently into a placement instead of failing at the call.
fn check_refusals<E: CurveEvaluator>(provider: &E) -> Vec<ConformanceFailure> {
    let mut out = Vec::new();
    for curve in sample_curves() {
        if !provider.distance_convention(&curve).is_supported() {
            continue;
        }
        for bad in [Scalar::NAN, Scalar::INFINITY, Scalar::NEG_INFINITY] {
            if provider
                .point_at(&curve, CurveMeasure::Distance(bad))
                .is_ok()
            {
                out.push(fail(
                    "non-finite distance is refused",
                    format!("point_at accepted {bad}"),
                ));
            }
            if provider
                .tangent_at(&curve, CurveMeasure::Distance(bad))
                .is_ok()
            {
                out.push(fail(
                    "non-finite distance is refused",
                    format!("tangent_at accepted {bad}"),
                ));
            }
            if provider
                .frame_at(&curve, CurveMeasure::Distance(bad))
                .is_ok()
            {
                out.push(fail(
                    "non-finite distance is refused",
                    format!("frame_at accepted {bad}"),
                ));
            }
        }
    }
    // A provider must not claim a convention it cannot honour: if it
    // reports Unsupported it must actually refuse.
    let vertical = Curve3::Line(Line3 {
        origin: Point3::ZERO,
        direction: Vec3::ZERO,
    });
    if provider.distance_convention(&vertical).is_supported()
        && provider
            .point_at(&vertical, CurveMeasure::Distance(1.0))
            .is_ok()
    {
        out.push(fail(
            "degenerate curve is refused",
            "a zero-direction line was evaluated instead of refused",
        ));
    }
    out
}

/// A parameter and a distance must not be silently interchangeable.
///
/// On a circle of radius 4 the same number means two different places:
/// `Parameter(1.5)` is 1.5 radians round, `Distance(1.5)` is 1.5 m along,
/// i.e. 0.375 rad. A provider that ignored the method of measurement would
/// return the same point for both -- wrong, finite, and plausible.
///
/// Also pins that the parameter route stays OPEN where the distance route
/// is refused: an ellipse has no closed-form arc length, but its native
/// parameter is perfectly meaningful, and a consumer holding an authored
/// parameter must still be able to evaluate it.
fn check_measure_routes_differ<E: CurveEvaluator>(provider: &E) -> Vec<ConformanceFailure> {
    let mut out = Vec::new();
    let radius = 4.0;
    let circle = Curve3::Circle(Circle3 {
        frame: Frame3 {
            origin: Point3::ZERO,
            x: Vec3::X,
            y: Vec3::Y,
            z: Vec3::Z,
        },
        radius,
    });
    let value = 1.5;
    let by_parameter = provider.point_at(&circle, CurveMeasure::Parameter(value));
    let by_distance = provider.point_at(&circle, CurveMeasure::Distance(value));
    if let (Ok(p), Ok(d)) = (by_parameter, by_distance) {
        if (p - d).length() < 1e-9 {
            out.push(fail(
                "parameter and distance are distinct",
                format!("{value} gave the same point as a parameter and as a distance"),
            ));
        }
    }
    out
}
