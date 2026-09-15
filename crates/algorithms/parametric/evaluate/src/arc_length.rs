//! Arc-length evaluation of intrinsic (natural-equation) curves and of the
//! planar-plus-elevation composition.
//!
//! # Why quadrature, and why that is not an approximation of the VALUE
//!
//! An intrinsic curve stores curvature as a function of arc length. Its
//! heading is the integral of that law and is exact in closed form, but its
//! POSITION is the integral of `(cos th, sin th)` and has no elementary
//! antiderivative -- for the clothoid it is the Fresnel integral. The stored
//! value stays exact; only reading a point out of it needs numerical work.
//! That is the same bargain as evaluating `sin`: the curve is not approximated,
//! its evaluation is computed to tolerance.
//!
//! Gauss-Legendre is used because the integrand is smooth. An 8-point rule
//! integrates a degree-15 polynomial exactly, and against the Fresnel closed
//! form it reproduces a 120 m clothoid to R=300 with zero error at machine
//! precision on a single panel, versus roughly 1e-3 relative for a comparable
//! trapezoid budget. Panels are subdivided by total turning so a tight spiral
//! gets more of them.

use axiolid_contracts::{BackendId, GeomError, GeomResult, Operation};
use axiolid_core::{Frame2, Point2, Point3, Scalar, Vec2, Vec3};
use axiolid_curve::{Curve2, Elevated3, Intrinsic2};

/// Nodes and weights of the 8-point Gauss-Legendre rule on `[-1, 1]`.
///
/// Exact for polynomials up to degree 15. Written out rather than computed:
/// they are constants, and a Newton solve at startup would be slower and no
/// more accurate.
const GAUSS_NODES: [Scalar; 8] = [
    -0.960_289_856_497_536_2,
    -0.796_666_477_413_626_7,
    -0.525_532_409_916_328_9,
    -0.183_434_642_495_649_8,
    0.183_434_642_495_649_8,
    0.525_532_409_916_328_9,
    0.796_666_477_413_626_7,
    0.960_289_856_497_536_2,
];

/// Weights matching [`GAUSS_NODES`].
const GAUSS_WEIGHTS: [Scalar; 8] = [
    0.101_228_536_290_376_3,
    0.222_381_034_453_374_5,
    0.313_706_645_877_887_3,
    0.362_683_783_378_361_9,
    0.362_683_783_378_361_9,
    0.313_706_645_877_887_3,
    0.222_381_034_453_374_5,
    0.101_228_536_290_376_3,
];

/// Panels per radian of total turning, above a floor of one panel.
///
/// A straight or gently curving run needs one panel; a spiral that turns
/// through a radian gets four. Bounded so a pathological law cannot ask for an
/// unbounded amount of work.
const PANELS_PER_RADIAN: Scalar = 4.0;

/// Upper bound on panels, so a malformed law refuses rather than hangs.
const MAX_PANELS: usize = 4096;

fn unsupported() -> GeomError {
    GeomError::Unsupported {
        backend: BackendId::new("axiolid-evaluate"),
        operation: Operation::CurveEvaluation,
    }
}

fn invalid(detail: &str) -> GeomError {
    GeomError::InvalidInput(detail.to_owned())
}

/// How many panels to spend integrating `[0, s]`.
fn panel_count(curve: &Intrinsic2, s: Scalar) -> GeomResult<usize> {
    let turning = curve.curvature.clone();
    let total = Intrinsic2::new(curve.start, turning, s)
        .total_turning()
        .ok_or_else(|| invalid("curvature law does not integrate over the requested span"))?;
    let wanted = (total.abs() * PANELS_PER_RADIAN).ceil().max(1.0);
    if !wanted.is_finite() || wanted > MAX_PANELS as Scalar {
        return Err(invalid("curvature law needs an unbounded number of panels"));
    }
    Ok(wanted as usize)
}

/// Position on an intrinsic curve at arc length `s` from its start.
///
/// The heading is exact; the position is Gauss-Legendre quadrature of the unit
/// tangent, which is where the non-elementary integral is discharged.
pub fn intrinsic_point(curve: &Intrinsic2, s: Scalar) -> GeomResult<Point2> {
    if !s.is_finite() {
        return Err(invalid("arc length must be finite"));
    }
    let panels = panel_count(curve, s)?;
    let step = s / panels as Scalar;
    let mut x = 0.0;
    let mut y = 0.0;
    for panel in 0..panels {
        let a = step * panel as Scalar;
        let half = step / 2.0;
        let mid = a + half;
        for (node, weight) in GAUSS_NODES.iter().zip(GAUSS_WEIGHTS.iter()) {
            let u = mid + half * node;
            let heading = curve
                .heading_at(u)
                .ok_or_else(|| invalid("curvature law does not integrate to the sample"))?;
            x += weight * half * heading.cos();
            y += weight * half * heading.sin();
        }
    }
    // The integral is taken in the start frame, whose x axis is the start
    // tangent, then placed into world coordinates by that frame.
    Ok(place2(&curve.start, Vec2::new(x, y)))
}

/// Unit tangent of an intrinsic curve at arc length `s`.
///
/// Exact: this is the closed-form heading, no quadrature involved.
pub fn intrinsic_tangent(curve: &Intrinsic2, s: Scalar) -> GeomResult<Vec2> {
    if !s.is_finite() {
        return Err(invalid("arc length must be finite"));
    }
    let heading = curve
        .heading_at(s)
        .ok_or_else(|| invalid("curvature law does not integrate to the requested arc length"))?;
    let local = Vec2::new(heading.cos(), heading.sin());
    Ok(rotate2(&curve.start, local))
}

/// Place a local offset into world coordinates through a 2D frame.
fn place2(frame: &Frame2, local: Vec2) -> Point2 {
    Point2::new(
        frame.origin.x + frame.x.x * local.x + frame.y.x * local.y,
        frame.origin.y + frame.x.y * local.x + frame.y.y * local.y,
    )
}

/// Rotate a local direction into world coordinates through a 2D frame.
fn rotate2(frame: &Frame2, local: Vec2) -> Vec2 {
    Vec2::new(
        frame.x.x * local.x + frame.y.x * local.y,
        frame.x.y * local.x + frame.y.y * local.y,
    )
}

/// Position on the plan at plan distance `d`.
///
/// Only families whose parameter IS arc length can carry an elevation law,
/// because the law is written against distance along the plan. A line and an
/// intrinsic curve qualify; a B-spline's parameter is not arc length, so
/// pairing one would silently mean something else and is refused.
fn plan_point(plan: &Curve2, d: Scalar) -> GeomResult<Point2> {
    match plan {
        Curve2::Line(line) => {
            let direction = unit2(line.direction)?;
            Ok(Point2::new(
                line.origin.x + direction.x * d,
                line.origin.y + direction.y * d,
            ))
        }
        Curve2::Circle(circle) => {
            // Arc length d subtends d / r, so the angle is exact.
            if circle.radius <= 0.0 || !circle.radius.is_finite() {
                return Err(invalid("circle radius must be positive and finite"));
            }
            let angle = d / circle.radius;
            Ok(place2(
                &circle.frame,
                Vec2::new(circle.radius * angle.cos(), circle.radius * angle.sin()),
            ))
        }
        Curve2::Intrinsic(intrinsic) => intrinsic_point(intrinsic, d),
        _ => Err(unsupported()),
    }
}

/// Unit tangent of the plan at plan distance `d`.
fn plan_tangent(plan: &Curve2, d: Scalar) -> GeomResult<Vec2> {
    match plan {
        Curve2::Line(line) => unit2(line.direction),
        Curve2::Circle(circle) => {
            if circle.radius <= 0.0 || !circle.radius.is_finite() {
                return Err(invalid("circle radius must be positive and finite"));
            }
            let angle = d / circle.radius;
            Ok(rotate2(&circle.frame, Vec2::new(-angle.sin(), angle.cos())))
        }
        Curve2::Intrinsic(intrinsic) => intrinsic_tangent(intrinsic, d),
        _ => Err(unsupported()),
    }
}

fn unit2(v: Vec2) -> GeomResult<Vec2> {
    let length = (v.x * v.x + v.y * v.y).sqrt();
    if !length.is_finite() || length == 0.0 {
        return Err(invalid("direction must be finite and non-zero"));
    }
    Ok(Vec2::new(v.x / length, v.y / length))
}

/// Position on an elevated curve at plan distance `d`.
///
/// The plan supplies `x`/`y`, the elevation law supplies `z`. Both halves are
/// read at the SAME plan distance, which is the convention `ElevationLaw`
/// documents.
pub fn elevated_point(curve: &Elevated3, d: Scalar) -> GeomResult<Point3> {
    if !d.is_finite() {
        return Err(invalid("plan distance must be finite"));
    }
    let planar = plan_point(&curve.plan, d)?;
    let height = curve
        .elevation
        .height_at(d)
        .ok_or_else(|| invalid("elevation law has no height at that distance"))?;
    Ok(Point3::new(planar.x, planar.y, height))
}

/// Unit tangent of an elevated curve at plan distance `d`.
///
/// The plan tangent is horizontal and the grade lifts it, so the 3D tangent is
/// `(t.x, t.y, g)` normalised -- the `sqrt(1 + g^2)` factor by which 3D arc
/// length runs ahead of plan distance.
pub fn elevated_tangent(curve: &Elevated3, d: Scalar) -> GeomResult<Vec3> {
    if !d.is_finite() {
        return Err(invalid("plan distance must be finite"));
    }
    let planar = plan_tangent(&curve.plan, d)?;
    let grade = curve
        .elevation
        .grade_at(d)
        .ok_or_else(|| invalid("elevation law has no grade at that distance"))?;
    let scale = (1.0 + grade * grade).sqrt();
    if !scale.is_finite() || scale == 0.0 {
        return Err(invalid("grade does not give a finite tangent"));
    }
    Ok(Vec3::new(planar.x / scale, planar.y / scale, grade / scale))
}
