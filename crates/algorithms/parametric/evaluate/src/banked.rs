//! Evaluation of banked curves: point, tangent and the rolled section frame.
//!
//! A [`Banked3`] is evaluated by PLAN distance `d`, the parameter of its
//! centreline. The point is `base(d) + e(d) z`, `e` the pivot law. The
//! tangent is that point's own: its derivative in `d` is
//! `(p, g + e')`, `p` the unit plan tangent, `g` the profile grade and `e'`
//! the pivot's rate, so a pivot that moves (rotation about the low rail
//! through a transition) tilts the tangent by exactly what it adds to the
//! grade. With a constant pivot it is the centreline's tangent.
//!
//! The section is rolled about that tangent `t`. With `n` the horizontal
//! left normal and `u = t x n`, the roll `rho` gives the lateral axis
//! `l = cos(rho) n + sin(rho) u` (towards the left rail head) and the
//! section up `v = -sin(rho) n + cos(rho) u`; `rho` comes from
//! [`BankConvention`](axiolid_curve::BankConvention): `asin(D / b)` for a
//! rotation about the tangent, `asin(D / (b cos theta))` for an exact
//! vertical rise `D`. The rail heads are at `point +- (b / 2) l`. Every
//! refusal of the curve crate's [`BankError`] reaches the caller by name.
//!
//! Position is as exact as the centreline's; the frame adds closed-form
//! trigonometry only.

use axiolid_contracts::{GeomError, GeomResult};
use axiolid_core::{Frame3, Point3, Scalar, Vec3};
use axiolid_curve::{BankError, Banked3, CantValue};

use crate::arc_length::{elevated_point, elevated_tangent};

fn refused(error: BankError) -> GeomError {
    GeomError::InvalidInput(format!("banked curve: {error}"))
}

fn invalid(detail: &str) -> GeomError {
    GeomError::InvalidInput(detail.to_owned())
}

/// The section of a banked curve at one station.
#[non_exhaustive]
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct BankedSection {
    /// The rotation point: the centreline point raised by the pivot.
    pub point: Point3,
    /// Unit tangent of the banked curve.
    pub tangent: Vec3,
    /// Unit lateral axis, square to the tangent, towards the left rail head.
    pub lateral: Vec3,
    /// Unit section up, `tangent x lateral`.
    pub up: Vec3,
    /// Cant `D`: the left rail head above the right, as the law states it.
    pub cant: Scalar,
    /// Nominal bank angle `psi = asin(D / b)`, radians.
    pub bank_angle: Scalar,
    /// The roll about the tangent the convention gives, radians: `psi`
    /// itself for a rotation about the tangent.
    pub roll: Scalar,
    /// Elevation of the rotation point above the profile.
    pub pivot: Scalar,
    /// Grade of the banked curve, `dz/dd`: the profile's plus the pivot's
    /// rate.
    pub grade: Scalar,
    /// Rail-head distance `b`.
    pub rail_head_distance: Scalar,
}

impl BankedSection {
    /// The section as a frame, in the curve-evaluation contract's layout:
    /// `x` the tangent, `y` the section up, `z = x x y` (to the right).
    #[must_use]
    pub fn frame(&self) -> Frame3 {
        Frame3 {
            origin: self.point,
            x: self.tangent,
            y: self.up,
            z: -self.lateral,
        }
    }

    /// Left and right rail heads, `b / 2` either side of the rotation point
    /// along the lateral axis.
    #[must_use]
    pub fn rail_heads(&self) -> (Point3, Point3) {
        let half = 0.5 * self.rail_head_distance * self.lateral;
        (self.point + half, self.point - half)
    }
}

/// Unit plan tangent, horizontal, at plan distance `d`.
fn plan_direction(curve: &Banked3, d: Scalar) -> GeomResult<Vec3> {
    let tangent = elevated_tangent(&curve.base, d)?;
    let horizontal = Vec3::new(tangent.x, tangent.y, 0.0);
    let length = horizontal.length();
    if !length.is_finite() || length <= 0.0 {
        return Err(invalid("banked curve: centreline has no plan direction"));
    }
    Ok(horizontal / length)
}

/// The cant law covers `d`, or the station has no section.
fn covered(curve: &Banked3, d: Scalar) -> GeomResult<()> {
    curve.cant.value_at(d).map(|_| ()).ok_or_else(|| {
        refused(BankError::OutsideLaw {
            law: "cant",
            distance: d,
        })
    })
}

/// Position on a banked curve at plan distance `d`: `base(d) + e(d) z`.
///
/// # Errors
///
/// The centreline's refusals, and the pivot and cant laws' by name
/// outside their span.
pub fn banked_point(curve: &Banked3, d: Scalar) -> GeomResult<Point3> {
    if !d.is_finite() {
        return Err(invalid("plan distance must be finite"));
    }
    covered(curve, d)?;
    let (pivot, _) = curve.pivot_at(d).map_err(refused)?;
    let base = elevated_point(&curve.base, d)?;
    Ok(base + pivot * Vec3::Z)
}

/// Derivative of [`banked_point`] in plan distance: `(p, g + e')`. Not
/// unit: its length is `sqrt(1 + (g + e')^2)`.
///
/// # Errors
///
/// As [`banked_point`].
pub fn banked_derivative(curve: &Banked3, d: Scalar) -> GeomResult<Vec3> {
    if !d.is_finite() {
        return Err(invalid("plan distance must be finite"));
    }
    covered(curve, d)?;
    let (_, pivot_rate) = curve.pivot_at(d).map_err(refused)?;
    let plan = plan_direction(curve, d)?;
    let grade = curve
        .base
        .elevation
        .grade_at(d)
        .ok_or_else(|| invalid("elevation law has no grade at that distance"))?;
    let derivative = Vec3::new(plan.x, plan.y, grade + pivot_rate);
    if !derivative.is_finite() {
        return Err(invalid("banked curve: grade is not finite"));
    }
    Ok(derivative)
}

/// Unit tangent of a banked curve at plan distance `d`.
///
/// # Errors
///
/// As [`banked_point`].
pub fn banked_tangent(curve: &Banked3, d: Scalar) -> GeomResult<Vec3> {
    Ok(banked_derivative(curve, d)?.normalize())
}

/// The section at plan distance `d`: rotation point, rolled frame, cant,
/// bank angle and roll.
///
/// # Errors
///
/// As [`banked_point`], and by name: a cant beyond the rail-head distance,
/// beyond the vertical span `b cos(theta)` under a vertical-rise
/// convention, an angle piece beyond a quarter turn, an angle piece in the
/// pivot law, a bad rail-head distance.
pub fn banked_section(curve: &Banked3, d: Scalar) -> GeomResult<BankedSection> {
    let point = banked_point(curve, d)?;
    let derivative = banked_derivative(curve, d)?;
    let (pivot, _) = curve.pivot_at(d).map_err(refused)?;
    let b = curve.rail_head_distance;
    let tangent = derivative.normalize();
    // `cos(theta)` from the derivative `(p, G)`: `1 / sqrt(1 + G^2)`.
    let grade = derivative.z;
    let grade_cosine = 1.0 / grade.hypot(1.0);
    let cant = curve.cant_at(d).map_err(refused)?;
    let bank_angle = curve.bank_angle_at(d).map_err(refused)?;
    let roll = match curve.cant.value_at(d) {
        Some(CantValue::Angle(psi)) => curve.convention.roll_from_angle(psi, b, grade_cosine),
        _ => curve.convention.roll(cant, b, grade_cosine),
    }
    .map_err(refused)?;
    let plan = Vec3::new(derivative.x, derivative.y, 0.0);
    let normal = Vec3::Z.cross(plan).normalize();
    let square_up = tangent.cross(normal);
    let (sin, cos) = roll.sin_cos();
    let lateral = cos * normal + sin * square_up;
    let up = -sin * normal + cos * square_up;
    Ok(BankedSection {
        point,
        tangent,
        lateral,
        up,
        cant,
        bank_angle,
        roll,
        pivot,
        grade,
        rail_head_distance: b,
    })
}
