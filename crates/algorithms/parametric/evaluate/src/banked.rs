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
use axiolid_curve::{BankError, Banked3, CantForm, CantLaw, CantPiece, CantValue};

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
    // Every elevation law, including an intrinsic profile (#238).
    let grade = crate::elevation::elevation_grade(&curve.base.elevation, d)?;
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

/// Second derivative of [`banked_point`] in plan distance (#252):
/// the centreline's `(p'', z'')` plus the pivot's `e''` vertically.
///
/// # Errors
///
/// As [`banked_point`].
pub fn banked_second_derivative(curve: &Banked3, d: Scalar) -> GeomResult<Vec3> {
    if !d.is_finite() {
        return Err(invalid("plan distance must be finite"));
    }
    covered(curve, d)?;
    if curve.pivot.has_angle_pieces() {
        return Err(refused(BankError::AngleInPivot));
    }
    let [_, _, bend] = pivot_jet(&curve.pivot, d, false).ok_or_else(|| {
        refused(BankError::OutsideLaw {
            law: "pivot",
            distance: d,
        })
    })?;
    let base = crate::elevated::elevated_second_derivative(&curve.base, d)?;
    Ok(base + bend * Vec3::Z)
}

/// The height piece of a law at `d` and the distance into it: the piece
/// starting there, or with `before` the one ending there.
fn piece_at(law: &CantLaw, d: Scalar, before: bool) -> Option<(&CantPiece, Scalar)> {
    if !d.is_finite() || d < 0.0 || !law.is_well_formed() {
        return None;
    }
    let mut start = 0.0;
    let last = law.pieces.len() - 1;
    for (index, piece) in law.pieces.iter().enumerate() {
        let end = start + piece.length;
        let inside = if before { d <= end } else { d < end };
        if inside || (index == last && d <= end) {
            return Some((piece, d - start));
        }
        start = end;
    }
    None
}

/// `(e, e', e'')` of a height law at `d`, in plan distance; `None` off
/// the law or for an angle piece.
pub(crate) fn pivot_jet(law: &CantLaw, d: Scalar, before: bool) -> Option<[Scalar; 3]> {
    let (piece, s) = piece_at(law, d, before)?;
    let length = piece.length;
    let xi = s / length;
    let (value, slope, bend) = match &piece.form {
        CantForm::Polynomial { coefficients } => {
            let horner =
                |k: usize| -> Scalar {
                    coefficients.iter().enumerate().skip(k).rev().fold(
                        0.0,
                        |accumulated, (i, c)| {
                            let falling: Scalar = (0..k).map(|j| (i - j) as Scalar).product();
                            accumulated * xi + c * falling
                        },
                    )
                };
            (horner(0), horner(1), horner(2))
        }
        CantForm::Cosine { start, change } => {
            let pi = core::f64::consts::PI;
            let (sin, cos) = (pi * xi).sin_cos();
            (
                start + change * 0.5 * (1.0 - cos),
                change * 0.5 * pi * sin,
                change * 0.5 * pi * pi * cos,
            )
        }
        CantForm::Sine { start, change } => {
            let tau = core::f64::consts::TAU;
            let (sin, cos) = (tau * xi).sin_cos();
            (
                start + change * (xi - sin / tau),
                change * (1.0 - cos),
                change * tau * sin,
            )
        }
        _ => return None,
    };
    let out = [value, slope / length, bend / (length * length)];
    out.iter().all(|v| v.is_finite()).then_some(out)
}

/// Certified `(sup |e'|, sup |e''|, sup |e'''|)` of a height law over
/// `[lo, hi]` inside ONE piece (#252); `None` across a seam, off the law,
/// or for an angle piece.
///
/// In the piece's `xi = s / L`, the `k`-th derivative in plan distance is
/// the form's `k`-th in `xi` over `L^k`: a polynomial's by
/// `sum i!/(i-k)! |c_i| m^(i-k)` with `m` the larger `|xi|`, the
/// half-cosine's `|dD| / 2 pi^k`, the sine transition's `2 |dD|`,
/// `2 pi |dD|`, `4 pi^2 |dD|`.
pub(crate) fn pivot_bounds(law: &CantLaw, lo: Scalar, hi: Scalar) -> Option<[Scalar; 3]> {
    let (lo, hi) = (lo.min(hi), lo.max(hi));
    let (piece, s_lo) = piece_at(law, lo, false)?;
    let s_hi = s_lo + (hi - lo);
    if s_hi > piece.length * (1.0 + 1e-12) {
        return None;
    }
    let length = piece.length;
    let in_xi = match &piece.form {
        CantForm::Polynomial { coefficients } => {
            let m = (s_lo / length).abs().max((s_hi / length).abs());
            let term = |k: usize| -> Scalar {
                coefficients
                    .iter()
                    .enumerate()
                    .skip(k)
                    .map(|(i, c)| {
                        let falling: Scalar = (0..k).map(|j| (i - j) as Scalar).product();
                        falling * c.abs() * m.powi((i - k) as i32)
                    })
                    .sum()
            };
            [term(1), term(2), term(3)]
        }
        CantForm::Cosine { change, .. } => {
            let pi = core::f64::consts::PI;
            let half = 0.5 * change.abs();
            [half * pi, half * pi * pi, half * pi * pi * pi]
        }
        CantForm::Sine { change, .. } => {
            let tau = core::f64::consts::TAU;
            let c = change.abs();
            [2.0 * c, tau * c, tau * tau * c]
        }
        _ => return None,
    };
    let out = [
        in_xi[0] / length,
        in_xi[1] / (length * length),
        in_xi[2] / (length * length * length),
    ]
    .map(|v| v * (1.0 + crate::bound::ROUNDING));
    out.iter().all(|v| v.is_finite()).then_some(out)
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
