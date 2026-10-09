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
//! A held-rail pivot piece (`CantForm::AboutRail`, #279) is derived from
//! the cant law: `e = e0 + s D / 2`, `s = +-1`. Over a height piece of the
//! cant its derivatives are `s D^(k) / 2`; over an angle piece, with
//! `D = b sin(psi)`,
//! `e' = s (b / 2) cos(psi) psi'`,
//! `e'' = s (b / 2) (cos(psi) psi'' - sin(psi) psi'^2)` and
//! `e''' = s (b / 2) (cos(psi) psi''' - 3 sin(psi) psi' psi'' - cos(psi) psi'^3)`,
//! so with `P_k = sup |psi^(k)|` the certified bounds are
//! `(b / 2) P_1`, `(b / 2) (P_2 + P_1^2)` and
//! `(b / 2) (P_3 + 3 P_1 P_2 + P_1^3)`. The Viennese bend's `P_k` are
//! exact over the span asked for: each derivative of its shape is taken at
//! the span's ends and at its critical points inside (`viennese_sup`);
//! over a whole piece they are `140 / 64`, `16.8 / sqrt(5)` and `52.5`
//! times `|dpsi| / L^k`. Where the pivot is derived, the cant law's seams
//! are the pivot's too (`derived_cant_seams`).
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
    if curve.cant.has_rail_pieces() {
        return Err(refused(BankError::RailInCant));
    }
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
    // Every refusal of the laws by name first, then the jet.
    curve.pivot_at(d).map_err(refused)?;
    let [_, _, bend] = pivot_jet(curve, d, false).ok_or_else(|| {
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

/// `(e, e', e'')` of a banked curve's pivot at `d`, in plan distance,
/// read from the piece starting there or with `before` the one ending
/// there; `None` off the law, for an angle piece, or where a held-rail
/// piece's cant has no value.
pub(crate) fn pivot_jet(curve: &Banked3, d: Scalar, before: bool) -> Option<[Scalar; 3]> {
    let (piece, s) = piece_at(&curve.pivot, d, before)?;
    if let CantForm::AboutRail { rail, elevation } = piece.form {
        let half = 0.5 * rail.pivot_sign();
        let [cant, rate, bend] = cant_jet(curve, d, before)?;
        let out = [elevation + half * cant, half * rate, half * bend];
        return out.iter().all(|v| v.is_finite()).then_some(out);
    }
    height_jet(piece, s)
}

/// `(D, D', D'')` of a banked curve's cant at `d`, in plan distance: a
/// height piece's own, an angle piece's through `D = b sin(psi)`. `None`
/// off the law, for a held-rail piece, or where `|D| > b` or
/// `|psi| > pi / 2`.
fn cant_jet(curve: &Banked3, d: Scalar, before: bool) -> Option<[Scalar; 3]> {
    let b = curve.rail_head_distance;
    if !(b.is_finite() && b > 0.0) {
        return None;
    }
    let (piece, s) = piece_at(&curve.cant, d, before)?;
    let out = if let CantForm::VienneseBend { start, change } = piece.form {
        let length = piece.length;
        let xi = s / length;
        let xi2 = xi * xi;
        let one_minus = 1.0 - xi;
        let shape = xi2 * xi2 * (35.0 + xi * (-84.0 + xi * (70.0 - 20.0 * xi)));
        let slope = 140.0 * xi2 * xi * one_minus * one_minus * one_minus;
        let bend = 420.0 * xi2 * one_minus * one_minus * (1.0 - 2.0 * xi);
        let psi = start + change * shape;
        if !psi.is_finite() || psi.abs() > core::f64::consts::FRAC_PI_2 {
            return None;
        }
        let rate = change * slope / length;
        let accel = change * bend / (length * length);
        let (sin, cos) = psi.sin_cos();
        [
            b * sin,
            b * cos * rate,
            b * (cos * accel - sin * rate * rate),
        ]
    } else {
        let jet = height_jet(piece, s)?;
        if jet[0].abs() > b {
            return None;
        }
        jet
    };
    out.iter().all(|v| v.is_finite()).then_some(out)
}

/// `(v, v', v'')` of a height piece `s` into it, in plan distance; `None`
/// for an angle or held-rail piece.
fn height_jet(piece: &CantPiece, s: Scalar) -> Option<[Scalar; 3]> {
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

/// Certified `(sup |e'|, sup |e''|, sup |e'''|)` of a banked curve's pivot
/// over `[lo, hi]` inside ONE piece (#252), and for a held-rail piece
/// inside one piece of the cant law too (#279: half the cant's bounds, see
/// the [module documentation](self)); `None` across a seam, off the law,
/// or for an angle piece.
pub(crate) fn pivot_bounds(curve: &Banked3, lo: Scalar, hi: Scalar) -> Option<[Scalar; 3]> {
    let (lo, hi) = (lo.min(hi), lo.max(hi));
    let (piece, s_lo) = within_one_piece(&curve.pivot, lo, hi)?;
    if let CantForm::AboutRail { .. } = piece.form {
        return cant_bounds(curve, lo, hi).map(|bounds| bounds.map(|v| 0.5 * v));
    }
    height_bounds(piece, s_lo, s_lo + (hi - lo))
}

/// The piece of `law` holding all of `[lo, hi]`, and `lo`'s distance into
/// it.
fn within_one_piece(law: &CantLaw, lo: Scalar, hi: Scalar) -> Option<(&CantPiece, Scalar)> {
    let (piece, s_lo) = piece_at(law, lo, false)?;
    (s_lo + (hi - lo) <= piece.length * (1.0 + 1e-12)).then_some((piece, s_lo))
}

/// Certified `(sup |D'|, sup |D''|, sup |D'''|)` of a banked curve's cant
/// over `[lo, hi]` inside one of its pieces: a height piece's own, an
/// angle piece's through `D = b sin(psi)` (module documentation).
fn cant_bounds(curve: &Banked3, lo: Scalar, hi: Scalar) -> Option<[Scalar; 3]> {
    let b = curve.rail_head_distance;
    if !(b.is_finite() && b > 0.0) {
        return None;
    }
    let (piece, s_lo) = within_one_piece(&curve.cant, lo, hi)?;
    let out = if let CantForm::VienneseBend { change, .. } = piece.form {
        let length = piece.length;
        let c = change.abs();
        let span = (s_lo / length, (s_lo + (hi - lo)) / length);
        let p1 = c * viennese_sup(1, span) / length;
        let p2 = c * viennese_sup(2, span) / (length * length);
        let p3 = c * viennese_sup(3, span) / (length * length * length);
        [
            b * p1,
            b * (p2 + p1 * p1),
            b * (p3 + 3.0 * p1 * p2 + p1 * p1 * p1),
        ]
        .map(|v| v * (1.0 + crate::bound::ROUNDING))
    } else {
        height_bounds(piece, s_lo, s_lo + (hi - lo))?
    };
    out.iter().all(|v| v.is_finite()).then_some(out)
}

/// `sup |f^(k)|` over `xi` in `[lo, hi]` (clamped to `[0, 1]`) of the
/// Viennese bend's shape `f = xi^4 (35 - 84 xi + 70 xi^2 - 20 xi^3)`, for
/// `k` in `1..=3`:
/// `f' = 140 w^3`, `f'' = 420 w^2 (1 - 2 xi)`, `f''' = 840 w (1 - 5 w)`
/// with `w = xi (1 - xi)`. Each is smooth, so its supremum over the span
/// is at an end or at a critical point inside: `1/2` for `f'`,
/// `1/2 +- 1 / sqrt(20)` for `f''`, `1/2` and `1/2 +- sqrt(0.15)` for
/// `f'''`. Over the whole piece that is `140 / 64`, `16.8 / sqrt(5)` and
/// `52.5`. The values are inflated by [`crate::bound::ROUNDING`] relative
/// and `1e-15` absolute for the rounding of their own evaluation.
fn viennese_sup(k: usize, (lo, hi): (Scalar, Scalar)) -> Scalar {
    let (lo, hi) = (lo.clamp(0.0, 1.0), hi.clamp(0.0, 1.0));
    let derivative = |xi: Scalar| -> Scalar {
        let w = xi * (1.0 - xi);
        match k {
            1 => 140.0 * w * w * w,
            2 => 420.0 * w * w * (1.0 - 2.0 * xi),
            _ => 840.0 * w * (1.0 - 5.0 * w),
        }
    };
    let critical: &[Scalar] = match k {
        1 => &[0.5],
        2 => &[0.5 - 0.223_606_797_749_979, 0.5 + 0.223_606_797_749_979],
        _ => &[
            0.5 - 0.387_298_334_620_741_7,
            0.5,
            0.5 + 0.387_298_334_620_741_7,
        ],
    };
    let sup = critical
        .iter()
        .copied()
        .filter(|&xi| xi > lo && xi < hi)
        .chain([lo, hi])
        .map(|xi| derivative(xi).abs())
        .fold(0.0, Scalar::max);
    sup * (1.0 + crate::bound::ROUNDING) + 1e-15
}

/// Certified `(sup |v'|, sup |v''|, sup |v'''|)` of a height piece over
/// `[s_lo, s_hi]` into it; `None` for an angle or held-rail piece.
///
/// In the piece's `xi = s / L`, the `k`-th derivative in plan distance is
/// the form's `k`-th in `xi` over `L^k`: a polynomial's by
/// `sum i!/(i-k)! |c_i| m^(i-k)` with `m` the larger `|xi|`, the
/// half-cosine's `|dD| / 2 pi^k`, the sine transition's `2 |dD|`,
/// `2 pi |dD|`, `4 pi^2 |dD|`.
fn height_bounds(piece: &CantPiece, s_lo: Scalar, s_hi: Scalar) -> Option<[Scalar; 3]> {
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

/// The cant law's seams inside or bounding a held-rail piece of the
/// pivot law, ascending: where the pivot derived from the cant may lose
/// smoothness (#279). Empty when the pivot has no held-rail piece.
pub(crate) fn derived_cant_seams(curve: &Banked3) -> Vec<Scalar> {
    let mut spans = Vec::new();
    let mut start = 0.0;
    for piece in &curve.pivot.pieces {
        let end = start + piece.length;
        if matches!(piece.form, CantForm::AboutRail { .. }) {
            spans.push((start, end));
        }
        start = end;
    }
    if spans.is_empty() {
        return Vec::new();
    }
    curve
        .cant
        .seams()
        .into_iter()
        .filter(|&seam| spans.iter().any(|&(a, b)| seam >= a && seam <= b))
        .collect()
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
