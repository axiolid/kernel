//! A centreline carrying a roll about its tangent: cant, bank angle and the
//! point the section rotates about.
//!
//! A track is authored as three independent laws over plan distance: the
//! plan, the vertical profile ([`Elevated3`]), and a cant law giving how far
//! one rail head stands above the other. [`Banked3`] pairs the first two
//! with the third, so the section frame an object is placed by is part of
//! the curve value instead of something every consumer re-derives.
//!
//! # The laws
//!
//! [`CantLaw`] is a run of [`CantPiece`]s laid end to end from plan distance
//! zero. Each piece is written in its OWN normalised coordinate
//! `xi = s / length`, `s` being the distance into the piece, because that is
//! how transition forms are tabulated: a Bloss transition is
//! `D1 + dD (3 - 2 xi) xi^2` whatever its length. The forms ([`CantForm`])
//! give the cant `D` directly (polynomial in `xi`, half-cosine, sine), or the
//! bank angle `psi` directly (the seventh-degree Viennese bend). A law value
//! says which it is ([`CantValue`]).
//!
//! The rail-head distance `b` relates them: `psi = asin(D / b)` and
//! `D = b sin(psi)`. A cant larger than `b` has no bank angle and is refused
//! by name ([`BankError::CantExceedsRailHeadDistance`]), never clamped.
//!
//! The pivot law `e(d)` is a second [`CantLaw`] of height pieces: the
//! elevation of the point the section rotates about above the profile. Zero
//! rotates about the centreline, `D / 2` (for a left-high cant) about the
//! low rail. The banked curve's point is `base(d) + e(d) z`.
//!
//! # Positive cant raises the LEFT rail
//!
//! `D` is the left rail head's height minus the right one's, left being
//! the side to the left of the direction of travel along the plan. A
//! positive bank angle rolls the section anticlockwise seen from behind.
//!
//! # Two conventions, named, never defaulted
//!
//! On a grade, "the rail heads differ by `D`" is ambiguous, and the
//! difference is systematic, not noise. [`BankConvention`] names the two
//! readings and a [`Banked3`] carries one; there is no default (ADR 0081).
//! With `t` the unit tangent of the banked curve, `theta` its grade angle,
//! `n` the horizontal left normal and `u = t x n` (the up direction square
//! to the tangent), the section is rolled about `t` by an angle `rho`:
//! lateral `l = cos(rho) n + sin(rho) u`, section up
//! `v = -sin(rho) n + cos(rho) u`. Both frames are exactly orthonormal and
//! the rail heads lie `b` apart along `l` in the plane square to `t`.
//!
//! - [`BankConvention::TangentRotation`]: `rho = psi = asin(D / b)`. The
//!   cross-section is the nominal level-track section turned about the
//!   tangent. The rail heads then differ in height by `D cos(theta)`.
//! - [`BankConvention::VerticalRise`]: `rho = asin(D / (b cos theta))`, the
//!   one roll for which the rail heads differ in height by exactly `D`.
//!   Refused when `|D| > b cos(theta)`.
//!
//! The rises differ by `D (1 - cos theta)`: `1.9994e-4 D` at a 2% grade,
//! 0.03 mm on a 150 mm cant. On level track the two agree exactly.

use core::fmt;

use axiolid_core::Scalar;

use crate::Elevated3;

/// The shape of one cant piece over its own `xi = s / length` in `[0, 1]`.
///
/// Height forms give the cant `D` (the left rail head above the right) in
/// length units; the angle form gives the bank angle `psi` in radians.
#[non_exhaustive]
#[derive(Debug, Clone, PartialEq)]
pub enum CantForm {
    /// `D(xi) = coefficients[0] + coefficients[1] xi + coefficients[2] xi^2 + ...`
    ///
    /// Covers a constant cant (degree 0), a linear transition (degree 1),
    /// the Bloss transition `D1 + dD (3 - 2 xi) xi^2` (degree 3) and each
    /// half of a Helmert transition (degree 2; see [`CantPiece::helmert`]).
    /// An empty list is zero.
    Polynomial {
        /// Coefficients in ascending powers of `xi`.
        coefficients: Vec<Scalar>,
    },
    /// `D(xi) = start + change (1 - cos(pi xi)) / 2`: the half-cosine
    /// transition, tangent to its neighbours' slopes at both ends.
    Cosine {
        /// Cant at `xi = 0`.
        start: Scalar,
        /// Cant gained over the piece.
        change: Scalar,
    },
    /// `D(xi) = start + change (xi - sin(2 pi xi) / (2 pi))`: the sine
    /// transition, whose slope AND curvature vanish at both ends.
    Sine {
        /// Cant at `xi = 0`.
        start: Scalar,
        /// Cant gained over the piece.
        change: Scalar,
    },
    /// `psi(xi) = start + change xi^4 (35 - 84 xi + 70 xi^2 - 20 xi^3)`: the
    /// Viennese bend, written for the bank angle itself. Its first three
    /// derivatives vanish at both ends. The cant follows as `b sin(psi)`.
    VienneseBend {
        /// Bank angle at `xi = 0`, radians.
        start: Scalar,
        /// Bank angle gained over the piece, radians.
        change: Scalar,
    },
}

/// What a cant law gives at a station: a cant height or a bank angle.
#[non_exhaustive]
#[derive(Debug, Clone, Copy, PartialEq)]
pub enum CantValue {
    /// Cant `D`, or a rate of it, in length units (per unit plan distance
    /// for a rate).
    Height(Scalar),
    /// Bank angle `psi`, or a rate of it, in radians (per unit plan
    /// distance for a rate).
    Angle(Scalar),
}

/// One piece of a [`CantLaw`]: a form over a plan length.
#[derive(Debug, Clone, PartialEq)]
pub struct CantPiece {
    /// Plan length the piece covers. Positive and finite when well formed.
    pub length: Scalar,
    /// The law over the piece's own `xi = s / length`.
    pub form: CantForm,
}

impl CantPiece {
    /// A piece of `form` over `length`.
    #[must_use]
    pub const fn new(length: Scalar, form: CantForm) -> Self {
        Self { length, form }
    }

    /// A constant cant.
    #[must_use]
    pub fn constant(length: Scalar, cant: Scalar) -> Self {
        Self::new(
            length,
            CantForm::Polynomial {
                coefficients: vec![cant],
            },
        )
    }

    /// A linear transition from `start` to `end`.
    #[must_use]
    pub fn linear(length: Scalar, start: Scalar, end: Scalar) -> Self {
        Self::new(
            length,
            CantForm::Polynomial {
                coefficients: vec![start, end - start],
            },
        )
    }

    /// A Bloss transition: `start + (end - start) (3 - 2 xi) xi^2`.
    #[must_use]
    pub fn bloss(length: Scalar, start: Scalar, end: Scalar) -> Self {
        let change = end - start;
        Self::new(
            length,
            CantForm::Polynomial {
                coefficients: vec![start, 0.0, 3.0 * change, -2.0 * change],
            },
        )
    }

    /// A Helmert (two-parabola) transition from `start` to `end` over
    /// `length`, as its two halves.
    ///
    /// Over the whole transition, with `x = s / length`, the cant is
    /// `start + 2 dD x^2` up to `x = 1/2` and `end - 2 dD (1 - x)^2` after.
    /// Each half is written in its own `xi`:
    /// `start + dD/2 xi^2`, then `start + dD/2 + dD xi - dD/2 xi^2`.
    #[must_use]
    pub fn helmert(length: Scalar, start: Scalar, end: Scalar) -> [Self; 2] {
        let change = end - start;
        let half = 0.5 * length;
        [
            Self::new(
                half,
                CantForm::Polynomial {
                    coefficients: vec![start, 0.0, 0.5 * change],
                },
            ),
            Self::new(
                half,
                CantForm::Polynomial {
                    coefficients: vec![start + 0.5 * change, change, -0.5 * change],
                },
            ),
        ]
    }

    /// A half-cosine transition from `start` to `end`.
    #[must_use]
    pub fn cosine(length: Scalar, start: Scalar, end: Scalar) -> Self {
        Self::new(
            length,
            CantForm::Cosine {
                start,
                change: end - start,
            },
        )
    }

    /// A sine transition from `start` to `end`.
    #[must_use]
    pub fn sine(length: Scalar, start: Scalar, end: Scalar) -> Self {
        Self::new(
            length,
            CantForm::Sine {
                start,
                change: end - start,
            },
        )
    }

    /// A Viennese bend from bank angle `start` to `end`, radians.
    #[must_use]
    pub fn viennese_bend(length: Scalar, start: Scalar, end: Scalar) -> Self {
        Self::new(
            length,
            CantForm::VienneseBend {
                start,
                change: end - start,
            },
        )
    }

    /// Value and rate (per unit PLAN distance) at `s` into the piece.
    fn sample(&self, s: Scalar) -> Option<(CantValue, CantValue)> {
        let length = self.length;
        if !(length.is_finite() && length > 0.0) {
            return None;
        }
        let xi = s / length;
        let (value, slope, angle) = match &self.form {
            CantForm::Polynomial { coefficients } => {
                let value = coefficients
                    .iter()
                    .rev()
                    .fold(0.0, |accumulated, c| accumulated * xi + c);
                let slope = coefficients
                    .iter()
                    .enumerate()
                    .skip(1)
                    .rev()
                    .fold(0.0, |accumulated, (power, c)| {
                        accumulated * xi + c * power as Scalar
                    });
                (value, slope, false)
            }
            CantForm::Cosine { start, change } => {
                let phase = core::f64::consts::PI * xi;
                (
                    start + change * 0.5 * (1.0 - phase.cos()),
                    change * 0.5 * core::f64::consts::PI * phase.sin(),
                    false,
                )
            }
            CantForm::Sine { start, change } => {
                let phase = core::f64::consts::TAU * xi;
                (
                    start + change * (xi - phase.sin() / core::f64::consts::TAU),
                    change * (1.0 - phase.cos()),
                    false,
                )
            }
            CantForm::VienneseBend { start, change } => {
                let xi2 = xi * xi;
                let shape = xi2 * xi2 * (35.0 + xi * (-84.0 + xi * (70.0 - 20.0 * xi)));
                let one_minus = 1.0 - xi;
                let rate = 140.0 * xi2 * xi * one_minus * one_minus * one_minus;
                (start + change * shape, change * rate, true)
            }
        };
        let rate = slope / length;
        if !(value.is_finite() && rate.is_finite()) {
            return None;
        }
        Some(if angle {
            (CantValue::Angle(value), CantValue::Angle(rate))
        } else {
            (CantValue::Height(value), CantValue::Height(rate))
        })
    }

    fn is_well_formed(&self) -> bool {
        let finite = match &self.form {
            CantForm::Polynomial { coefficients } => coefficients.iter().all(|c| c.is_finite()),
            CantForm::Cosine { start, change }
            | CantForm::Sine { start, change }
            | CantForm::VienneseBend { start, change } => start.is_finite() && change.is_finite(),
        };
        finite && self.length.is_finite() && self.length > 0.0
    }

    fn is_angle(&self) -> bool {
        matches!(self.form, CantForm::VienneseBend { .. })
    }
}

/// A scalar law over plan distance: pieces laid end to end from zero.
///
/// Piece `i` spans `[sum of earlier lengths, + length_i)`; a seam belongs to
/// the piece that starts there, and the last piece is closed at its far end.
/// A distance before zero or past the last piece has no value: the law is
/// not extrapolated, because an extrapolated transition is not a cant
/// anyone authored.
#[derive(Debug, Clone, PartialEq, Default)]
pub struct CantLaw {
    /// Pieces in order of plan distance.
    pub pieces: Vec<CantPiece>,
}

impl CantLaw {
    /// A law from its pieces.
    #[must_use]
    pub const fn new(pieces: Vec<CantPiece>) -> Self {
        Self { pieces }
    }

    /// Zero over `length`: the pivot law of a rotation about the
    /// centreline.
    #[must_use]
    pub fn zero(length: Scalar) -> Self {
        Self::new(vec![CantPiece::constant(length, 0.0)])
    }

    /// Total plan length the law covers.
    #[must_use]
    pub fn length(&self) -> Scalar {
        self.pieces.iter().map(|piece| piece.length).sum()
    }

    /// Interior seams, in plan distance, ascending.
    #[must_use]
    pub fn seams(&self) -> Vec<Scalar> {
        let mut at = 0.0;
        let mut seams = Vec::with_capacity(self.pieces.len().saturating_sub(1));
        for piece in self.pieces.iter().take(self.pieces.len().saturating_sub(1)) {
            at += piece.length;
            seams.push(at);
        }
        seams
    }

    /// Whether every piece has a positive finite length and finite
    /// coefficients, and there is at least one piece.
    #[must_use]
    pub fn is_well_formed(&self) -> bool {
        !self.pieces.is_empty() && self.pieces.iter().all(CantPiece::is_well_formed)
    }

    /// Whether any piece gives a bank angle rather than a height.
    #[must_use]
    pub fn has_angle_pieces(&self) -> bool {
        self.pieces.iter().any(CantPiece::is_angle)
    }

    /// Value at plan distance `distance`.
    ///
    /// `None` when the law is malformed or `distance` is outside it.
    #[must_use]
    pub fn value_at(&self, distance: Scalar) -> Option<CantValue> {
        self.sample(distance).map(|(value, _)| value)
    }

    /// Rate of change per unit plan distance at `distance`, from the piece
    /// that owns it (the one starting there, at a seam).
    #[must_use]
    pub fn rate_at(&self, distance: Scalar) -> Option<CantValue> {
        self.sample(distance).map(|(_, rate)| rate)
    }

    fn sample(&self, distance: Scalar) -> Option<(CantValue, CantValue)> {
        if !distance.is_finite() || distance < 0.0 || !self.is_well_formed() {
            return None;
        }
        let mut start = 0.0;
        let last = self.pieces.len() - 1;
        for (index, piece) in self.pieces.iter().enumerate() {
            let end = start + piece.length;
            if distance < end || (index == last && distance <= end) {
                return piece.sample(distance - start);
            }
            start = end;
        }
        None
    }
}

/// Which height difference a cant `D` states on a grade (ADR 0081).
///
/// See the [module docs](self) for both frames. No variant is a default:
/// a [`Banked3`] names one.
#[non_exhaustive]
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub enum BankConvention {
    /// The level-track section turned about the 3D tangent by
    /// `psi = asin(D / b)`. The rail heads differ in height by
    /// `D cos(theta)` at grade angle `theta`.
    TangentRotation,
    /// The section rolled about the 3D tangent by `asin(D / (b cos theta))`,
    /// so the rail heads, `b` apart square to the tangent, differ in height
    /// by exactly `D`.
    VerticalRise,
}

impl BankConvention {
    /// The roll about the tangent, in radians, for cant `cant` over
    /// rail-head distance `rail_head_distance` at a grade whose angle has
    /// cosine `grade_cosine`.
    ///
    /// # Errors
    ///
    /// [`BankError::RailHeadDistance`] for a non-positive or non-finite
    /// `b`; [`BankError::CantExceedsRailHeadDistance`] when `|D| > b`;
    /// [`BankError::CantExceedsVerticalSpan`] when, under
    /// [`VerticalRise`](Self::VerticalRise), `|D| > b cos(theta)`;
    /// [`BankError::NonFinite`] for a non-finite cant or cosine.
    pub fn roll(
        self,
        cant: Scalar,
        rail_head_distance: Scalar,
        grade_cosine: Scalar,
    ) -> Result<Scalar, BankError> {
        let psi = bank_angle(cant, rail_head_distance)?;
        self.roll_from_angle(psi, rail_head_distance, grade_cosine)
    }

    /// [`roll`](Self::roll) from a nominal bank angle `psi` with
    /// `sin(psi) = D / b`, as an angle-form piece gives it.
    ///
    /// Under [`TangentRotation`](Self::TangentRotation) this is `psi`
    /// itself, with no round trip through the cant.
    ///
    /// # Errors
    ///
    /// As [`roll`](Self::roll), and [`BankError::AngleOutOfRange`] for
    /// `|psi| > pi / 2`, which no cant gives.
    pub fn roll_from_angle(
        self,
        psi: Scalar,
        rail_head_distance: Scalar,
        grade_cosine: Scalar,
    ) -> Result<Scalar, BankError> {
        check_rail_head_distance(rail_head_distance)?;
        if !psi.is_finite() || !grade_cosine.is_finite() {
            return Err(BankError::NonFinite);
        }
        if psi.abs() > core::f64::consts::FRAC_PI_2 {
            return Err(BankError::AngleOutOfRange { angle: psi });
        }
        match self {
            Self::TangentRotation => Ok(psi),
            Self::VerticalRise => {
                let cant = rail_head_distance * psi.sin();
                let span = rail_head_distance * grade_cosine;
                if span <= 0.0 || cant.abs() > span {
                    return Err(BankError::CantExceedsVerticalSpan { cant, span });
                }
                Ok((cant / span).asin())
            }
        }
    }
}

/// The nominal bank angle `psi = asin(D / b)`.
///
/// # Errors
///
/// [`BankError::RailHeadDistance`] for a non-positive or non-finite `b`,
/// [`BankError::NonFinite`] for a non-finite cant, and
/// [`BankError::CantExceedsRailHeadDistance`] when `|D| > b`.
pub fn bank_angle(cant: Scalar, rail_head_distance: Scalar) -> Result<Scalar, BankError> {
    check_rail_head_distance(rail_head_distance)?;
    if !cant.is_finite() {
        return Err(BankError::NonFinite);
    }
    if cant.abs() > rail_head_distance {
        return Err(BankError::CantExceedsRailHeadDistance {
            cant,
            rail_head_distance,
        });
    }
    Ok((cant / rail_head_distance).asin())
}

fn check_rail_head_distance(rail_head_distance: Scalar) -> Result<(), BankError> {
    if rail_head_distance.is_finite() && rail_head_distance > 0.0 {
        Ok(())
    } else {
        Err(BankError::RailHeadDistance(rail_head_distance))
    }
}

/// Why a banked curve has no section at a station.
#[non_exhaustive]
#[derive(Debug, Clone, Copy, PartialEq)]
pub enum BankError {
    /// The rail-head distance is not positive and finite.
    RailHeadDistance(Scalar),
    /// `|D| > b`: the rail heads cannot differ by more than they are apart.
    CantExceedsRailHeadDistance {
        /// The cant asked for.
        cant: Scalar,
        /// The rail-head distance `b`.
        rail_head_distance: Scalar,
    },
    /// Under [`BankConvention::VerticalRise`], `|D| > b cos(theta)`: the
    /// grade leaves the rail heads too little vertical reach for the cant.
    CantExceedsVerticalSpan {
        /// The cant asked for.
        cant: Scalar,
        /// `b cos(theta)`.
        span: Scalar,
    },
    /// A bank angle beyond a quarter turn.
    AngleOutOfRange {
        /// The angle, radians.
        angle: Scalar,
    },
    /// The named law is malformed or has no value at the distance.
    OutsideLaw {
        /// `"cant"` or `"pivot"`.
        law: &'static str,
        /// The plan distance asked for.
        distance: Scalar,
    },
    /// The pivot law has an angle-form piece; a pivot is an elevation.
    AngleInPivot,
    /// A non-finite cant, angle or grade.
    NonFinite,
}

impl fmt::Display for BankError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::RailHeadDistance(b) => {
                write!(f, "rail-head distance {b} is not positive and finite")
            }
            Self::CantExceedsRailHeadDistance {
                cant,
                rail_head_distance,
            } => write!(
                f,
                "cant {cant} exceeds the rail-head distance {rail_head_distance}, so it has no bank angle"
            ),
            Self::CantExceedsVerticalSpan { cant, span } => write!(
                f,
                "cant {cant} exceeds the vertical span {span} of the rail heads on this grade"
            ),
            Self::AngleOutOfRange { angle } => {
                write!(f, "bank angle {angle} is beyond a quarter turn")
            }
            Self::OutsideLaw { law, distance } => write!(
                f,
                "{law} law is malformed or has no value at plan distance {distance}"
            ),
            Self::AngleInPivot => {
                write!(f, "pivot law has an angle-form piece; a pivot is an elevation")
            }
            Self::NonFinite => write!(f, "cant, bank angle or grade is not finite"),
        }
    }
}

impl std::error::Error for BankError {}

/// A centreline with a roll law: the composition a canted track is
/// authored as. See the [module docs](self).
///
/// Evaluated by PLAN distance `d`, like its base. The point is
/// `base(d) + pivot(d) z`; the section frame is rolled about that point's
/// own unit tangent by the angle [`convention`](Self::convention) gives.
/// The curve's span is the cant law's: no station outside it has a section.
#[derive(Debug, Clone, PartialEq)]
pub struct Banked3 {
    /// Centreline: plan and vertical profile.
    pub base: Elevated3,
    /// Cant `D(d)` (left rail head above right), or bank angle, by piece.
    pub cant: CantLaw,
    /// Elevation `e(d)` of the rotation point above the profile. Height
    /// pieces only.
    pub pivot: CantLaw,
    /// Rail-head distance `b`, from which `psi = asin(D / b)`.
    pub rail_head_distance: Scalar,
    /// Which height difference `D` states on a grade.
    pub convention: BankConvention,
}

impl Banked3 {
    /// Pair a centreline with cant and pivot laws under a named
    /// convention.
    #[must_use]
    pub const fn new(
        base: Elevated3,
        cant: CantLaw,
        pivot: CantLaw,
        rail_head_distance: Scalar,
        convention: BankConvention,
    ) -> Self {
        Self {
            base,
            cant,
            pivot,
            rail_head_distance,
            convention,
        }
    }

    /// Plan length the curve is defined over: the cant law's.
    #[must_use]
    pub fn span(&self) -> Scalar {
        self.cant.length()
    }

    /// Cant `D` at `distance`: the law's height, or `b sin(psi)` for an
    /// angle piece.
    ///
    /// # Errors
    ///
    /// [`BankError::OutsideLaw`] off the law, [`BankError::RailHeadDistance`]
    /// for a bad `b`, and the angle checks of [`bank_angle_at`](Self::bank_angle_at).
    pub fn cant_at(&self, distance: Scalar) -> Result<Scalar, BankError> {
        check_rail_head_distance(self.rail_head_distance)?;
        match self.cant_value(distance)? {
            CantValue::Height(cant) => {
                bank_angle(cant, self.rail_head_distance)?;
                Ok(cant)
            }
            CantValue::Angle(psi) => {
                check_angle(psi)?;
                Ok(self.rail_head_distance * psi.sin())
            }
        }
    }

    /// Nominal bank angle `psi = asin(D / b)` at `distance`; an angle piece
    /// gives it directly.
    ///
    /// # Errors
    ///
    /// [`BankError::CantExceedsRailHeadDistance`] when `|D| > b`,
    /// [`BankError::AngleOutOfRange`] for an angle piece beyond a quarter
    /// turn, [`BankError::OutsideLaw`] off the law.
    pub fn bank_angle_at(&self, distance: Scalar) -> Result<Scalar, BankError> {
        match self.cant_value(distance)? {
            CantValue::Height(cant) => bank_angle(cant, self.rail_head_distance),
            CantValue::Angle(psi) => {
                check_rail_head_distance(self.rail_head_distance)?;
                check_angle(psi)?;
                Ok(psi)
            }
        }
    }

    /// Elevation `e` of the rotation point above the profile at `distance`,
    /// and its rate per unit plan distance.
    ///
    /// # Errors
    ///
    /// [`BankError::AngleInPivot`] for a pivot law with an angle piece,
    /// [`BankError::OutsideLaw`] off the law.
    pub fn pivot_at(&self, distance: Scalar) -> Result<(Scalar, Scalar), BankError> {
        if self.pivot.has_angle_pieces() {
            return Err(BankError::AngleInPivot);
        }
        let outside = BankError::OutsideLaw {
            law: "pivot",
            distance,
        };
        match (self.pivot.value_at(distance), self.pivot.rate_at(distance)) {
            (Some(CantValue::Height(value)), Some(CantValue::Height(rate))) => Ok((value, rate)),
            _ => Err(outside),
        }
    }

    fn cant_value(&self, distance: Scalar) -> Result<CantValue, BankError> {
        self.cant.value_at(distance).ok_or(BankError::OutsideLaw {
            law: "cant",
            distance,
        })
    }
}

fn check_angle(psi: Scalar) -> Result<(), BankError> {
    if !psi.is_finite() {
        return Err(BankError::NonFinite);
    }
    if psi.abs() > core::f64::consts::FRAC_PI_2 {
        return Err(BankError::AngleOutOfRange { angle: psi });
    }
    Ok(())
}
