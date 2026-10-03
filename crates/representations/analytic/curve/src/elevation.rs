//! Elevation laws and the composition of a planar curve with one.
//!
//! A road or rail centreline is authored as two independent laws: a planar
//! layout carrying the transition spirals, and a vertical profile giving
//! height as a function of distance along that layout. The 3D centreline is
//! their composition. This module holds the vertical half and the pairing;
//! the planar half is an ordinary [`Curve2`].
//!
//! # Height is a function of PLAN distance, not 3D arc length
//!
//! The vertical profile a surveyor writes is indexed by chainage measured on
//! the horizontal layout, so `height_at` takes plan distance. The two differ
//! whenever the grade is non-zero, because `ds3 = sqrt(1 + g^2) ds_plan`.
//! Over a 120 m curve at a 2% entry grade the 3D length exceeds the plan
//! length by 12.5 mm, which reads the profile 0.35 mm off if the distances
//! are confused. That is small but it is a systematic misreading, not noise,
//! and it compounds along a chain of segments.
//!
//! Naming the convention here means a consumer never has to guess which
//! distance a law is written in.
//!
//! # Closed forms here, quadrature in an evaluator
//!
//! A polynomial and a circular arc in the `(d, z)` plane have closed-form
//! heights and grades, so [`ElevationLaw::height_at`] and
//! [`ElevationLaw::grade_at`] answer them here. A profile given by its
//! curvature against its own arc length ([`ElevationLaw::Intrinsic`], the
//! clothoid between two grades) does not: its height is the Fresnel-type
//! integral the [`CurvatureLaw`] module refuses to compute, followed by an
//! inversion of plan distance against arc length. This crate stores that
//! law exactly and reports `None` for it; an evaluator that can state its
//! quadrature tolerance reads it (`axiolid-evaluate`'s `elevation` module).

use axiolid_core::Scalar;

use crate::{CurvatureLaw, Curve2};

/// Height as a function of distance along the plan.
///
/// Mirrors [`CurvatureLaw`] in shape so the two halves of
/// an alignment read the same way, but stays a separate type: a curvature law
/// is a property of a planar curve and an elevation law is not, and sharing one
/// enum would make a meaningless pairing representable.
///
/// Dirty imported data stays representable, as everywhere else in this crate.
/// Mismatched piece lists report `None` rather than guessing.
#[non_exhaustive]
#[derive(Debug, Clone, PartialEq)]
pub enum ElevationLaw {
    /// `z(d) = coefficients[0] + coefficients[1] * d + coefficients[2] * d^2 + ...`
    ///
    /// Degree 1 is a constant gradient, degree 2 the parabolic vertical curve
    /// used to join two grades. Those are the two vertical segment kinds that
    /// carry most alignment data, and both are exact here.
    ///
    /// An empty coefficient list is the zero polynomial: height zero.
    Polynomial {
        /// Coefficients in ascending powers of plan distance.
        coefficients: Vec<Scalar>,
    },
    /// Pieces laid end to end along plan distance, each with its own law.
    ///
    /// `breaks` holds the INTERIOR seam positions measured from the start, so
    /// `laws.len() == breaks.len() + 1` and piece `i` spans
    /// `breaks[i - 1] .. breaks[i]`.
    ///
    /// Each piece's law is written in its OWN distance, restarting at zero at
    /// its seam, so moving a piece never rewrites its coefficients. This
    /// matches [`CurvatureLaw::Piecewise`]
    /// deliberately: a vertical profile is authored as a run of segments and
    /// the seams are observable data.
    Piecewise {
        /// Interior seam positions in plan distance, ascending.
        breaks: Vec<Scalar>,
        /// One law per piece; `laws.len() == breaks.len() + 1`.
        laws: Vec<ElevationLaw>,
    },
    /// A circular arc in the `(d, z)` plane: the grade angle changes at a
    /// constant rate along the profile's own arc length.
    ///
    /// With `t0 = atan(grade)` and the signed `radius` `R` (positive is a
    /// sag, turning counter-clockwise in `(d, z)`; negative is a crest):
    ///
    /// - `sin t(d) = sin t0 + d / R`
    /// - `z(d) = height + R (cos t0 - cos t(d))`
    /// - `grade(d) = tan t(d)`
    ///
    /// This is the circle itself, not the parabola `d^2 / (2R)` that
    /// approximates it near the vertex; the two differ by `O(d^4 / R^3)`.
    /// The arc has heights only while `|sin t0 + d / R| < 1`: past that the
    /// circle turns vertical, and `height_at` and `grade_at` report `None`.
    ///
    /// Evaluated without cancellation: `R (cos t0 - cos t(d))` is computed
    /// as `d (sin t(d) + sin t0) / (cos t0 + cos t(d))`, which removes the
    /// difference of nearly equal cosines that `d / R -> 0` would otherwise
    /// multiply by a large `R`.
    CircularArc {
        /// Height at the piece start.
        height: Scalar,
        /// Grade `dz/dd` at the piece start.
        grade: Scalar,
        /// Signed radius: positive sag, negative crest. Finite and non-zero.
        radius: Scalar,
    },
    /// A profile given by its curvature against its OWN arc length.
    ///
    /// The vertical counterpart of [`Curve2::Intrinsic`]: a planar curve in
    /// `(d, z)` starting at `(0, height)` in the direction `atan(grade)`,
    /// whose curvature is `curvature(s)` with `s` the arc length of the
    /// profile (not plan distance), positive counter-clockwise in `(d, z)`
    /// -- a sag. A linear law is the clothoid joining two grades; a constant
    /// law is the circular arc that [`Self::CircularArc`] holds in closed
    /// form.
    ///
    /// Height at a plan distance has no closed form (see the module
    /// documentation), so `height_at` and `grade_at` report `None`, and an
    /// evaluator integrates the curve and inverts `d(s)` to a stated
    /// tolerance. Only `d >= 0` is defined, since the law is integrated
    /// forward from its start, and only while `d(s)` keeps increasing.
    Intrinsic {
        /// Height at the piece start.
        height: Scalar,
        /// Grade `dz/dd` at the piece start.
        grade: Scalar,
        /// Curvature against the profile's own arc length from its start.
        curvature: CurvatureLaw,
    },
}

impl ElevationLaw {
    /// A constant height.
    #[must_use]
    pub fn level(height: Scalar) -> Self {
        Self::Polynomial {
            coefficients: vec![height],
        }
    }

    /// A constant gradient: `z(d) = height + grade * d`.
    #[must_use]
    pub fn constant_grade(height: Scalar, grade: Scalar) -> Self {
        Self::Polynomial {
            coefficients: vec![height, grade],
        }
    }

    /// A parabolic vertical curve joining `entry_grade` to `exit_grade` over
    /// `length`.
    ///
    /// The rate of change of grade is `(exit - entry) / length`, so
    /// `z(d) = height + entry * d + (exit - entry) / (2 * length) * d^2`.
    /// A non-positive length is storable; naming it is a validator's job.
    #[must_use]
    pub fn parabolic(
        height: Scalar,
        entry_grade: Scalar,
        exit_grade: Scalar,
        length: Scalar,
    ) -> Self {
        Self::Polynomial {
            coefficients: vec![
                height,
                entry_grade,
                (exit_grade - entry_grade) / (2.0 * length),
            ],
        }
    }

    /// A circular vertical arc from `height` at `grade`, of signed `radius`
    /// (positive sag, negative crest). See [`Self::CircularArc`].
    #[must_use]
    pub const fn circular_arc(height: Scalar, grade: Scalar, radius: Scalar) -> Self {
        Self::CircularArc {
            height,
            grade,
            radius,
        }
    }

    /// A profile from `height` at `grade` whose curvature runs along its own
    /// arc length. See [`Self::Intrinsic`].
    #[must_use]
    pub const fn intrinsic(height: Scalar, grade: Scalar, curvature: CurvatureLaw) -> Self {
        Self::Intrinsic {
            height,
            grade,
            curvature,
        }
    }

    /// Whether the piece lists agree and the seams ascend, and every stored
    /// number is finite (an arc's radius also non-zero).
    #[must_use]
    pub fn is_well_formed(&self) -> bool {
        match self {
            Self::Polynomial { coefficients } => coefficients.iter().all(|c| c.is_finite()),
            Self::CircularArc {
                height,
                grade,
                radius,
            } => height.is_finite() && grade.is_finite() && radius.is_finite() && *radius != 0.0,
            Self::Intrinsic {
                height,
                grade,
                curvature,
            } => height.is_finite() && grade.is_finite() && curvature.is_well_formed(),
            Self::Piecewise { breaks, laws } => {
                laws.len() == breaks.len() + 1
                    && breaks.iter().all(|b| b.is_finite())
                    && breaks.windows(2).all(|w| w[0] <= w[1])
                    && laws.iter().all(Self::is_well_formed)
            }
        }
    }

    /// Height at `distance` along the plan.
    ///
    /// Returns `None` when the law is malformed or the distance is not finite,
    /// rather than extrapolating off a piece that does not exist.
    #[must_use]
    pub fn height_at(&self, distance: Scalar) -> Option<Scalar> {
        if !distance.is_finite() {
            return None;
        }
        match self {
            Self::Polynomial { coefficients } => {
                Some(horner(coefficients, distance)).filter(|z| z.is_finite())
            }
            Self::Piecewise { breaks, laws } => {
                let (law, local) = piece_at(breaks, laws, distance)?;
                law.height_at(local)
            }
            Self::CircularArc {
                height,
                grade,
                radius,
            } => {
                let arc = ArcState::at(*grade, *radius, distance)?;
                // R (cos t0 - cos t) = R (sin^2 t - sin^2 t0) / (cos t0 + cos t)
                // and R (sin t - sin t0) = d, so the large R cancels exactly.
                Some(height + distance * (arc.sin + arc.sin0) / (arc.cos0 + arc.cos))
                    .filter(|z| z.is_finite())
            }
            // No closed form: an evaluator integrates it (module docs).
            Self::Intrinsic { .. } => None,
        }
    }

    /// Grade -- `dz/dd` -- at `distance` along the plan.
    ///
    /// This is the slope the 3D tangent needs, and it is why the polynomial is
    /// differentiated exactly rather than differenced.
    #[must_use]
    pub fn grade_at(&self, distance: Scalar) -> Option<Scalar> {
        if !distance.is_finite() {
            return None;
        }
        match self {
            Self::Polynomial { coefficients } => {
                let derivative: Vec<Scalar> = coefficients
                    .iter()
                    .enumerate()
                    .skip(1)
                    .map(|(power, c)| c * power as Scalar)
                    .collect();
                Some(horner(&derivative, distance)).filter(|g| g.is_finite())
            }
            Self::Piecewise { breaks, laws } => {
                let (law, local) = piece_at(breaks, laws, distance)?;
                law.grade_at(local)
            }
            Self::CircularArc { grade, radius, .. } => {
                let arc = ArcState::at(*grade, *radius, distance)?;
                Some(arc.sin / arc.cos).filter(|g| g.is_finite())
            }
            // No closed form: an evaluator integrates it (module docs).
            Self::Intrinsic { .. } => None,
        }
    }

    /// The innermost piece covering `distance`, and that distance rebased
    /// to the piece's own start.
    ///
    /// Descends through nested [`Self::Piecewise`] laws, so the returned law
    /// is never piecewise; any other law is its own piece at the unchanged
    /// distance. Seams belong to the piece that starts there, as in
    /// `height_at`. `None` when the distance is not finite or a piece list
    /// is malformed.
    ///
    /// This is what an evaluator needs to read a law it evaluates itself
    /// (an [`Self::Intrinsic`] piece) inside a composed profile.
    #[must_use]
    pub fn piece_at(&self, distance: Scalar) -> Option<(&Self, Scalar)> {
        if !distance.is_finite() {
            return None;
        }
        match self {
            Self::Piecewise { breaks, laws } => {
                let (law, local) = piece_at(breaks, laws, distance)?;
                law.piece_at(local)
            }
            _ => Some((self, distance)),
        }
    }
}

/// The trigonometry of a circular vertical arc at one plan distance.
struct ArcState {
    /// `sin t0` of the start direction.
    sin0: Scalar,
    /// `cos t0` of the start direction.
    cos0: Scalar,
    /// `sin t(d)` at the distance.
    sin: Scalar,
    /// `cos t(d)` at the distance, positive inside the domain.
    cos: Scalar,
}

impl ArcState {
    /// `None` outside the arc's domain `|sin t0 + d / R| < 1`, or for a
    /// non-finite grade or a zero or non-finite radius.
    fn at(grade: Scalar, radius: Scalar, distance: Scalar) -> Option<Self> {
        if !(grade.is_finite() && radius.is_finite() && radius != 0.0) {
            return None;
        }
        // sin and cos of atan(grade) without forming the angle: no rounding
        // through atan and back, and exact at grade zero.
        let norm = grade.hypot(1.0);
        let sin0 = grade / norm;
        let cos0 = 1.0 / norm;
        let sin = sin0 + distance / radius;
        if sin.is_nan() || sin.abs() >= 1.0 {
            return None;
        }
        // (1 - s)(1 + s) keeps its relative accuracy as |s| -> 1, where
        // 1 - s^2 would not.
        let cos = ((1.0 - sin) * (1.0 + sin)).sqrt();
        if cos <= 0.0 {
            return None;
        }
        Some(Self {
            sin0,
            cos0,
            sin,
            cos,
        })
    }
}

/// Evaluate ascending-power coefficients at `x`.
fn horner(coefficients: &[Scalar], x: Scalar) -> Scalar {
    coefficients
        .iter()
        .rev()
        .fold(0.0, |accumulated, c| accumulated * x + c)
}

/// The piece covering `distance`, and that distance rebased to the piece start.
///
/// Pieces are half-open so a seam belongs to the piece that starts there; the
/// final piece is closed at its far end so the curve's last point evaluates.
fn piece_at<'a>(
    breaks: &[Scalar],
    laws: &'a [ElevationLaw],
    distance: Scalar,
) -> Option<(&'a ElevationLaw, Scalar)> {
    if laws.len() != breaks.len() + 1 {
        return None;
    }
    let index = breaks.partition_point(|b| *b <= distance);
    let start = if index == 0 { 0.0 } else { breaks[index - 1] };
    laws.get(index).map(|law| (law, distance - start))
}

/// A planar curve carrying an independent elevation law.
///
/// This is the composition an alignment centreline actually is: the plan is
/// exact -- including the transition spirals a [`Curve2::Intrinsic`] holds --
/// and the vertical profile is exact, and neither is approximated to pair
/// them. Evaluation is parameterised by PLAN DISTANCE, which is the parameter
/// both halves are authored against.
///
/// Composition rather than re-encoding: a `Curve3::BSpline` fitted through the
/// pair would lose both. Storing the two laws keeps each one's own exactness
/// and lets a consumer recover either half unchanged.
///
/// Only a plan whose parameter is arc length (line, circle, intrinsic,
/// [`Chain2`](crate::Chain2)) can carry a law. A B-spline's parameter is not
/// a distance, so evaluators refuse an elevated B-spline plan rather than
/// reinterpret the law (ADR 0060); a B-spline read by arc length is a
/// chain piece.
#[derive(Debug, Clone, PartialEq)]
pub struct Elevated3 {
    /// Horizontal layout. Boxed to keep [`Curve3`](crate::Curve3) small.
    pub plan: Box<Curve2>,
    /// Height as a function of distance along `plan`.
    pub elevation: ElevationLaw,
}

impl Elevated3 {
    /// Pair a planar curve with an elevation law.
    #[must_use]
    pub fn new(plan: Curve2, elevation: ElevationLaw) -> Self {
        Self {
            plan: Box::new(plan),
            elevation,
        }
    }
}
