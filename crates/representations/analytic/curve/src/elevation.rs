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

use axiolid_core::Scalar;

use crate::Curve2;

/// Height as a function of distance along the plan.
///
/// Mirrors [`CurvatureLaw`](crate::CurvatureLaw) in shape so the two halves of
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
    /// matches [`CurvatureLaw::Piecewise`](crate::CurvatureLaw::Piecewise)
    /// deliberately: a vertical profile is authored as a run of segments and
    /// the seams are observable data.
    Piecewise {
        /// Interior seam positions in plan distance, ascending.
        breaks: Vec<Scalar>,
        /// One law per piece; `laws.len() == breaks.len() + 1`.
        laws: Vec<ElevationLaw>,
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

    /// Whether the piece lists agree and the seams ascend.
    #[must_use]
    pub fn is_well_formed(&self) -> bool {
        match self {
            Self::Polynomial { coefficients } => coefficients.iter().all(|c| c.is_finite()),
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
        }
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
