//! Intrinsic (natural) plane curve data: curvature as a function of arc length.
//!
//! A plane curve is fixed up to rigid motion by its curvature law k(s). Anchoring
//! that law to a start frame fixes it absolutely. This is the *natural equation*
//! (Cesaro/Whewell) rather than a parametric map, and it is the honest way to
//! carry a clothoid: the spiral has no elementary parametric form, so a kernel
//! that only speaks parametrically has to approximate it before it has even
//! stored it.
//!
//! # Exact and symbolic, deliberately not evaluated
//!
//! Everything here is closed-form symbolic data and closed-form symbolic
//! operations on it: differentiate the law, negate it, measure the total turning
//! it accumulates. No module in this crate turns a law into points.
//!
//! That boundary is not laziness, it is the mathematics. Recovering position from
//! k(s) requires integrating the tangent angle and then integrating the tangent,
//! which for a linear law is the Fresnel integral -- not an elementary function.
//! Any point you get is a quadrature result carrying a tolerance, so it belongs
//! to an evaluator that can state that tolerance, never to the representation.
//! Storing the law exactly and refusing to fake points is what keeps a clothoid
//! a clothoid instead of a polyline that used to be one.

use axiolid_core::{Frame2, Scalar};

/// Curvature as a closed-form function of arc length.
///
/// Sign follows the usual plane convention: positive curvature turns the
/// tangent counter-clockwise. `s` is arc length measured from the curve start,
/// so a law is meaningful on `[0, length]` of the curve that carries it.
#[non_exhaustive]
#[derive(Debug, Clone, PartialEq)]
pub enum CurvatureLaw {
    /// `k(s) = curvature`.
    ///
    /// Zero is a straight line; any other value is a circular arc of radius
    /// `1 / curvature`.
    Constant {
        /// The constant curvature.
        curvature: Scalar,
    },
    /// `k(s) = coefficients[0] + coefficients[1] * s + coefficients[2] * s^2 + ...`
    ///
    /// Degree 1 is the clothoid (Euler spiral), whose curvature is linear in
    /// arc length; that is the transition spiral used between straight and
    /// circular track so lateral acceleration ramps linearly instead of
    /// stepping. Higher degrees cover the Bloss and cubic-parabola families.
    ///
    /// An empty coefficient list is the zero polynomial: a straight line.
    Polynomial {
        /// Coefficients in ascending powers of arc length.
        coefficients: Vec<Scalar>,
    },
    /// `k(s) = mean + amplitude * sin(angular_frequency * s + phase)`.
    ///
    /// The sinusoidal transition family (Klein, cosine ramps). Kept distinct
    /// from `Polynomial` because it is exactly representable this way and a
    /// truncated series would not be.
    Sinusoid {
        /// Curvature the oscillation is centred on.
        mean: Scalar,
        /// Peak deviation from `mean`.
        amplitude: Scalar,
        /// Radians of phase per unit arc length.
        angular_frequency: Scalar,
        /// Phase offset at `s = 0`, in radians.
        phase: Scalar,
    },
}

impl CurvatureLaw {
    /// The zero law: a straight line.
    #[must_use]
    pub const fn straight() -> Self {
        Self::Constant { curvature: 0.0 }
    }

    /// A circular arc of the given signed curvature.
    #[must_use]
    pub const fn circular(curvature: Scalar) -> Self {
        Self::Constant { curvature }
    }

    /// A clothoid whose curvature runs from `start` to `end` over `length`.
    ///
    /// This is the transition-spiral constructor: the rate is derived rather
    /// than asked for, because the two endpoint curvatures and the length are
    /// what an alignment actually specifies.
    ///
    /// A non-finite or zero `length` cannot define a rate, so the result is a
    /// constant `start` law -- the honest degenerate answer, not a division by
    /// zero smuggled into a coefficient.
    #[must_use]
    pub fn clothoid(start: Scalar, end: Scalar, length: Scalar) -> Self {
        if !length.is_finite() || length == 0.0 {
            return Self::Constant { curvature: start };
        }
        Self::Polynomial {
            coefficients: vec![start, (end - start) / length],
        }
    }

    /// Whether the law is identically zero, i.e. a straight line.
    ///
    /// Exact: this is a structural test on the stored coefficients, not a
    /// sampled one.
    #[must_use]
    pub fn is_straight(&self) -> bool {
        match self {
            Self::Constant { curvature } => *curvature == 0.0,
            Self::Polynomial { coefficients } => coefficients.iter().all(|c| *c == 0.0),
            Self::Sinusoid {
                mean,
                amplitude,
                angular_frequency,
                ..
            } => {
                // A zero frequency freezes the sine at its phase value, so the
                // law is constant but not necessarily zero; that case is only
                // straight when the frozen value cancels the mean, which needs
                // evaluation. Report false rather than guess.
                *mean == 0.0 && *amplitude == 0.0 && angular_frequency.is_finite()
            }
        }
    }

    /// Whether the law is constant in arc length.
    #[must_use]
    pub fn is_constant(&self) -> bool {
        match self {
            Self::Constant { .. } => true,
            Self::Polynomial { coefficients } => coefficients.iter().skip(1).all(|c| *c == 0.0),
            Self::Sinusoid { amplitude, .. } => *amplitude == 0.0,
        }
    }

    /// The derivative law `dk/ds`, in closed form.
    ///
    /// Exact symbolic differentiation. The sharpness of a clothoid is the
    /// constant this returns.
    #[must_use]
    pub fn derivative(&self) -> Self {
        match self {
            Self::Constant { .. } => Self::Constant { curvature: 0.0 },
            Self::Polynomial { coefficients } => Self::Polynomial {
                coefficients: coefficients
                    .iter()
                    .enumerate()
                    .skip(1)
                    .map(|(power, c)| *c * power as Scalar)
                    .collect(),
            },
            // d/ds [m + A sin(w s + p)] = A w cos(w s + p)
            //                           = A w sin(w s + p + pi/2)
            // Folding the cosine back into a sine keeps the family closed under
            // differentiation, so no new variant is needed.
            Self::Sinusoid {
                amplitude,
                angular_frequency,
                phase,
                ..
            } => Self::Sinusoid {
                mean: 0.0,
                amplitude: amplitude * angular_frequency,
                angular_frequency: *angular_frequency,
                phase: phase + core::f64::consts::FRAC_PI_2,
            },
        }
    }

    /// The mirrored law `-k(s)`: the same curve reflected.
    #[must_use]
    pub fn reversed_orientation(&self) -> Self {
        match self {
            Self::Constant { curvature } => Self::Constant {
                curvature: -curvature,
            },
            Self::Polynomial { coefficients } => Self::Polynomial {
                coefficients: coefficients.iter().map(|c| -c).collect(),
            },
            Self::Sinusoid {
                mean,
                amplitude,
                angular_frequency,
                phase,
            } => Self::Sinusoid {
                mean: -mean,
                amplitude: -amplitude,
                angular_frequency: *angular_frequency,
                phase: *phase,
            },
        }
    }
}

/// A plane curve given by its natural equation: a curvature law anchored to a
/// start frame and run for a finite arc length.
///
/// The frame supplies the rigid motion that k(s) alone cannot: its origin is the
/// curve start, and its `x` axis is the start tangent direction.
///
/// Like every other value in this crate, dirty imported data stays
/// representable -- a non-positive length is storable, and naming it is the
/// job of a validator, not of the type.
#[derive(Debug, Clone, PartialEq)]
pub struct Intrinsic2 {
    /// Start frame: origin at the curve start, `x` along the start tangent.
    pub start: Frame2,
    /// Curvature as a function of arc length from `start`.
    pub curvature: CurvatureLaw,
    /// Arc length the law is defined over.
    pub length: Scalar,
}

impl Intrinsic2 {
    /// Anchor a curvature law to a start frame over an arc length.
    #[must_use]
    pub const fn new(start: Frame2, curvature: CurvatureLaw, length: Scalar) -> Self {
        Self {
            start,
            curvature,
            length,
        }
    }

    /// Whether the curve is a straight segment.
    #[must_use]
    pub fn is_straight(&self) -> bool {
        self.curvature.is_straight()
    }

    /// Total tangent turning over the curve, in radians, in closed form.
    ///
    /// This is the integral of k(s) over `[0, length]` -- exact for every law in
    /// the family, because each one has an elementary antiderivative. It is the
    /// *angle* that integrates in closed form; position does not, which is why
    /// this method exists and an `evaluate` does not.
    ///
    /// Returns `None` when the length is not finite, since the integral is then
    /// undefined rather than merely large.
    #[must_use]
    pub fn total_turning(&self) -> Option<Scalar> {
        if !self.length.is_finite() {
            return None;
        }
        let s = self.length;
        match &self.curvature {
            CurvatureLaw::Constant { curvature } => Some(curvature * s),
            // Integral of sum c_i s^i is sum c_i s^(i+1) / (i+1).
            CurvatureLaw::Polynomial { coefficients } => Some(
                coefficients
                    .iter()
                    .enumerate()
                    .map(|(power, c)| c * s.powi(power as i32 + 1) / (power as Scalar + 1.0))
                    .sum(),
            ),
            CurvatureLaw::Sinusoid {
                mean,
                amplitude,
                angular_frequency,
                phase,
            } => {
                // Integral of m + A sin(w s + p) is m s - (A/w)[cos(w s + p) - cos(p)].
                // A zero frequency makes that undefined, but the integrand is then
                // the constant m + A sin(p), which integrates directly.
                if *angular_frequency == 0.0 {
                    return Some((mean + amplitude * phase.sin()) * s);
                }
                let turn = mean * s
                    - (amplitude / angular_frequency)
                        * ((angular_frequency * s + phase).cos() - phase.cos());
                Some(turn)
            }
        }
    }
}
