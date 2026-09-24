//! The height wave a plane cuts across a cylinder, in the cylinder's own
//! parameters.
//!
//! A cylinder is parameterised by angle `u` and height `v`. A plane
//! `z = c0 + c1 x + c2 y` meets it where `x` and `y` are themselves `r cos u`
//! and `r sin u` (or `a cos u` and `b sin u` on an elliptical cylinder), so
//! the cut, read in `(u, v)`, is
//!
//! ```text
//! v(u) = mean + cosine * cos(u) + sine * sin(u)
//! ```
//!
//! That curve is what a trimmed cylinder face needs as the pcurve of its
//! sloped rim. No other [`Curve2`](crate::Curve2) variant holds it exactly:
//! it is not a conic in `(u, v)`, and a B-spline only approximates it,
//! because `u` is an angle and the wave is transcendental in it.

use axiolid_core::Scalar;

/// The graph `v = mean + cosine cos(t) + sine sin(t)`, parameterised by `t`.
///
/// The point at parameter `t` is `(t, v(t))`: the parameter IS the first
/// coordinate. On a cylinder that coordinate is the angle, so a pcurve use
/// states the angle span it covers directly as its interval, and a span
/// longer than a full turn is representable (it simply wraps the surface).
///
/// The curve is defined for every finite `t` and repeats with period `2 pi`
/// in `v`; its conventional domain is one full turn.
///
/// Dirty imported data stays representable, as everywhere else in this crate:
/// validation rejects non-finite coefficients, construction does not.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Sinusoid2 {
    /// Mean height: the constant term.
    pub mean: Scalar,
    /// Amplitude of the `cos(t)` term.
    pub cosine: Scalar,
    /// Amplitude of the `sin(t)` term.
    pub sine: Scalar,
}

impl Sinusoid2 {
    /// Height at parameter `t`.
    #[must_use]
    pub fn height(&self, t: Scalar) -> Scalar {
        let (s, c) = t.sin_cos();
        self.mean + self.cosine * c + self.sine * s
    }

    /// Peak deviation from the mean, `sqrt(cosine^2 + sine^2)`.
    ///
    /// The wave spans `mean - amplitude ..= mean + amplitude` over a full
    /// turn, which is what a caller needs to check a cut stays clear of
    /// another boundary.
    #[must_use]
    pub fn amplitude(&self) -> Scalar {
        self.cosine.hypot(self.sine)
    }

    /// Whether every coefficient is finite.
    #[must_use]
    pub fn is_finite(&self) -> bool {
        self.mean.is_finite() && self.cosine.is_finite() && self.sine.is_finite()
    }
}
