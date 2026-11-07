//! Outward-rounded interval arithmetic for the certified chord bounds (#414).
//!
//! Each operation widens its result by one ulp each way, and the
//! trigonometric ranges by [`TRIG_PAD`], so an enclosure computed in floating
//! point still holds the exact value.

use axiolid_core::Scalar;

use super::{contains_lattice, cos_range};

/// Absolute padding of a sine or cosine range, covering the library's
/// rounding of the endpoint values.
pub(super) const TRIG_PAD: Scalar = 8.0 * Scalar::EPSILON;

/// A closed interval whose operations round outward.
#[derive(Debug, Clone, Copy)]
pub(super) struct Iv {
    pub(super) lo: Scalar,
    pub(super) hi: Scalar,
}

impl Iv {
    pub(super) fn point(x: Scalar) -> Self {
        Self { lo: x, hi: x }
    }

    pub(super) fn widen(self) -> Self {
        Self {
            lo: self.lo.next_down(),
            hi: self.hi.next_up(),
        }
    }

    pub(super) fn add(self, o: Self) -> Self {
        Self {
            lo: self.lo + o.lo,
            hi: self.hi + o.hi,
        }
        .widen()
    }

    pub(super) fn neg(self) -> Self {
        Self {
            lo: -self.hi,
            hi: -self.lo,
        }
    }

    pub(super) fn mul(self, o: Self) -> Self {
        let p = [
            self.lo * o.lo,
            self.lo * o.hi,
            self.hi * o.lo,
            self.hi * o.hi,
        ];
        Self {
            lo: p.iter().copied().fold(Scalar::INFINITY, Scalar::min),
            hi: p.iter().copied().fold(Scalar::NEG_INFINITY, Scalar::max),
        }
        .widen()
    }

    pub(super) fn scale(self, k: Scalar) -> Self {
        self.mul(Self::point(k))
    }

    pub(super) fn sq(self) -> Self {
        let (a, b) = (self.lo * self.lo, self.hi * self.hi);
        let lo = if self.lo <= 0.0 && self.hi >= 0.0 {
            0.0
        } else {
            a.min(b).next_down().max(0.0)
        };
        Self {
            lo,
            hi: a.max(b).next_up(),
        }
    }

    /// `self / o` for an `o` holding no zero, `NaN` ends otherwise.
    pub(super) fn div(self, o: Self) -> Self {
        if o.straddles_zero() {
            return Self {
                lo: Scalar::NAN,
                hi: Scalar::NAN,
            };
        }
        let p = [
            self.lo / o.lo,
            self.lo / o.hi,
            self.hi / o.lo,
            self.hi / o.hi,
        ];
        Self {
            lo: p.iter().copied().fold(Scalar::INFINITY, Scalar::min),
            hi: p.iter().copied().fold(Scalar::NEG_INFINITY, Scalar::max),
        }
        .widen()
    }

    /// The square root of the non-negative part.
    pub(super) fn sqrt(self) -> Self {
        Self {
            lo: self.lo.max(0.0).sqrt().next_down().max(0.0),
            hi: self.hi.max(0.0).sqrt().next_up(),
        }
    }

    pub(super) fn straddles_zero(self) -> bool {
        self.lo <= 0.0 && self.hi >= 0.0
    }

    /// `sup |x|`.
    pub(super) fn mag(self) -> Scalar {
        self.lo.abs().max(self.hi.abs())
    }

    /// `inf |x|`, for an interval holding no zero.
    pub(super) fn mig(self) -> Scalar {
        self.lo.abs().min(self.hi.abs())
    }

    pub(super) fn finite(self) -> bool {
        self.lo.is_finite() && self.hi.is_finite()
    }
}

/// The range of `cos` over `[lo, hi]`.
pub(super) fn cos_iv(lo: Scalar, hi: Scalar) -> Iv {
    let (min, max) = cos_range(lo, hi);
    Iv {
        lo: min - TRIG_PAD,
        hi: max + TRIG_PAD,
    }
}

/// The range of `sin` over `[lo, hi]`.
pub(super) fn sin_iv(lo: Scalar, hi: Scalar) -> Iv {
    use core::f64::consts::{FRAC_PI_2, TAU};
    let max = if contains_lattice(lo, hi, FRAC_PI_2, TAU) {
        1.0
    } else {
        lo.sin().max(hi.sin())
    };
    let min = if contains_lattice(lo, hi, -FRAC_PI_2, TAU) {
        -1.0
    } else {
        lo.sin().min(hi.sin())
    };
    Iv {
        lo: min - TRIG_PAD,
        hi: max + TRIG_PAD,
    }
}
