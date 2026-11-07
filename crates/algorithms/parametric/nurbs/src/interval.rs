//! Closed intervals, rounded outward (ADR 0077).
//!
//! Every operation widens its result by a few ulps, so an enclosure
//! computed with [`I`] contains the exact value for every point of its
//! operands. The certificates of the section tracer and the root isolation
//! along a traced curve are built on it.

use axiolid_core::Scalar;

/// A closed interval, rounded outward.
#[derive(Debug, Clone, Copy, PartialEq)]
pub(crate) struct I {
    pub(crate) lo: Scalar,
    pub(crate) hi: Scalar,
}

fn down(x: Scalar) -> Scalar {
    x - x.abs() * 4.0 * Scalar::EPSILON - Scalar::MIN_POSITIVE
}

fn up(x: Scalar) -> Scalar {
    x + x.abs() * 4.0 * Scalar::EPSILON + Scalar::MIN_POSITIVE
}

impl I {
    pub(crate) fn new(a: Scalar, b: Scalar) -> Self {
        Self {
            lo: down(a.min(b)),
            hi: up(a.max(b)),
        }
    }

    pub(crate) fn point(x: Scalar) -> Self {
        Self { lo: x, hi: x }
    }

    /// `x` give or take `r`.
    pub(crate) fn around(x: Scalar, r: Scalar) -> Self {
        Self::new(x - r.abs(), x + r.abs())
    }

    pub(crate) fn add(self, o: Self) -> Self {
        Self::new(self.lo + o.lo, self.hi + o.hi)
    }

    pub(crate) fn sub(self, o: Self) -> Self {
        Self::new(self.lo - o.hi, self.hi - o.lo)
    }

    pub(crate) fn mul(self, o: Self) -> Self {
        let c = [
            self.lo * o.lo,
            self.lo * o.hi,
            self.hi * o.lo,
            self.hi * o.hi,
        ];
        Self::new(
            c.iter().copied().fold(Scalar::INFINITY, Scalar::min),
            c.iter().copied().fold(Scalar::NEG_INFINITY, Scalar::max),
        )
    }

    pub(crate) fn scale(self, k: Scalar) -> Self {
        Self::new(self.lo * k, self.hi * k)
    }

    /// `self / o` for `o` bounded away from zero.
    pub(crate) fn div(self, o: Self) -> Option<Self> {
        if o.lo <= 0.0 && o.hi >= 0.0 {
            return None;
        }
        let c = [
            self.lo / o.lo,
            self.lo / o.hi,
            self.hi / o.lo,
            self.hi / o.hi,
        ];
        Some(Self::new(
            c.iter().copied().fold(Scalar::INFINITY, Scalar::min),
            c.iter().copied().fold(Scalar::NEG_INFINITY, Scalar::max),
        ))
    }

    pub(crate) fn hull(self, o: Self) -> Self {
        Self {
            lo: self.lo.min(o.lo),
            hi: self.hi.max(o.hi),
        }
    }

    pub(crate) fn meets(self, o: Self) -> bool {
        self.lo <= o.hi && o.lo <= self.hi
    }

    /// Strictly inside `o`.
    pub(crate) fn within(self, o: Self) -> bool {
        self.lo > o.lo && self.hi < o.hi
    }

    pub(crate) fn mid(self) -> Scalar {
        0.5 * (self.lo + self.hi)
    }

    pub(crate) fn is_finite(self) -> bool {
        self.lo.is_finite() && self.hi.is_finite()
    }
}
