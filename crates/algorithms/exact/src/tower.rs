//! Nested square roots: values in a tower of adjoined radicals.
//!
//! [`crate::Root2`] holds one square root over plain numbers. Arc
//! constructions need more: a hit on a circle is `a + b*sqrt(D)`, and a
//! distance or a second construction from that hit takes another square
//! root of an expression that already contains `sqrt(D)`. This module
//! represents such values exactly.
//!
//! # Representation
//!
//! A [`Tower`] is a list of radicands `r_1, ..., r_k`, each built only
//! from the radicals before it. A [`Nested`] value at level `k` is
//! `a + b*sqrt(r_k)` with `a` and `b` at level `k - 1`, stored flat as a
//! coefficient vector of length `2^k` (low half `a`, high half `b`).
//!
//! # Why no canonical form is needed
//!
//! Radicals need not be independent: `sqrt(8)` and `sqrt(2)` may both be
//! adjoined, and then `sqrt(8) - 2*sqrt(2)` has non-zero coefficients but
//! value zero. That is fine. Both the product rule
//! `(a + b√r)(c + d√r) = (ac + bd·r) + (ad + bc)√r` and the sign rule below
//! are identities about real numbers, true whatever the coefficients are.
//! So signs are exact, including exact zeros, without ever reducing to a
//! basis, which is what keeps this module small.
//!
//! # Sign
//!
//! `sign(a + b*sqrt(r))`, recursively by level: if `a` and `b*sqrt(r)`
//! agree in sign (or one is zero) that is the answer; otherwise it is
//! `sign(a) * sign(a^2 - b^2*r)`, a value one level down. Level 0 asks the
//! arithmetic `T` directly, so the same code runs as the interval filter
//! and as the exact fallback (see [`crate::certify()`]).
//!
//! # Cost
//!
//! Each level squares once, so the polynomial degree in the inputs, and
//! with it exact-tier mantissa length, doubles per level; a product costs
//! five products one level down. Depth is capped at [`MAX_DEPTH`].

use axiolid_guarantees::Sign;

use crate::arith::{sign_product, Arith};
use crate::certify::ExactError;

/// Most radicals one tower may hold. A product at depth `k` costs `5^k`
/// base products, and exact mantissas grow `2^k`-fold in degree.
pub const MAX_DEPTH: usize = 6;

/// A value in a [`Tower`]: `2^level` coefficients over `T`.
///
/// Only meaningful together with the tower that made it.
#[derive(Debug, Clone, PartialEq)]
pub struct Nested<T> {
    coeffs: Vec<T>,
}

impl<T> Nested<T> {
    /// How many radicals this value can involve.
    #[must_use]
    pub fn level(&self) -> usize {
        self.coeffs.len().trailing_zeros() as usize
    }

    /// The coefficients, low half first at every level.
    #[must_use]
    pub fn coeffs(&self) -> &[T] {
        &self.coeffs
    }
}

/// The radicals adjoined so far.
#[derive(Debug, Clone, Default)]
pub struct Tower<T> {
    /// Radicand `i` (0-based), lifted to exactly level `i`.
    radicands: Vec<Nested<T>>,
}

impl<T: Arith> Tower<T> {
    /// An empty tower: values are plain `T`.
    #[must_use]
    pub fn new() -> Self {
        Self {
            radicands: Vec::new(),
        }
    }

    /// Number of radicals adjoined.
    #[must_use]
    pub fn depth(&self) -> usize {
        self.radicands.len()
    }

    /// A plain value.
    #[must_use]
    pub fn value(&self, value: T) -> Nested<T> {
        Nested {
            coeffs: vec![value],
        }
    }

    /// The value of a finite `f64` (see [`Arith::from_f64`]).
    #[must_use]
    pub fn from_f64(&self, value: f64) -> Nested<T> {
        self.value(T::from_f64(value))
    }

    /// Adjoin `sqrt(radicand)` and return it as a value.
    ///
    /// The radicand must not be negative. That is not checked here, since
    /// the interval tier may be unable to tell; [`Tower::sign`] of any
    /// value involving a negative radicand returns `None`, which the exact
    /// tier reports as [`ExactError::Undefined`].
    ///
    /// # Errors
    ///
    /// [`ExactError::TooDeep`] beyond [`MAX_DEPTH`] radicals.
    ///
    /// # Panics
    ///
    /// If `radicand` came from a different tower (its level exceeds this
    /// tower's depth).
    pub fn sqrt(&mut self, radicand: &Nested<T>) -> Result<Nested<T>, ExactError> {
        let level = self.depth();
        if level >= MAX_DEPTH {
            return Err(ExactError::TooDeep);
        }
        self.own(radicand);
        self.radicands.push(self.lift(radicand, level));
        let mut coeffs = vec![T::from_f64(0.0); 1 << (level + 1)];
        coeffs[1 << level] = T::from_f64(1.0);
        Ok(Nested { coeffs })
    }

    /// `x + y`.
    #[must_use]
    pub fn add(&self, x: &Nested<T>, y: &Nested<T>) -> Nested<T> {
        self.zip(x, y, T::add)
    }

    /// `x - y`.
    #[must_use]
    pub fn sub(&self, x: &Nested<T>, y: &Nested<T>) -> Nested<T> {
        self.zip(x, y, T::sub)
    }

    /// `-x`.
    #[must_use]
    pub fn neg(&self, x: &Nested<T>) -> Nested<T> {
        self.own(x);
        Nested {
            coeffs: x.coeffs.iter().map(T::neg).collect(),
        }
    }

    /// `x * y`.
    #[must_use]
    pub fn mul(&self, x: &Nested<T>, y: &Nested<T>) -> Nested<T> {
        self.own(x);
        self.own(y);
        let level = x.level().max(y.level());
        let (x, y) = (self.lift(x, level), self.lift(y, level));
        Nested {
            coeffs: self.mul_at(level, &x.coeffs, &y.coeffs),
        }
    }

    /// The sign of `x`, or `None` when `T` cannot decide or a radicand
    /// involved is negative.
    ///
    /// Approximate arithmetics first evaluate `x` numerically through
    /// [`Arith::sqrt_enclosure`]; that decides whenever `x` is visibly away
    /// from zero. Otherwise, and always for exact arithmetics, signs are
    /// decided by case analysis on the coefficients.
    #[must_use]
    pub fn sign(&self, x: &Nested<T>) -> Option<Sign> {
        self.own(x);
        // An interval sign is sound whenever it is decided. A numeric
        // Zero can only come from [0, 0], i.e. from coefficients that are
        // literally zero, where Zero is also the exact answer.
        if let Some(sign) = self
            .numeric(x.level(), &x.coeffs)
            .and_then(|value| value.sign())
        {
            return Some(sign);
        }
        self.sign_at(x.level(), &x.coeffs)
    }

    /// Numerical value of a coefficient vector, for arithmetics that can
    /// enclose square roots. `None` for exact arithmetics, or when some
    /// radicand might be negative.
    fn numeric(&self, level: usize, x: &[T]) -> Option<T> {
        if level == 0 {
            return Some(x[0].clone());
        }
        let half = 1 << (level - 1);
        let (a, b) = x.split_at(half);
        let root = self
            .numeric(level - 1, &self.radicands[level - 1].coeffs)?
            .sqrt_enclosure()?;
        let a = self.numeric(level - 1, a)?;
        let b = self.numeric(level - 1, b)?;
        Some(a.add(&b.mul(&root)))
    }

    /// The sign of `x - y`.
    #[must_use]
    pub fn cmp(&self, x: &Nested<T>, y: &Nested<T>) -> Option<Sign> {
        self.sign(&self.sub(x, y))
    }

    fn own(&self, x: &Nested<T>) {
        assert!(
            x.level() <= self.depth(),
            "a Nested value is only meaningful in the tower that made it"
        );
    }

    /// `x` padded with zero coefficients to `level`.
    fn lift(&self, x: &Nested<T>, level: usize) -> Nested<T> {
        let mut coeffs = x.coeffs.clone();
        coeffs.resize(1 << level, T::from_f64(0.0));
        Nested { coeffs }
    }

    fn zip(&self, x: &Nested<T>, y: &Nested<T>, op: impl Fn(&T, &T) -> T) -> Nested<T> {
        self.own(x);
        self.own(y);
        let level = x.level().max(y.level());
        let (x, y) = (self.lift(x, level), self.lift(y, level));
        Nested {
            coeffs: x
                .coeffs
                .iter()
                .zip(&y.coeffs)
                .map(|(a, b)| op(a, b))
                .collect(),
        }
    }

    /// Product of two coefficient vectors of length `2^level`.
    fn mul_at(&self, level: usize, x: &[T], y: &[T]) -> Vec<T> {
        if level == 0 {
            return vec![x[0].mul(&y[0])];
        }
        let half = 1 << (level - 1);
        let (a, b) = x.split_at(half);
        let (c, d) = y.split_at(half);
        let radicand = &self.radicands[level - 1].coeffs;
        let ac = self.mul_at(level - 1, a, c);
        let bd = self.mul_at(level - 1, b, d);
        let bdr = self.mul_at(level - 1, &bd, radicand);
        let ad = self.mul_at(level - 1, a, d);
        let bc = self.mul_at(level - 1, b, c);
        let mut out: Vec<T> = ac.iter().zip(&bdr).map(|(p, q)| p.add(q)).collect();
        out.extend(ad.iter().zip(&bc).map(|(p, q)| p.add(q)));
        out
    }

    fn sign_at(&self, level: usize, x: &[T]) -> Option<Sign> {
        if level == 0 {
            return x[0].sign();
        }
        let half = 1 << (level - 1);
        let (a, b) = x.split_at(half);
        let radicand = &self.radicands[level - 1].coeffs;
        let sr = self.sign_at(level - 1, radicand)?;
        if sr == Sign::Negative {
            return None;
        }
        let sa = self.sign_at(level - 1, a)?;
        // b*sqrt(r) has the sign of b, unless r is zero.
        let sb = if sr == Sign::Zero {
            Sign::Zero
        } else {
            self.sign_at(level - 1, b)?
        };
        if sb == Sign::Zero {
            return Some(sa);
        }
        if sa == Sign::Zero || sa == sb {
            return Some(sb);
        }
        let a2 = self.mul_at(level - 1, a, a);
        let b2 = self.mul_at(level - 1, b, b);
        let b2r = self.mul_at(level - 1, &b2, radicand);
        let diff: Vec<T> = a2.iter().zip(&b2r).map(|(p, q)| p.sub(q)).collect();
        let dominance = self.sign_at(level - 1, &diff)?;
        Some(sign_product(sa, dominance))
    }
}
