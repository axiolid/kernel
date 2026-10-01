//! Real polynomials in the power basis and all their real roots in an
//! interval.
//!
//! # Method
//!
//! Roots of `p` are isolated by the roots of `p'`, recursively (the
//! approach of OCCT `math_DirectPolynomialRoots` generalised, and of
//! Collins and Loos' derivative sequence): between consecutive critical
//! points `p` is monotone, so a sign change there is exactly one root,
//! found by [`find_root`].
//!
//! Every evaluation carries a rigorous bound on its own rounding error
//! (Horner with `gamma(2n) * sum |a_i| |x|^i`, Higham, *Accuracy and
//! Stability of Numerical Algorithms*, 5.1). A sign is trusted only where
//! `|p(x)|` exceeds that bound. Where it does not -- at a multiple root, a
//! cluster closer than `f64` can separate, or a near-touch of the axis --
//! the answer is [`RootKind::Unresolved`] with the region and the most
//! roots it can hold, rather than a guessed count.
//!
//! The roots are those of the polynomial whose coefficients are the given
//! `f64` values. Monotonicity between critical points holds to working
//! precision: the critical points come from the rounded derivative.

use crate::error::{all_finite, finite, NumericError, NumericResult};
use crate::root::find_root;

/// A real polynomial `a[0] + a[1] x + ... + a[n] x^n` with finite
/// coefficients.
#[derive(Debug, Clone, PartialEq)]
pub struct Polynomial {
    coefficients: Vec<f64>,
}

/// How certain a reported root is.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum RootKind {
    /// The polynomial is monotone on `[lower, upper]` and its sign
    /// certifiably changes across it: exactly one simple root lies inside.
    Simple,
    /// `|p|` is within its rounding bound near `x`. Between zero and
    /// `max_count` roots, counted with multiplicity, lie in
    /// `[lower, upper]`: a multiple root, a cluster, or a near-touch that
    /// `f64` evaluation cannot tell apart. A root at an end of the
    /// searched interval is reported this way, with `max_count` 1, since
    /// it may lie just outside.
    Unresolved {
        /// Most roots the region can hold (by Rolle's theorem: one more
        /// than the critical points inside it).
        max_count: usize,
    },
}

/// A real root, or an unresolved root region, of a polynomial.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct PolynomialRoot {
    /// Best estimate.
    pub x: f64,
    /// Lower end of the enclosure.
    pub lower: f64,
    /// Upper end of the enclosure. `|x - root| <= upper - lower` for every
    /// root the entry accounts for.
    pub upper: f64,
    /// Whether the entry is one certified simple root.
    pub kind: RootKind,
}

impl PolynomialRoot {
    /// A bound on the distance from `x` to each root this entry accounts
    /// for: the enclosure width.
    pub fn error_bound(&self) -> f64 {
        self.upper - self.lower
    }
}

impl Polynomial {
    /// A polynomial from coefficients in ascending order of power.
    /// Trailing zero coefficients are dropped.
    ///
    /// # Errors
    ///
    /// Refuses a non-finite coefficient.
    pub fn new(coefficients: &[f64]) -> NumericResult<Self> {
        all_finite(coefficients, "polynomial coefficient")?;
        let mut coefficients = coefficients.to_vec();
        while coefficients.last() == Some(&0.0) {
            coefficients.pop();
        }
        Ok(Self { coefficients })
    }

    /// `(x - r_0)(x - r_1)...`, expanded. A convenience for building test
    /// and reference polynomials.
    ///
    /// # Errors
    ///
    /// Refuses a non-finite root or an expansion that overflows.
    pub fn from_roots(roots: &[f64]) -> NumericResult<Self> {
        all_finite(roots, "root")?;
        let mut c = vec![1.0];
        for &r in roots {
            let mut next = vec![0.0; c.len() + 1];
            for (i, &ci) in c.iter().enumerate() {
                next[i + 1] += ci;
                next[i] -= r * ci;
            }
            c = next;
        }
        if c.iter().any(|v| !v.is_finite()) {
            return Err(NumericError::InvalidArgument {
                name: "roots",
                reason: "expand to coefficients that overflow",
            });
        }
        Self::new(&c)
    }

    /// Coefficients in ascending order of power, without trailing zeros.
    pub fn coefficients(&self) -> &[f64] {
        &self.coefficients
    }

    /// Degree, or `None` for the zero polynomial.
    pub fn degree(&self) -> Option<usize> {
        self.coefficients.len().checked_sub(1)
    }

    /// `p(x)` by Horner's rule.
    pub fn evaluate(&self, x: f64) -> f64 {
        self.coefficients
            .iter()
            .rev()
            .fold(0.0, |acc, &a| acc * x + a)
    }

    /// `p(x)` and a rigorous bound on the rounding error of that value.
    ///
    /// The bound is `gamma(2n + 1) * sum |a_i| |x|^i` plus an underflow
    /// allowance, where `gamma(k) = k u / (1 - k u)`.
    pub fn evaluate_with_bound(&self, x: f64) -> (f64, f64) {
        let n = self.coefficients.len();
        let mut value = 0.0f64;
        let mut magnitude = 0.0f64;
        let ax = x.abs();
        for &a in self.coefficients.iter().rev() {
            value = value * x + a;
            magnitude = magnitude * ax + a.abs();
        }
        let u = f64::EPSILON / 2.0;
        let k = (2 * n + 1) as f64;
        let gamma = k * u / (1.0 - k * u);
        // Products may underflow: allow one smallest normal per operation.
        let bound = 2.0 * gamma * magnitude + (2 * n) as f64 * f64::MIN_POSITIVE;
        (value, bound)
    }

    /// The derivative.
    pub fn derivative(&self) -> Polynomial {
        let coefficients: Vec<f64> = self
            .coefficients
            .iter()
            .enumerate()
            .skip(1)
            .map(|(i, &a)| a * i as f64)
            .collect();
        Polynomial { coefficients }
    }

    /// Every real root in `[lower, upper]`, in increasing order.
    ///
    /// Each entry is either one certified [`RootKind::Simple`] root with an
    /// enclosure, or an [`RootKind::Unresolved`] region with the most
    /// roots it can hold. No root of the polynomial in the interval is left
    /// out: a root not reported as simple lies in an unresolved region.
    ///
    /// # Errors
    ///
    /// Refuses non-finite bounds, `lower > upper`, and the zero polynomial
    /// ([`NumericError::ZeroPolynomial`]).
    pub fn real_roots(&self, lower: f64, upper: f64) -> NumericResult<Vec<PolynomialRoot>> {
        finite(lower, "lower bound")?;
        finite(upper, "upper bound")?;
        if lower > upper {
            return Err(NumericError::InvalidArgument {
                name: "interval",
                reason: "must have lower <= upper",
            });
        }
        if self.coefficients.is_empty() {
            return Err(NumericError::ZeroPolynomial);
        }
        self.isolate(lower, upper)
    }

    fn sign(&self, x: f64) -> i8 {
        let (v, bound) = self.evaluate_with_bound(x);
        if v.abs() <= bound {
            0
        } else if v > 0.0 {
            1
        } else {
            -1
        }
    }

    fn isolate(&self, lower: f64, upper: f64) -> NumericResult<Vec<PolynomialRoot>> {
        if self.coefficients.len() <= 1 {
            // A non-zero constant has no roots.
            return Ok(Vec::new());
        }
        let critical = self.derivative().isolate(lower, upper)?;

        // Breakpoints: the interval ends and each critical point, with the
        // number of critical points (with multiplicity) each carries.
        let mut points: Vec<(f64, usize)> = vec![(lower, 0)];
        for c in &critical {
            let count = match c.kind {
                RootKind::Simple => 1,
                RootKind::Unresolved { max_count } => max_count,
            };
            let x = c.x.clamp(lower, upper);
            match points.last_mut() {
                Some(last) if last.0 == x => last.1 += count,
                _ => points.push((x, count)),
            }
        }
        match points.last_mut() {
            Some(last) if last.0 == upper => {}
            _ => points.push((upper, 0)),
        }
        let signs: Vec<i8> = points.iter().map(|&(x, _)| self.sign(x)).collect();

        let mut roots = Vec::new();
        let mut i = 0;
        while i < points.len() {
            if signs[i] == 0 {
                // A run of breakpoints where p is indistinguishable from zero.
                let start = i;
                let mut crit_count = 0;
                while i < points.len() && signs[i] == 0 {
                    crit_count += points[i].1;
                    i += 1;
                }
                let end = i - 1;
                let left_limit = if start > 0 {
                    points[start - 1].0
                } else {
                    lower
                };
                let right_limit = if i < points.len() { points[i].0 } else { upper };
                let x = if start == end {
                    points[start].0
                } else {
                    0.5 * points[start].0 + 0.5 * points[end].0
                };
                roots.push(PolynomialRoot {
                    x,
                    lower: self.boundary(points[start].0, left_limit),
                    upper: self.boundary(points[end].0, right_limit),
                    kind: RootKind::Unresolved {
                        max_count: crit_count + 1,
                    },
                });
                continue;
            }
            if i + 1 < points.len() && signs[i + 1] != 0 && signs[i] != signs[i + 1] {
                let (a, b) = (points[i].0, points[i + 1].0);
                let root = find_root(|x| self.evaluate(x), a, b, 0.0)?;
                let lo = if self.sign(root.lower) != 0 {
                    root.lower
                } else {
                    self.boundary(root.lower, a)
                };
                let hi = if self.sign(root.upper) != 0 {
                    root.upper
                } else {
                    self.boundary(root.upper, b)
                };
                roots.push(PolynomialRoot {
                    x: root.x,
                    lower: lo,
                    upper: hi,
                    kind: RootKind::Simple,
                });
            }
            i += 1;
        }
        Ok(roots)
    }

    /// From `inside`, where the sign is uncertain, move toward `limit` to
    /// the nearest point found whose sign is certain. Returns `limit` if
    /// the sign is uncertain there too.
    fn boundary(&self, inside: f64, limit: f64) -> f64 {
        if inside == limit || self.sign(limit) == 0 {
            return limit;
        }
        let (mut uncertain, mut certain) = (inside, limit);
        for _ in 0..200 {
            let mid = 0.5 * uncertain + 0.5 * certain;
            if mid == uncertain || mid == certain {
                break;
            }
            if self.sign(mid) == 0 {
                uncertain = mid;
            } else {
                certain = mid;
            }
        }
        certain
    }
}
