#![forbid(unsafe_code)]

//! Filtered exact arithmetic (ADR 0068).
//!
//! Exact *constructions* -- where two segments cross, where a line meets a
//! circle -- produce numbers `f64` cannot hold. This crate answers each sign
//! question in at most two passes over the same expression:
//!
//! 1. [`Interval`] arithmetic: `f64` with every operation rounded outward,
//!    so the true value is provably inside `[lo, hi]`. If zero is outside,
//!    the sign is proven. This decides almost every real input.
//! 2. Otherwise [`Dyadic`] arithmetic: a big-integer mantissa (`num-bigint`)
//!    times a power of two. Every finite `f64` is one, and `+ - *` stay
//!    exact, so the sign is exact.
//!
//! There is no division, by design: constructions clear denominators, and a
//! quotient's sign is `sign(numerator) * sign(denominator)`. Square roots
//! appear only inside [`Root2`], the value `(a + b*sqrt(c)) / d`, whose sign
//! and order are decided by squaring with case analysis, never by
//! evaluating the root.
//!
//! Expressions are written once against the [`Arith`] trait and run in both
//! tiers, so the fast path and the exact path cannot compute different
//! polynomials.

pub mod arith;
pub mod certify;
pub mod construct;
pub mod dyadic;
pub mod interval;
pub mod root;

pub use arith::Arith;
pub use certify::{certify, filter, require_finite, ExactError, SignExpr};
pub use construct::{
    compare_along, crossing_orientation, crossing_orientation_filter, line_circle_hits, Branch,
    Circle, HitCount, Line, LineHit,
};
pub use dyadic::Dyadic;
pub use interval::Interval;
pub use root::{sign_root, sign_two_roots, Root2};
