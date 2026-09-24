//! Soundness of the fast tier: the exact value always lies in the interval.
//!
//! Random expressions are evaluated in both tiers from the same inputs. The
//! dyadic result is exact, so "contains" is checked exactly, by comparing
//! dyadic values, never in floating point. Inputs span the whole `f64`
//! range: subnormals, huge magnitudes, exact zeros, mixed exponents.

use axiolid_exact::{Arith, Dyadic, Interval};
use axiolid_guarantees::Sign;

struct Rng(u64);

impl Rng {
    fn next(&mut self) -> u64 {
        self.0 ^= self.0 << 13;
        self.0 ^= self.0 >> 7;
        self.0 ^= self.0 << 17;
        self.0
    }

    fn below(&mut self, n: u64) -> u64 {
        self.next() % n
    }

    /// Finite doubles from every regime, not just "nice" ones.
    fn float(&mut self) -> f64 {
        match self.below(6) {
            0 => 0.0,
            1 => (self.below(2001) as f64) - 1000.0,
            // Any finite bit pattern: subnormals and huge values included.
            2 => loop {
                let value = f64::from_bits(self.next());
                if value.is_finite() {
                    break value;
                }
            },
            3 => f64::from_bits(self.below(1 << 20)) * if self.below(2) == 0 { 1.0 } else { -1.0 },
            _ => ((self.next() >> 11) as f64 / (1u64 << 53) as f64) * 2000.0 - 1000.0,
        }
    }
}

#[derive(Clone)]
enum Expr {
    Leaf(f64),
    Add(Box<Expr>, Box<Expr>),
    Sub(Box<Expr>, Box<Expr>),
    Mul(Box<Expr>, Box<Expr>),
    Neg(Box<Expr>),
    Square(Box<Expr>),
}

fn random_expr(rng: &mut Rng, depth: u32) -> Expr {
    if depth == 0 || rng.below(4) == 0 {
        return Expr::Leaf(rng.float());
    }
    let choice = rng.below(5);
    let mut child = || Box::new(random_expr(rng, depth - 1));
    match choice {
        0 => Expr::Add(child(), child()),
        1 => Expr::Sub(child(), child()),
        2 => Expr::Mul(child(), child()),
        3 => Expr::Neg(child()),
        _ => Expr::Square(child()),
    }
}

fn eval<T: Arith>(expr: &Expr) -> T {
    match expr {
        Expr::Leaf(value) => T::from_f64(*value),
        Expr::Add(l, r) => eval::<T>(l).add(&eval::<T>(r)),
        Expr::Sub(l, r) => eval::<T>(l).sub(&eval::<T>(r)),
        Expr::Mul(l, r) => eval::<T>(l).mul(&eval::<T>(r)),
        Expr::Neg(e) => eval::<T>(e).neg(),
        Expr::Square(e) => eval::<T>(e).square(),
    }
}

/// `bound <= exact` (lower = true) or `exact <= bound`, decided exactly.
/// Infinite bounds hold trivially in their own direction.
fn bound_holds(bound: f64, exact: &Dyadic, lower: bool) -> bool {
    if bound.is_nan() {
        return false;
    }
    if bound.is_infinite() {
        return (bound < 0.0) == lower;
    }
    let difference = Dyadic::from_f64(bound).sub(exact).sign();
    if lower {
        difference != Some(Sign::Positive)
    } else {
        difference != Some(Sign::Negative)
    }
}

#[test]
fn every_interval_contains_the_exact_value() {
    let mut rng = Rng(0x1234_5678_9ABC_DEF1);
    let (mut decided, mut total) = (0, 0);
    for _ in 0..20_000 {
        let expr = random_expr(&mut rng, 4);
        let interval: Interval = eval(&expr);
        let exact: Dyadic = eval(&expr);
        total += 1;
        assert!(
            bound_holds(interval.lo(), &exact, true) && bound_holds(interval.hi(), &exact, false),
            "interval [{}, {}] misses the exact value (mantissa bits {})",
            interval.lo(),
            interval.hi(),
            exact.bits()
        );
        // A decided interval sign must be the exact sign.
        if let Some(sign) = interval.sign() {
            decided += 1;
            assert_eq!(Some(sign), exact.sign(), "filter certified a wrong sign");
        }
    }
    // Guard against a vacuous pass: most random expressions must decide.
    assert!(
        decided * 2 > total,
        "filter decided only {decided} of {total}"
    );
}

#[test]
fn nan_and_infinity_never_decide_a_sign() {
    assert_eq!(Interval::point(f64::NAN).sign(), None);
    assert_eq!(Interval::point(f64::INFINITY).sign(), None);
    // inf - inf is NaN inside the arithmetic; it must widen, not decide.
    let huge = Interval::point(f64::MAX);
    let overflowed = huge.mul(&huge);
    assert!(
        overflowed.lo() > 0.0,
        "an overflowed positive product stays positive"
    );
    assert_eq!(overflowed.sub(&overflowed).sign(), None);
}

#[test]
fn the_filter_proves_nothing_about_a_cancelling_sum() {
    // 0.1 + 0.2 - 0.3 is not exactly zero in doubles (it is 2^-54), and
    // the interval must not claim zero: it widens past it.
    let value = Interval::point(0.1)
        .add(&Interval::point(0.2))
        .sub(&Interval::point(0.3));
    let exact = Dyadic::from_f64(0.1)
        .add(&Dyadic::from_f64(0.2))
        .sub(&Dyadic::from_f64(0.3));
    assert_eq!(exact.sign(), Some(Sign::Positive));
    assert!(value.sign() != Some(Sign::Zero));
    assert!(value.sign().is_none() || value.sign() == exact.sign());
}

#[test]
fn square_root_enclosures_contain_the_true_root() {
    // lo <= sqrt(x) <= hi  <=>  lo^2 <= x <= hi^2 for non-negative bounds,
    // and squares of doubles are exact in dyadic arithmetic. So containment
    // is checked exactly, with no rounded root anywhere in the test.
    let mut rng = Rng(0x0DDB_A11C_0FFE_E123);
    let mut checked = 0;
    for _ in 0..20_000 {
        let value = rng.float().abs();
        let Some(root) = Interval::point(value).sqrt_enclosure() else {
            panic!("a non-negative point must have a root enclosure");
        };
        let exact = Dyadic::from_f64(value);
        let lo = Dyadic::from_f64(root.lo());
        let hi = Dyadic::from_f64(root.hi());
        assert!(root.lo() >= 0.0, "a root enclosure never goes below zero");
        assert_ne!(
            lo.square().sub(&exact).sign(),
            Some(Sign::Positive),
            "lo^2 <= x for {value:e}"
        );
        assert_ne!(
            hi.square().sub(&exact).sign(),
            Some(Sign::Negative),
            "x <= hi^2 for {value:e}"
        );
        checked += 1;
    }
    assert_eq!(checked, 20_000);
    // A radicand that might be negative has no enclosure: the filter must
    // not give a sign to a value that may not be real.
    let straddling = Interval::point(1.0).sub(&Interval::point(1.0));
    assert!(straddling.sqrt_enclosure().is_none());
    assert!(Interval::point(-4.0).sqrt_enclosure().is_none());
    assert!(Interval::WHOLE.sqrt_enclosure().is_none());
    // Exact tiers decide by case analysis, never by a rounded root.
    assert!(Dyadic::from_f64(4.0).sqrt_enclosure().is_none());
}
