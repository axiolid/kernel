//! `Dyadic::enclosure` is sound: the exact value lies in the interval.
//!
//! Checked exactly (dyadic comparisons, never floats) over products and
//! sums that produce mantissas far longer than 64 bits, where the floor
//! shift and the final rounding both matter.

use axiolid_exact::{Arith, Dyadic};
use axiolid_guarantees::Sign;

struct Rng(u64);

impl Rng {
    fn next(&mut self) -> u64 {
        self.0 ^= self.0 << 13;
        self.0 ^= self.0 >> 7;
        self.0 ^= self.0 << 17;
        self.0
    }

    /// Finite doubles over a wide but not overflowing range.
    fn float(&mut self) -> f64 {
        let mantissa = (self.next() >> 11) as f64 / (1u64 << 53) as f64;
        let exponent = (self.next() % 120) as i32 - 60;
        let sign = if self.next() % 2 == 0 { 1.0 } else { -1.0 };
        sign * (0.5 + mantissa) * 2f64.powi(exponent)
    }
}

fn below_or_equal(bound: f64, value: &Dyadic) -> bool {
    Dyadic::from_f64(bound).sub(value).sign() != Some(Sign::Positive)
}

fn above_or_equal(bound: f64, value: &Dyadic) -> bool {
    Dyadic::from_f64(bound).sub(value).sign() != Some(Sign::Negative)
}

#[test]
fn every_enclosure_contains_the_exact_value() {
    let mut rng = Rng(0xE7C1_05E5_0000_0001);
    let mut long = 0;
    for _ in 0..20_000 {
        let (a, b, c, e) = (rng.float(), rng.float(), rng.float(), rng.float());
        // a*b + c*e has a mantissa of up to ~120 bits after alignment.
        let value = Dyadic::from_f64(a)
            .mul(&Dyadic::from_f64(b))
            .add(&Dyadic::from_f64(c).mul(&Dyadic::from_f64(e)));
        long += usize::from(value.bits() > 64);
        let enclosure = value.enclosure();
        assert!(
            below_or_equal(enclosure.lo(), &value),
            "lo {} above value",
            enclosure.lo()
        );
        assert!(
            above_or_equal(enclosure.hi(), &value),
            "hi {} below value",
            enclosure.hi()
        );
        // Tight enough to decide: a few ulps, not the whole line.
        assert!(enclosure.lo().is_finite() && enclosure.hi().is_finite());
    }
    // Guard against a vacuous pass: most cases exercise the floor shift.
    assert!(long > 10_000, "only {long} long mantissas");
}

#[test]
fn exact_small_values_enclose_themselves_tightly() {
    // Doubles in the normal range and their short-mantissa products are
    // exactly representable, so their enclosure is the point itself. A
    // loose enclosure is still sound, which is why only a width check
    // catches it; a loose one made the arc overlay's filter fail on inputs
    // like 10.0, whose enclosure used to be [10, 12].
    for v in [
        1.0,
        -3.5,
        10.0,
        4.0,
        1.5,
        0.001,
        123_456.789,
        1e-280,
        -2.5e270,
    ] {
        let enclosure = Dyadic::from_f64(v).enclosure();
        assert!(enclosure.contains(v), "{v} not in its own enclosure");
        assert_eq!(
            (enclosure.lo(), enclosure.hi()),
            (v, v),
            "{v} enclosed loosely"
        );
    }
    let product = Dyadic::from_f64(10.0).mul(&Dyadic::from_f64(-0.375));
    let enclosure = product.enclosure();
    assert_eq!((enclosure.lo(), enclosure.hi()), (-3.75, -3.75));
    // Far outside the normal range the enclosure may be loose, never wrong.
    for v in [1e-300, -2.5e300] {
        assert!(
            Dyadic::from_f64(v).enclosure().contains(v),
            "{v} not enclosed"
        );
    }
    // Zero is exact: [0, 0].
    let zero = Dyadic::from_f64(0.0).enclosure();
    assert_eq!((zero.lo(), zero.hi()), (0.0, 0.0));
}
