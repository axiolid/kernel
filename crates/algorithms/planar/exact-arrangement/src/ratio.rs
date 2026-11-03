//! Rounding of an exact ratio to binary64.

use axiolid_exact::Dyadic;
use num_bigint::{BigUint, Sign as BigSign};

/// `numerator / denominator` rounded once to the nearest binary64, ties to
/// even, subnormals included. `None` for a zero denominator or a result
/// beyond the finite range.
pub fn round_ratio(numerator: &Dyadic, denominator: &Dyadic) -> Option<f64> {
    let n = numerator.mantissa().magnitude();
    let d = denominator.mantissa().magnitude();
    if d.bits() == 0 {
        return None;
    }
    if n.bits() == 0 {
        return Some(0.0);
    }
    let negative = (numerator.mantissa().sign() == BigSign::Minus)
        != (denominator.mantissa().sign() == BigSign::Minus);

    // q = floor(n·2^k / d) has 55 or 56 bits; the value is (q + frac)·2^e.
    let k = 55 + d.bits() as i64 - n.bits() as i64;
    let (q, sticky) = if k >= 0 {
        let scaled = n << k as u64;
        (&scaled / d, (&scaled % d).bits() != 0)
    } else {
        let scaled = d << (-k) as u64;
        (n / &scaled, (n % &scaled).bits() != 0)
    };
    let e = numerator.exponent() - denominator.exponent() - k;
    let q_bits = q.bits() as i64;
    // The lowest kept bit: 53 significant bits, or the subnormal grid.
    let lsb = (e + q_bits - 1 - 52).max(-1074);
    let drop = lsb - e;
    if drop > q_bits {
        return Some(0.0);
    }
    let kept = &q >> drop as u64;
    let rest = &q - (&kept << drop as u64);
    let half = BigUint::from(1u8) << (drop - 1) as u64;
    let round_up = match rest.cmp(&half) {
        std::cmp::Ordering::Greater => true,
        std::cmp::Ordering::Equal => sticky || kept.bit(0),
        std::cmp::Ordering::Less => false,
    };
    let mantissa = u64::try_from(&kept).ok()? + u64::from(round_up);
    if mantissa == 0 {
        return Some(0.0);
    }
    let top = lsb + (64 - mantissa.leading_zeros() as i64);
    if top > 1024 {
        return None;
    }
    let magnitude = scale_pow2(mantissa as f64, lsb);
    Some(if negative { -magnitude } else { magnitude })
}

/// `x · 2^e` for an `x` and `e` whose product is exactly representable.
fn scale_pow2(mut x: f64, mut e: i64) -> f64 {
    let pow2 = |s: i64| f64::from_bits(((1023 + s) as u64) << 52);
    while e != 0 {
        let step = e.clamp(-500, 500);
        x *= pow2(step);
        e -= step;
    }
    x
}

#[cfg(test)]
mod tests {
    use axiolid_exact::Arith;

    use super::*;

    fn dyadic(x: f64) -> Dyadic {
        Dyadic::try_from_f64(x).expect("finite")
    }

    /// A small seeded generator; the tests need spread, not quality.
    struct Lcg(u64);

    impl Lcg {
        fn next(&mut self) -> u64 {
            self.0 = self
                .0
                .wrapping_mul(6364136223846793005)
                .wrapping_add(1442695040888963407);
            self.0 >> 11
        }
    }

    #[test]
    fn division_of_doubles_matches_ieee() {
        let mut rng = Lcg(7);
        for _ in 0..4000 {
            let scale = |r: u64| {
                let mantissa = (r >> 3) as f64 + 1.0;
                let shift = (r % 61) as i32 - 30;
                mantissa * 2f64.powi(shift) * if r & 1 == 0 { 1.0 } else { -1.0 }
            };
            let (n, d) = (scale(rng.next()), scale(rng.next()));
            let got = round_ratio(&dyadic(n), &dyadic(d)).expect("in range");
            assert!(got.to_bits() == (n / d).to_bits(), "{n} / {d}");
        }
    }

    #[test]
    fn ties_round_to_even_and_sticky_bits_round_up() {
        let ulp = f64::EPSILON;
        // 1 + 2^-53 is the tie between 1 and 1 + 2^-52: even is 1.
        let tie = dyadic(1.0).add(&dyadic(ulp / 2.0));
        assert_eq!(round_ratio(&tie, &dyadic(1.0)), Some(1.0));
        // 1 + 3·2^-53 is the tie between 1 + 2^-52 and 1 + 2^-51: even is the latter.
        let tie = dyadic(1.0).add(&dyadic(1.5 * ulp));
        assert_eq!(round_ratio(&tie, &dyadic(1.0)), Some(1.0 + 2.0 * ulp));
        // Just above the first tie rounds up.
        let above = dyadic(1.0)
            .add(&dyadic(ulp / 2.0))
            .add(&dyadic(2f64.powi(-200)));
        assert_eq!(round_ratio(&above, &dyadic(1.0)), Some(1.0 + ulp));
        // The same tie plus a remainder only the sticky bit sees: the
        // quotient needs a division that does not terminate.
        let sticky = dyadic(3.0)
            .add(&dyadic(1.5 * ulp))
            .add(&dyadic(2f64.powi(-200)));
        assert_eq!(round_ratio(&sticky, &dyadic(3.0)), Some(1.0 + ulp));
        // Without the remainder it is the tie again, and rounds to even.
        let tie = dyadic(3.0).add(&dyadic(1.5 * ulp));
        assert_eq!(round_ratio(&tie, &dyadic(3.0)), Some(1.0));
    }

    #[test]
    fn subnormals_and_underflow_round_correctly() {
        let min = f64::from_bits(1);
        assert_eq!(round_ratio(&dyadic(min), &dyadic(1.0)), Some(min));
        // Exactly half the smallest subnormal ties to even (zero)...
        let half = Dyadic::from_parts(1.into(), -1075);
        assert_eq!(round_ratio(&half, &dyadic(1.0)), Some(0.0));
        // ...and anything above rounds up to it.
        let above = half.add(&Dyadic::from_parts(1.into(), -1200));
        assert_eq!(round_ratio(&above, &dyadic(1.0)), Some(min));
        // Three halves of the smallest subnormal ties to the even 2·min.
        let three_halves = Dyadic::from_parts(3.into(), -1075);
        assert_eq!(round_ratio(&three_halves, &dyadic(1.0)), Some(2.0 * min));
        // Far below the grid.
        let tiny = Dyadic::from_parts(1.into(), -3000);
        assert_eq!(round_ratio(&tiny, &dyadic(1.0)), Some(0.0));
    }

    #[test]
    fn overflow_and_zero_denominator_are_refused() {
        let huge = Dyadic::from_parts(1.into(), 1024);
        assert_eq!(round_ratio(&huge, &dyadic(1.0)), None);
        assert_eq!(round_ratio(&dyadic(1.0), &Dyadic::zero()), None);
        let max = dyadic(f64::MAX);
        assert_eq!(round_ratio(&max, &dyadic(1.0)), Some(f64::MAX));
        // Rounds up past MAX.
        let past = max.add(&Dyadic::from_parts(1.into(), 970));
        assert_eq!(round_ratio(&past, &dyadic(1.0)), None);
    }

    #[test]
    fn signs_combine_and_zero_is_positive_zero() {
        assert_eq!(round_ratio(&dyadic(-3.0), &dyadic(4.0)), Some(-0.75));
        assert_eq!(round_ratio(&dyadic(3.0), &dyadic(-4.0)), Some(-0.75));
        assert_eq!(round_ratio(&dyadic(-3.0), &dyadic(-4.0)), Some(0.75));
        let zero = round_ratio(&Dyadic::zero(), &dyadic(-4.0)).expect("zero");
        assert_eq!(zero.to_bits(), 0.0f64.to_bits());
    }
}
