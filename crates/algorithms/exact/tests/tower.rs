//! Nested square roots against an independent fixed-point oracle.

use axiolid_exact::{Arith, Dyadic, Interval, Nested, Tower};
use axiolid_guarantees::Sign;
use num_bigint::BigInt;

/// Fixed-point bits for the oracle: far beyond anything the tests need.
const F: u64 = 2000;

/// A value in the oracle: `round(x * 2^F)`, computed independently of the
/// crate (plain big-integer arithmetic and integer square roots).
#[derive(Clone)]
struct Fx(BigInt);

impl Fx {
    fn int(value: i64) -> Self {
        Fx(BigInt::from(value) << F)
    }
    fn add(&self, o: &Self) -> Self {
        Fx(&self.0 + &o.0)
    }
    fn sub(&self, o: &Self) -> Self {
        Fx(&self.0 - &o.0)
    }
    fn mul(&self, o: &Self) -> Self {
        Fx((&self.0 * &o.0) >> F)
    }
    fn sqrt(&self) -> Self {
        assert!(self.0 >= BigInt::from(0), "oracle radicand negative");
        Fx((&self.0 << F).sqrt())
    }
    /// Sign, or `None` inside the oracle's own error band.
    fn sign(&self) -> Option<Sign> {
        let band = BigInt::from(1) << (F / 2);
        if self.0 > band {
            Some(Sign::Positive)
        } else if self.0 < -band {
            Some(Sign::Negative)
        } else {
            None
        }
    }
}

/// One expression built in any arithmetic, so the same construction runs
/// in the crate's two tiers and in the oracle.
trait Build {
    type V: Clone;
    fn int(&mut self, v: i64) -> Self::V;
    fn add(&mut self, x: &Self::V, y: &Self::V) -> Self::V;
    fn sub(&mut self, x: &Self::V, y: &Self::V) -> Self::V;
    fn mul(&mut self, x: &Self::V, y: &Self::V) -> Self::V;
    fn sqrt(&mut self, x: &Self::V) -> Self::V;
}

impl<T: Arith> Build for Tower<T> {
    type V = Nested<T>;
    fn int(&mut self, v: i64) -> Nested<T> {
        self.from_f64(v as f64)
    }
    fn add(&mut self, x: &Nested<T>, y: &Nested<T>) -> Nested<T> {
        Tower::add(self, x, y)
    }
    fn sub(&mut self, x: &Nested<T>, y: &Nested<T>) -> Nested<T> {
        Tower::sub(self, x, y)
    }
    fn mul(&mut self, x: &Nested<T>, y: &Nested<T>) -> Nested<T> {
        Tower::mul(self, x, y)
    }
    fn sqrt(&mut self, x: &Nested<T>) -> Nested<T> {
        Tower::sqrt(self, x).expect("test towers stay shallow")
    }
}

struct Oracle;

impl Build for Oracle {
    type V = Fx;
    fn int(&mut self, v: i64) -> Fx {
        Fx::int(v)
    }
    fn add(&mut self, x: &Fx, y: &Fx) -> Fx {
        x.add(y)
    }
    fn sub(&mut self, x: &Fx, y: &Fx) -> Fx {
        x.sub(y)
    }
    fn mul(&mut self, x: &Fx, y: &Fx) -> Fx {
        x.mul(y)
    }
    fn sqrt(&mut self, x: &Fx) -> Fx {
        x.sqrt()
    }
}

fn exact_sign(build: impl Fn(&mut Tower<Dyadic>) -> Nested<Dyadic>) -> Sign {
    let mut tower = Tower::new();
    let value = build(&mut tower);
    tower.sign(&value).expect("exact arithmetic always decides")
}

/// sqrt(3 + 2*sqrt(2)) - (1 + sqrt(2)): zero, by (1 + sqrt 2)^2.
fn denest_two<B: Build>(b: &mut B) -> B::V {
    let two = b.int(2);
    let three = b.int(3);
    let one = b.int(1);
    let s2 = b.sqrt(&two);
    let inner = b.mul(&two, &s2);
    let inner = b.add(&three, &inner);
    let t = b.sqrt(&inner);
    let rhs = b.add(&one, &s2);
    b.sub(&t, &rhs)
}

/// sqrt(5 + 2*sqrt(6)) - sqrt(2) - sqrt(3), with sqrt(6) written as
/// sqrt(2)*sqrt(3): zero, and it needs three radicals to say so.
fn denest_three<B: Build>(b: &mut B) -> B::V {
    let two = b.int(2);
    let three = b.int(3);
    let five = b.int(5);
    let s2 = b.sqrt(&two);
    let s3 = b.sqrt(&three);
    let s6 = b.mul(&s2, &s3);
    let inner = b.mul(&two, &s6);
    let inner = b.add(&five, &inner);
    let t = b.sqrt(&inner);
    let rhs = b.add(&s2, &s3);
    b.sub(&t, &rhs)
}

#[test]
fn denesting_identities_are_exact_zeros() {
    assert_eq!(exact_sign(denest_two), Sign::Zero);
    assert_eq!(exact_sign(denest_three), Sign::Zero);
    // sqrt(8) - 2*sqrt(2): dependent radicals, non-zero coefficients.
    let dependent = exact_sign(|t| {
        let eight = t.from_f64(8.0);
        let two = t.from_f64(2.0);
        let s8 = t.sqrt(&eight).unwrap();
        let s2 = t.sqrt(&two).unwrap();
        let twice = Tower::mul(t, &two, &s2);
        Tower::sub(t, &s8, &twice)
    });
    assert_eq!(dependent, Sign::Zero);
}

#[test]
fn a_perturbation_far_below_f64_resolution_is_seen() {
    // denest_two + 2^-100 and - 2^-100: doubles see both as zero.
    for (offset, want) in [(1.0, Sign::Positive), (-1.0, Sign::Negative)] {
        let got = exact_sign(|t| {
            let zero = denest_two(t);
            let tiny = t.from_f64(offset * 2f64.powi(-100));
            Tower::add(t, &zero, &tiny)
        });
        assert_eq!(got, want);
    }
    // The filter must not decide these (it cannot see 2^-100 next to
    // values of size ~2.4), and must never decide wrongly.
    let mut tower = Tower::<Interval>::new();
    let zero = denest_two(&mut tower);
    let tiny = tower.from_f64(2f64.powi(-100));
    let near = Tower::add(&tower, &zero, &tiny);
    assert_eq!(tower.sign(&near), None);
}

/// A deterministic generator (xorshift), so failures reproduce.
struct Rng(u64);

impl Rng {
    fn next(&mut self) -> u64 {
        self.0 ^= self.0 << 13;
        self.0 ^= self.0 >> 7;
        self.0 ^= self.0 << 17;
        self.0
    }
    fn small(&mut self) -> i64 {
        (self.next() % 9) as i64 - 4
    }
}

/// A random program: a few radicals over earlier values, then a random
/// polynomial in everything. Radicands are forced non-negative as
/// `x*x + k` with `k >= 0`, so the program is valid in every arithmetic.
fn program<B: Build>(b: &mut B, seed: u64) -> B::V {
    let mut rng = Rng(seed);
    let mut pool = vec![b.int(rng.small()), b.int(rng.small())];
    let radicals = 1 + (rng.next() % 3) as usize;
    for _ in 0..radicals {
        let pick = pool[(rng.next() as usize) % pool.len()].clone();
        let k = b.int((rng.next() % 4) as i64);
        let square = b.mul(&pick, &pick);
        let radicand = b.add(&square, &k);
        let root = b.sqrt(&radicand);
        pool.push(root);
    }
    let mut acc = b.int(rng.small());
    for _ in 0..4 {
        let x = pool[(rng.next() as usize) % pool.len()].clone();
        let y = pool[(rng.next() as usize) % pool.len()].clone();
        let c = b.int(rng.small());
        let term = b.mul(&x, &y);
        let term = b.mul(&term, &c);
        acc = if rng.next() % 2 == 0 {
            b.add(&acc, &term)
        } else {
            b.sub(&acc, &term)
        };
    }
    acc
}

#[test]
fn random_programs_agree_with_the_oracle() {
    let (mut decided, mut zeros, mut filtered) = (0, 0, 0);
    for seed in 1..=3000u64 {
        let oracle = program(&mut Oracle, seed).sign();
        let mut exact_tower = Tower::<Dyadic>::new();
        let exact_value = program(&mut exact_tower, seed);
        let exact = exact_tower.sign(&exact_value).expect("exact decides");
        let mut fast_tower = Tower::<Interval>::new();
        let fast_value = program(&mut fast_tower, seed);
        if let Some(fast) = fast_tower.sign(&fast_value) {
            assert_eq!(fast, exact, "filter disagrees with exact, seed {seed}");
            filtered += 1;
        }
        match oracle {
            Some(sign) => {
                assert_eq!(exact, sign, "exact disagrees with oracle, seed {seed}");
                decided += 1;
            }
            // Inside the oracle's 2^-1000 band: only an exact zero is
            // plausible for integer programs this small.
            None => {
                assert_eq!(exact, Sign::Zero, "seed {seed}");
                zeros += 1;
            }
        }
    }
    // Guard against a vacuous pass, and against the filter regressing to
    // coefficient case analysis only (which decided 892 of 3000 here; the
    // numeric enclosure decides about 2860).
    assert!(decided > 2000, "only {decided} decided");
    assert!(zeros > 20, "only {zeros} exact zeros");
    assert!(filtered > 2700, "filter decided only {filtered}");
}

#[test]
fn depth_is_capped() {
    let mut tower = Tower::<Dyadic>::new();
    let two = tower.from_f64(2.0);
    for _ in 0..axiolid_exact::tower::MAX_DEPTH {
        tower.sqrt(&two).unwrap();
    }
    assert_eq!(tower.sqrt(&two), Err(axiolid_exact::ExactError::TooDeep));
}

#[test]
fn a_negative_radicand_is_undefined_not_a_sign() {
    let mut tower = Tower::<Dyadic>::new();
    let minus = tower.from_f64(-1.0);
    let root = tower.sqrt(&minus).unwrap();
    assert_eq!(tower.sign(&root), None);
}
