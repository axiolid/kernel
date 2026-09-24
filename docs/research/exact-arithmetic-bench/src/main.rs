//! Exact construction workload: intersect two segments with exact
//! rationals (inputs are f64, converted exactly), then take the sign of
//! orient2d(a, b, intersection). Same inputs for both libraries; the
//! signs must agree. Timing is the median of several full passes.
use std::time::Instant;

/// Deterministic xorshift so both libraries see identical inputs.
struct Rng(u64);
impl Rng {
    fn f(&mut self) -> f64 {
        self.0 ^= self.0 << 13;
        self.0 ^= self.0 >> 7;
        self.0 ^= self.0 << 17;
        // Coordinates in [-1000, 1000) with full 53-bit mantissas.
        ((self.0 >> 11) as f64 / (1u64 << 53) as f64) * 2000.0 - 1000.0
    }
}

fn inputs(n: usize) -> Vec<[f64; 10]> {
    let mut r = Rng(0x9E37_79B9_7F4A_7C15);
    (0..n).map(|_| std::array::from_fn(|_| r.f())).collect()
}

mod num {
    use num_bigint::BigInt;
    use num_rational::BigRational as Q;
    use num_traits::{Signed, Zero};

    fn q(x: f64) -> Q {
        Q::from_float(x).expect("finite")
    }

    /// Returns sign of orient2d(a, b, p) where p = seg(p1,p2) x seg(p3,p4).
    pub fn run(v: &[f64; 10]) -> i8 {
        let [x1, y1, x2, y2, x3, y3, x4, y4, ax, ay] = v.map(q);
        let d = (x1.clone() - &x2) * (y3.clone() - &y4) - (y1.clone() - &y2) * (x3.clone() - &x4);
        if d.is_zero() {
            return 0;
        }
        let c1 = x1.clone() * &y2 - y1.clone() * &x2;
        let c2 = x3.clone() * &y4 - y3.clone() * &x4;
        let px = (c1.clone() * (x3.clone() - &x4) - (x1.clone() - &x2) * &c2) / &d;
        let py = (c1 * (y3 - y4) - (y1.clone() - &y2) * c2) / &d;
        let o = (x1.clone() - &ax) * (py - &ay) - (y1 - &ay) * (px - &ax);
        let _ = BigInt::zero();
        if o.is_zero() {
            0
        } else if o.is_positive() {
            1
        } else {
            -1
        }
    }
}

mod dsh {
    use dashu_ratio::RBig as Q;

    fn q(x: f64) -> Q {
        Q::try_from(x).expect("finite")
    }

    pub fn run(v: &[f64; 10]) -> i8 {
        let [x1, y1, x2, y2, x3, y3, x4, y4, ax, ay] = v.map(q);
        let d = (&x1 - &x2) * (&y3 - &y4) - (&y1 - &y2) * (&x3 - &x4);
        if d == Q::ZERO {
            return 0;
        }
        let c1 = &x1 * &y2 - &y1 * &x2;
        let c2 = &x3 * &y4 - &y3 * &x4;
        let px = (&c1 * (&x3 - &x4) - (&x1 - &x2) * &c2) / &d;
        let py = (&c1 * (&y3 - &y4) - (&y1 - &y2) * &c2) / &d;
        let o = (&x1 - &ax) * (&py - &ay) - (&y1 - &ay) * (&px - &ax);
        if o == Q::ZERO {
            0
        } else if o > Q::ZERO {
            1
        } else {
            -1
        }
    }
}

fn time(label: &str, data: &[[f64; 10]], f: fn(&[f64; 10]) -> i8) -> Vec<i8> {
    let mut signs = Vec::new();
    let mut runs = Vec::new();
    for _ in 0..7 {
        let t = Instant::now();
        signs = data.iter().map(f).collect();
        runs.push(t.elapsed().as_secs_f64());
    }
    runs.sort_by(|a, b| a.partial_cmp(b).unwrap());
    let med = runs[runs.len() / 2];
    println!(
        "{label:<22} median {:>8.2} ms  ({:.2} us/op, min {:.2} ms, max {:.2} ms)",
        med * 1e3,
        med * 1e6 / data.len() as f64,
        runs[0] * 1e3,
        runs[runs.len() - 1] * 1e3
    );
    signs
}

fn main() {
    let data = inputs(20_000);
    let a = time("num-rational 0.4.2", &data, num::run);
    let b = time("dashu-ratio 0.6", &data, dsh::run);
    let c = time("num-bigint 0.4 (int)", &data, numi::run);
    let d = time("dashu-int 0.6 (int)", &data, dshi::run);
    let e = time("num-bigint 0.5 (int)", &data, numi5::run);
    assert_eq!(a, e, "num-bigint 0.5 integer form disagrees");
    assert_eq!(a, c, "num integer form disagrees with rational form");
    assert_eq!(a, d, "dashu integer form disagrees with rational form");
    let agree = a.iter().zip(&b).filter(|(x, y)| x == y).count();
    let pos = a.iter().filter(|s| **s > 0).count();
    println!("signs agree: {agree}/{}  (positive {pos})", a.len());
    assert_eq!(a, b, "libraries disagree on an exact sign");
}

/// Each input as (mantissa, shift) on one common grid: x = m * 2^(e - emin)
/// times 2^emin. Scaling every coordinate by the same positive factor
/// leaves every sign below unchanged, so the grid factor is dropped.
fn grid(v: &[f64; 10]) -> [(i64, usize); 10] {
    let parts = v.map(|x| {
        let bits = x.to_bits();
        let neg = bits >> 63 == 1;
        let exp = ((bits >> 52) & 0x7ff) as i32;
        let frac = (bits & ((1u64 << 52) - 1)) as i64;
        let (m, e) = if exp == 0 {
            (frac, -1074)
        } else {
            (frac | (1 << 52), exp - 1075)
        };
        (if neg { -m } else { m }, e)
    });
    let emin = parts.iter().map(|p| p.1).min().unwrap();
    parts.map(|(m, e)| (m, (e - emin) as usize))
}

/// Division-free exact sign: with p = (Nx/D, Ny/D),
/// sign(orient) = sign(D) * sign((x1-ax)(Ny - ay D) - (y1-ay)(Nx - ax D)).
/// No gcd, no rational normalisation: integers only.
mod numi {
    use num_bigint::BigInt;
    use num_traits::{Signed, Zero};

    pub fn run(v: &[f64; 10]) -> i8 {
        let [x1, y1, x2, y2, x3, y3, x4, y4, ax, ay] =
            super::grid(v).map(|(m, s)| BigInt::from(m) << s);
        let d = (&x1 - &x2) * (&y3 - &y4) - (&y1 - &y2) * (&x3 - &x4);
        if d.is_zero() {
            return 0;
        }
        let c1 = &x1 * &y2 - &y1 * &x2;
        let c2 = &x3 * &y4 - &y3 * &x4;
        let nx = &c1 * (&x3 - &x4) - (&x1 - &x2) * &c2;
        let ny = &c1 * (&y3 - &y4) - (&y1 - &y2) * &c2;
        let od = (&x1 - &ax) * (ny - &ay * &d) - (&y1 - &ay) * (nx - &ax * &d);
        let s = if od.is_zero() {
            0
        } else if od.is_positive() {
            1
        } else {
            -1
        };
        if d.is_positive() {
            s
        } else {
            -s
        }
    }
}

mod dshi {
    use dashu_int::IBig;

    pub fn run(v: &[f64; 10]) -> i8 {
        let [x1, y1, x2, y2, x3, y3, x4, y4, ax, ay] =
            super::grid(v).map(|(m, s)| IBig::from(m) << s);
        let d = (&x1 - &x2) * (&y3 - &y4) - (&y1 - &y2) * (&x3 - &x4);
        if d == IBig::ZERO {
            return 0;
        }
        let c1 = &x1 * &y2 - &y1 * &x2;
        let c2 = &x3 * &y4 - &y3 * &x4;
        let nx = &c1 * (&x3 - &x4) - (&x1 - &x2) * &c2;
        let ny = &c1 * (&y3 - &y4) - (&y1 - &y2) * &c2;
        let od = (&x1 - &ax) * (ny - &ay * &d) - (&y1 - &ay) * (nx - &ax * &d);
        let s = if od == IBig::ZERO {
            0
        } else if od > IBig::ZERO {
            1
        } else {
            -1
        };
        if d > IBig::ZERO {
            s
        } else {
            -s
        }
    }
}

mod numi5 {
    use nb5::BigInt;

    pub fn run(v: &[f64; 10]) -> i8 {
        let [x1, y1, x2, y2, x3, y3, x4, y4, ax, ay] =
            super::grid(v).map(|(m, s)| BigInt::from(m) << s);
        let d = (&x1 - &x2) * (&y3 - &y4) - (&y1 - &y2) * (&x3 - &x4);
        let zero = BigInt::from(0);
        if d == zero {
            return 0;
        }
        let c1 = &x1 * &y2 - &y1 * &x2;
        let c2 = &x3 * &y4 - &y3 * &x4;
        let nx = &c1 * (&x3 - &x4) - (&x1 - &x2) * &c2;
        let ny = &c1 * (&y3 - &y4) - (&y1 - &y2) * &c2;
        let od = (&x1 - &ax) * (ny - &ay * &d) - (&y1 - &ay) * (nx - &ax * &d);
        let s = if od == zero {
            0
        } else if od > zero {
            1
        } else {
            -1
        };
        if d > zero {
            s
        } else {
            -s
        }
    }
}
