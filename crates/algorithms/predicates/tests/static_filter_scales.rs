//! The static filter at extreme scales (#429): a derived bound that underflows
//! must be refused, and every sign it certifies must equal an exact `Dyadic`
//! reference.

use axiolid_core::{Point2, Point3};
use axiolid_exact::{Arith, Dyadic};
use axiolid_guarantees::Sign;
use axiolid_predicates::StaticFilter;

fn d(x: f64) -> Dyadic {
    Dyadic::from_f64(x)
}

fn orient2(a: Point2, b: Point2, c: Point2) -> Sign {
    let l = d(a.x).sub(&d(c.x)).mul(&d(b.y).sub(&d(c.y)));
    let r = d(a.y).sub(&d(c.y)).mul(&d(b.x).sub(&d(c.x)));
    l.sub(&r).sign().unwrap()
}

fn orient3(a: Point3, b: Point3, c: Point3, p: Point3) -> Sign {
    let v = |q: Point3| {
        [
            d(q.x).sub(&d(p.x)),
            d(q.y).sub(&d(p.y)),
            d(q.z).sub(&d(p.z)),
        ]
    };
    let (a, b, c) = (v(a), v(b), v(c));
    let m = |x: &[Dyadic; 3], y: &[Dyadic; 3]| x[0].mul(&y[1]).sub(&x[1].mul(&y[0]));
    let t0 = a[2].mul(&m(&b, &c));
    let t1 = b[2].mul(&m(&c, &a));
    let t2 = c[2].mul(&m(&a, &b));
    t0.add(&t1).add(&t2).sign().unwrap()
}

/// A range so small that a derived bound leaves the normal range is refused:
/// the bound underflowed to zero and certified every nonzero sign.
#[test]
fn a_range_whose_bound_underflows_is_refused() {
    assert!(StaticFilter::new(1.7440603504673385e-105).is_none());
    assert!(StaticFilter::new(1e-150).is_none());
    assert!(StaticFilter::new(f64::MIN_POSITIVE).is_none());
    assert!(StaticFilter::new(5e-324).is_none());
    assert!(StaticFilter::new(1e-90).is_some());
}

#[test]
fn the_issue_input_is_not_certified_wrong() {
    let a = Point3::new(
        6.585807585185483e-109,
        1.9329642585612275e-113,
        -4.5572251384406943e-113,
    );
    let b = Point3::new(
        8.377599435503294e-111,
        -3.250824966049973e-111,
        -5.626915003093131e-113,
    );
    let c = Point3::new(
        -1.7689163131676088e-109,
        2.425963146421869e-114,
        -1.6861898198154786e-112,
    );
    let p = Point3::new(
        -5.594950297017318e-109,
        1.3324259896260241e-107,
        -8.455405382043586e-107,
    );
    assert_eq!(orient3(a, b, c, p), Sign::Negative);
    assert!(StaticFilter::new(1.7440603504673385e-105).is_none());
}

/// Seeded fuzz over many scales: every certified sign equals the reference,
/// and the filter still answers at the usable scales.
#[test]
fn certified_signs_match_the_dyadic_reference_across_scales() {
    let mut state = 0xD1B5_4A32_D192_ED03u64;
    let mut next = move || {
        state ^= state << 13;
        state ^= state >> 7;
        state ^= state << 17;
        state
    };
    let (mut cases, mut answered, mut refused) = (0usize, 0usize, 0usize);
    for exponent in (-160..=100).step_by(4) {
        let bound = 10f64.powi(exponent);
        let Some(filter) = StaticFilter::new(bound) else {
            refused += 1;
            continue;
        };
        let mut coord = || {
            let unit = (next() >> 11) as f64 / (1u64 << 53) as f64 * 2.0 - 1.0;
            // Some coordinates far below the bound, so differences are tiny.
            let shrink = if next() % 3 == 0 {
                10f64.powi(-((next() % 6) as i32) * 40)
            } else {
                1.0
            };
            unit * bound * shrink
        };
        for _ in 0..1_500 {
            let q = [(); 3].map(|()| Point2::new(coord(), coord()));
            let r = [(); 4].map(|()| Point3::new(coord(), coord(), coord()));
            cases += 2;
            if let Some(s) = filter.orient2d(q[0], q[1], q[2]) {
                assert_eq!(s, orient2(q[0], q[1], q[2]), "orient2d {bound:e} {q:?}");
                answered += 1;
            }
            if let Some(s) = filter.orient3d(r[0], r[1], r[2], r[3]) {
                assert_eq!(
                    s,
                    orient3(r[0], r[1], r[2], r[3]),
                    "orient3d {bound:e} {r:?}"
                );
                answered += 1;
            }
        }
    }
    println!("static fuzz: {answered} certified of {cases}, {refused} ranges refused");
    assert!(answered > 10_000, "answered only {answered}/{cases}");
    assert!(refused > 0);
}
