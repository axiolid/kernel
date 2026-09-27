//! `incircle` and `insphere` against an exact dyadic reference, on nearly
//! cocircular and cospherical points whose coordinate differences do not
//! fit an `f64` (#190: the exact fallback rounded them first).

use axiolid_core::{Point2, Point3};
use axiolid_exact::{Arith, Dyadic};
use axiolid_guarantees::{Certified, Sign};
use axiolid_predicates::{incircle, insphere};

fn d(x: f64) -> Dyadic {
    Dyadic::from_f64(x)
}

fn sign(c: Certified) -> Sign {
    match c {
        Certified::Certain { sign, .. } => sign,
        other => panic!("{other:?}"),
    }
}

fn det3(m: [[Dyadic; 3]; 3]) -> Dyadic {
    let t = |a: &Dyadic, b: &Dyadic, c: &Dyadic, e: &Dyadic| a.mul(b).sub(&c.mul(e));
    m[0][0]
        .mul(&t(&m[1][1], &m[2][2], &m[1][2], &m[2][1]))
        .sub(&m[0][1].mul(&t(&m[1][0], &m[2][2], &m[1][2], &m[2][0])))
        .add(&m[0][2].mul(&t(&m[1][0], &m[2][1], &m[1][1], &m[2][0])))
}

fn incircle_reference(a: Point2, b: Point2, c: Point2, p: Point2) -> Sign {
    let row = |q: Point2| {
        let (x, y) = (d(q.x).sub(&d(p.x)), d(q.y).sub(&d(p.y)));
        let lift = x.mul(&x).add(&y.mul(&y));
        [x, y, lift]
    };
    det3([row(a), row(b), row(c)]).sign().unwrap()
}

fn insphere_reference(a: Point3, b: Point3, c: Point3, e: Point3, p: Point3) -> Sign {
    let row = |q: Point3| {
        let v = [
            d(q.x).sub(&d(p.x)),
            d(q.y).sub(&d(p.y)),
            d(q.z).sub(&d(p.z)),
        ];
        let lift = v[0].mul(&v[0]).add(&v[1].mul(&v[1])).add(&v[2].mul(&v[2]));
        [v[0].clone(), v[1].clone(), v[2].clone(), lift]
    };
    let m = [row(a), row(b), row(c), row(e)];
    // 4x4 determinant by expansion along the first row.
    let minor = |skip: usize| {
        let pick = |r: usize| {
            let v: Vec<Dyadic> = (0..4)
                .filter(|&k| k != skip)
                .map(|k| m[r][k].clone())
                .collect();
            [v[0].clone(), v[1].clone(), v[2].clone()]
        };
        det3([pick(1), pick(2), pick(3)])
    };
    let mut total = Dyadic::from_f64(0.0);
    for k in 0..4 {
        let term = m[0][k].mul(&minor(k));
        total = if k % 2 == 0 {
            total.add(&term)
        } else {
            total.sub(&term)
        };
    }
    total.sign().unwrap()
}

fn xorshift(s: &mut u64) -> f64 {
    *s ^= *s << 13;
    *s ^= *s >> 7;
    *s ^= *s << 17;
    (*s >> 11) as f64 / (1u64 << 53) as f64
}

#[test]
fn incircle_matches_the_exact_determinant_on_nearly_cocircular_points() {
    let mut s = 0x1234_5678_9abc_def1u64;
    let mut decided = 0;
    for _ in 0..20000 {
        // Points on a circle of decimal centre and radius, rounded: nearly
        // cocircular, so the filter hands over to the exact path.
        let (cx, cy, r) = (
            xorshift(&mut s) * 7.3,
            xorshift(&mut s) * 3.1,
            0.3 + xorshift(&mut s),
        );
        let mut on = || {
            let t = xorshift(&mut s) * std::f64::consts::TAU;
            Point2::new(cx + r * t.cos(), cy + r * t.sin())
        };
        let (a, b, c, p) = (on(), on(), on(), on());
        let orient = (b - a).perp_dot(c - a);
        let (a, b) = if orient > 0.0 { (a, b) } else { (b, a) };
        let expected = incircle_reference(a, b, c, p);
        assert_eq!(
            sign(incircle(a, b, c, p)),
            expected,
            "{a:?} {b:?} {c:?} {p:?}"
        );
        decided += 1;
    }
    assert_eq!(decided, 20000);
    // The square turned 45 degrees from #190: its samples are exactly
    // cocircular in fours.
    let q = |x: f64, y: f64| Point2::new(x, y);
    let (a, b, c, p) = (
        q(0.9655172413793104, -0.9655172413793104),
        q(1.0344827586206895, 1.0344827586206895),
        q(0.9655172413793103, 0.9655172413793103),
        q(3.0344827586206895, 0.9655172413793104),
    );
    let (a, b) = if (b - a).perp_dot(c - a) > 0.0 {
        (a, b)
    } else {
        (b, a)
    };
    assert_eq!(sign(incircle(a, b, c, p)), incircle_reference(a, b, c, p));
}

#[test]
fn insphere_matches_the_exact_determinant_on_nearly_cospherical_points() {
    let mut s = 0x0fed_cba9_8765_4321u64;
    for _ in 0..5000 {
        let (cx, cy, cz, r) = (
            xorshift(&mut s) * 5.1,
            xorshift(&mut s) * 2.3,
            xorshift(&mut s) * 1.7,
            0.3 + xorshift(&mut s),
        );
        let mut on = || {
            let (u, v) = (
                xorshift(&mut s) * std::f64::consts::TAU,
                xorshift(&mut s) * 2.0 - 1.0,
            );
            let w = (1.0 - v * v).sqrt();
            Point3::new(cx + r * w * u.cos(), cy + r * w * u.sin(), cz + r * v)
        };
        let (a, b, c, e, p) = (on(), on(), on(), on(), on());
        assert_eq!(
            sign(insphere(a, b, c, e, p)),
            insphere_reference(a, b, c, e, p),
            "{a:?} {b:?} {c:?} {e:?} {p:?}"
        );
    }
}
