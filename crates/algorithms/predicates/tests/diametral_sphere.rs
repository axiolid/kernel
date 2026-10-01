//! Differential gates for `in_diametral_sphere`.
//!
//! The oracle is the same degree-6 polynomial evaluated exactly in `i128`
//! over integer coordinates, independent of the expansion arithmetic under
//! test, plus a geometric oracle (squared distances to the circumcentre in
//! rational arithmetic) for the coplanar case the 3D Delaunay hull test uses.

use axiolid_core::Point3;
use axiolid_exact::{Arith, Dyadic};
use axiolid_guarantees::Sign;
use axiolid_predicates::{in_diametral_sphere, in_diametral_sphere_filter};

type V = [i128; 3];

fn sub(p: V, q: V) -> V {
    [p[0] - q[0], p[1] - q[1], p[2] - q[2]]
}

fn cross(p: V, q: V) -> V {
    [
        p[1] * q[2] - p[2] * q[1],
        p[2] * q[0] - p[0] * q[2],
        p[0] * q[1] - p[1] * q[0],
    ]
}

fn dot(p: V, q: V) -> i128 {
    p[0] * q[0] + p[1] * q[1] + p[2] * q[2]
}

/// Inside the diametral sphere iff `|d - centre|^2 < r^2`, written over the
/// common denominator `|n|^2` so it stays in integers.
fn oracle(a: V, b: V, c: V, d: V) -> Sign {
    let (u, v, w) = (sub(b, a), sub(c, a), sub(d, a));
    let n = cross(u, v);
    let power =
        dot(n, n) * dot(w, w) - dot(u, u) * dot(cross(w, v), n) - dot(v, v) * dot(cross(u, w), n);
    match (-power).signum() {
        1 => Sign::Positive,
        -1 => Sign::Negative,
        _ => Sign::Zero,
    }
}

fn rng(state: &mut u64) -> i64 {
    *state ^= *state << 13;
    *state ^= *state >> 7;
    *state ^= *state << 17;
    (*state >> 33) as i64
}

fn point(v: V) -> Point3 {
    Point3::new(v[0] as f64, v[1] as f64, v[2] as f64)
}

/// Random and constructed-cocircular cases agree with the integer oracle.
#[test]
fn in_diametral_sphere_agrees_with_an_independent_exact_oracle() {
    let mut state = 0x1234_5678_9ABC_DEF1;
    let bound = 200i64;
    let (mut zeros, mut positives, mut negatives) = (0usize, 0usize, 0usize);
    // Two orthogonal lattice vectors of equal length 3: centre +- s, +- t lie
    // exactly on one circle in a plane that is not axis-aligned.
    let s: V = [1, 2, 2];
    let t: V = [2, 1, -2];
    for i in 0..20_000 {
        let mut pt = || {
            [
                i128::from(rng(&mut state) % (2 * bound) - bound),
                i128::from(rng(&mut state) % (2 * bound) - bound),
                i128::from(rng(&mut state) % (2 * bound) - bound),
            ]
        };
        let (a, b, c, d) = match i % 3 {
            0 => (pt(), pt(), pt(), pt()),
            1 => {
                let centre = pt();
                let k = i128::from(rng(&mut state) % 5 + 1);
                let on = |sign: i128, v: V| {
                    [
                        centre[0] + sign * k * v[0],
                        centre[1] + sign * k * v[1],
                        centre[2] + sign * k * v[2],
                    ]
                };
                (on(1, s), on(1, t), on(-1, s), on(-1, t))
            }
            _ => {
                // Coplanar but generally not cocircular: d in the triangle's
                // plane at an integer combination of its edge vectors.
                let (a, b, c) = (pt(), pt(), pt());
                let (x, y) = (
                    i128::from(rng(&mut state) % 5 - 2),
                    i128::from(rng(&mut state) % 5 - 2),
                );
                let (u, v) = (sub(b, a), sub(c, a));
                let d = [
                    a[0] + x * u[0] + y * v[0],
                    a[1] + x * u[1] + y * v[1],
                    a[2] + x * u[2] + y * v[2],
                ];
                (a, b, c, d)
            }
        };
        let got = in_diametral_sphere(point(a), point(b), point(c), point(d))
            .sign()
            .expect("in range, so certified");
        let want = oracle(a, b, c, d);
        assert_eq!(got, want, "disagreement on {a:?} {b:?} {c:?} {d:?}");
        match want {
            Sign::Zero => zeros += 1,
            Sign::Positive => positives += 1,
            _ => negatives += 1,
        }
    }
    assert!(zeros > 1_000, "only {zeros} cocircular cases");
    assert!(positives > 1_000 && negatives > 1_000);
}

/// Coplanar `d` is decided against the circumcircle, independently of the
/// triangle's orientation.
#[test]
fn coplanar_points_are_decided_against_the_circumcircle() {
    let p = |x: f64, y: f64, z: f64| Point3::new(x, y, z);
    // Unit right triangle in z = 3: circumcentre (0.5, 0.5, 3), r^2 = 0.5.
    let (a, b, c) = (p(0.0, 0.0, 3.0), p(1.0, 0.0, 3.0), p(0.0, 1.0, 3.0));
    for (tri, name) in [((a, b, c), "abc"), ((a, c, b), "acb")] {
        let s = |d| in_diametral_sphere(tri.0, tri.1, tri.2, d).sign();
        assert_eq!(s(p(0.5, 0.5, 3.0)), Some(Sign::Positive), "{name}");
        assert_eq!(s(p(1.0, 1.0, 3.0)), Some(Sign::Zero), "{name}");
        assert_eq!(s(p(2.0, 2.0, 3.0)), Some(Sign::Negative), "{name}");
        // Off the plane: the sphere, not the infinite cylinder, decides.
        assert_eq!(s(p(0.5, 0.5, 3.5)), Some(Sign::Positive), "{name}");
        assert_eq!(s(p(0.5, 0.5, 4.0)), Some(Sign::Negative), "{name}");
    }
    // A degenerate (collinear) triangle has no sphere.
    assert_eq!(
        in_diametral_sphere(a, b, p(2.0, 0.0, 3.0), p(0.3, 0.1, 3.0)).sign(),
        Some(Sign::Zero)
    );
}

/// The filter must defer on exact cocircularity, and the exact path must then
/// find the sign one ULP away from it.
#[test]
fn the_exact_path_decides_near_cocircular_cases() {
    let ulp = |v: f64| f64::from_bits(v.to_bits() + 1);
    let (mut deferred, mut definite) = (0usize, 0usize);
    for k in 1..1_000i64 {
        let r = k as f64;
        // A rectangle in the tilted plane z = y: its corners are cocircular.
        let (a, b, c) = (
            Point3::new(-r, -r, -r),
            Point3::new(r, -r, -r),
            Point3::new(r, r, r),
        );
        let d_on = Point3::new(-r, r, r);
        assert_eq!(
            in_diametral_sphere(a, b, c, d_on).sign(),
            Some(Sign::Zero),
            "the rectangle's corners are cocircular"
        );
        if !in_diametral_sphere_filter(a, b, c, d_on).is_certain() {
            deferred += 1;
        }
        // One ULP further out, then one ULP further in.
        let d_out = Point3::new(-ulp(r), r, r);
        let d_in = Point3::new(-r.next_down(), r, r);
        if in_diametral_sphere(a, b, c, d_out).sign() == Some(Sign::Negative)
            && in_diametral_sphere(a, b, c, d_in).sign() == Some(Sign::Positive)
        {
            definite += 1;
        }
    }
    assert!(deferred > 900, "only {deferred} deferred");
    assert!(definite > 900, "only {definite} definite");
}

/// Coordinates outside the exact path's range are refused, not guessed.
#[test]
fn out_of_range_input_is_uncertain() {
    let p = |x: f64| Point3::new(x, 0.0, 0.0);
    let q = Point3::new(0.0, 1.0, 0.0);
    let r = Point3::new(0.0, 0.0, 1.0);
    assert!(!in_diametral_sphere(p(f64::NAN), q, r, p(0.5)).is_certain());
    // Cocircular, so the filter cannot settle it, and 2^200 is out of range.
    let big = 2f64.powi(200);
    let s = Point3::new(big, 0.0, 0.0);
    let t = Point3::new(-big, 0.0, 0.0);
    let u = Point3::new(0.0, big, 0.0);
    let v = Point3::new(0.0, -big, 0.0);
    assert!(!in_diametral_sphere(s, t, u, v).is_certain());
}

/// The same polynomial over exact dyadic rationals, for coordinates whose
/// differences do not fit an `f64`.
fn dyadic_reference(a: Point3, b: Point3, c: Point3, d: Point3) -> Sign {
    let vector = |p: Point3| {
        [
            Dyadic::from_f64(p.x).sub(&Dyadic::from_f64(a.x)),
            Dyadic::from_f64(p.y).sub(&Dyadic::from_f64(a.y)),
            Dyadic::from_f64(p.z).sub(&Dyadic::from_f64(a.z)),
        ]
    };
    let cross = |p: &[Dyadic; 3], q: &[Dyadic; 3]| {
        [
            p[1].mul(&q[2]).sub(&p[2].mul(&q[1])),
            p[2].mul(&q[0]).sub(&p[0].mul(&q[2])),
            p[0].mul(&q[1]).sub(&p[1].mul(&q[0])),
        ]
    };
    let dot = |p: &[Dyadic; 3], q: &[Dyadic; 3]| {
        p[0].mul(&q[0]).add(&p[1].mul(&q[1])).add(&p[2].mul(&q[2]))
    };
    let (u, v, w) = (vector(b), vector(c), vector(d));
    let n = cross(&u, &v);
    let power = dot(&n, &n)
        .mul(&dot(&w, &w))
        .sub(&dot(&u, &u).mul(&dot(&cross(&w, &v), &n)))
        .sub(&dot(&v, &v).mul(&dot(&cross(&u, &w), &n)));
    power.sign().expect("finite").flip()
}

/// Points on a tilted circle of decimal centre and radius, rounded: nearly
/// cocircular, with inexact coordinate differences, so the filter defers and
/// the exact path must keep every difference exact.
#[test]
fn in_diametral_sphere_matches_an_exact_dyadic_reference_near_cocircularity() {
    let mut state = 0x0BAD_C0DE_1234_5678u64;
    let mut unit = || (rng(&mut state) as f64) / (1u64 << 31) as f64;
    let mut deferred = 0usize;
    for _ in 0..5_000 {
        let centre = [unit() * 5.1, unit() * 2.3, unit() * 1.7];
        let r = 0.3 + unit();
        // An orthonormal pair spanning a tilted plane.
        let e1 = [2.0 / 3.0, 1.0 / 3.0, 2.0 / 3.0];
        let e2 = [1.0 / 2f64.sqrt(), 0.0, -1.0 / 2f64.sqrt()];
        let mut on = || {
            let t = unit() * std::f64::consts::TAU;
            let (cos, sin) = (r * t.cos(), r * t.sin());
            Point3::new(
                centre[0] + cos * e1[0] + sin * e2[0],
                centre[1] + cos * e1[1] + sin * e2[1],
                centre[2] + cos * e1[2] + sin * e2[2],
            )
        };
        let (a, b, c, d) = (on(), on(), on(), on());
        if !in_diametral_sphere_filter(a, b, c, d).is_certain() {
            deferred += 1;
        }
        assert_eq!(
            in_diametral_sphere(a, b, c, d).sign(),
            Some(dyadic_reference(a, b, c, d)),
            "{a:?} {b:?} {c:?} {d:?}"
        );
    }
    assert!(
        deferred > 1_000,
        "only {deferred} cases reached the exact path"
    );
}
