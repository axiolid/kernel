//! Differential gates for `incircle` and `insphere`.
//!
//! The oracles are exact i128 determinants over integer coordinates, entirely
//! independent of the expansion arithmetic under test. Coordinate bounds are
//! chosen so the lifted determinants cannot overflow i128: `incircle` is
//! degree 4 and `insphere` degree 6 in the coordinates, so the bounds differ.

use axiolid_core::{Point2, Point3};
use axiolid_guarantees::Sign;
use axiolid_predicates::{incircle, incircle_filter, insphere, insphere_filter};

fn sign_of(v: i128) -> Sign {
    match v.signum() {
        1 => Sign::Positive,
        -1 => Sign::Negative,
        _ => Sign::Zero,
    }
}

/// Exact `incircle` over integers: a 3x3 lifted determinant.
fn incircle_oracle(a: [i64; 2], b: [i64; 2], c: [i64; 2], d: [i64; 2]) -> Sign {
    let m = |v: [i64; 2]| {
        let (x, y) = (i128::from(v[0] - d[0]), i128::from(v[1] - d[1]));
        (x, y, x * x + y * y)
    };
    let (a, b, c) = (m(a), m(b), m(c));
    sign_of(
        a.2 * (b.0 * c.1 - c.0 * b.1)
            + b.2 * (c.0 * a.1 - a.0 * c.1)
            + c.2 * (a.0 * b.1 - b.0 * a.1),
    )
}

/// Exact `insphere` over integers: a 4x4 lifted determinant, expanded by
/// cofactors along the lifted column.
fn insphere_oracle(a: [i64; 3], b: [i64; 3], c: [i64; 3], d: [i64; 3], e: [i64; 3]) -> Sign {
    let m = |v: [i64; 3]| {
        let (x, y, z) = (
            i128::from(v[0] - e[0]),
            i128::from(v[1] - e[1]),
            i128::from(v[2] - e[2]),
        );
        (x, y, z, x * x + y * y + z * z)
    };
    let (a, b, c, d) = (m(a), m(b), m(c), m(d));
    let det3 =
        |p: (i128, i128, i128, i128), q: (i128, i128, i128, i128), r: (i128, i128, i128, i128)| {
            p.0 * (q.1 * r.2 - r.1 * q.2) - p.1 * (q.0 * r.2 - r.0 * q.2)
                + p.2 * (q.0 * r.1 - r.0 * q.1)
        };
    sign_of(-a.3 * det3(b, c, d) + b.3 * det3(a, c, d) - c.3 * det3(a, b, d) + d.3 * det3(a, b, c))
}

fn rng(state: &mut u64) -> i64 {
    *state ^= *state << 13;
    *state ^= *state >> 7;
    *state ^= *state << 17;
    (*state >> 33) as i64
}

/// `incircle` must match the independent oracle, including cocircular cases.
#[test]
fn incircle_agrees_with_an_independent_exact_oracle() {
    let mut state = 0x2545_F491_4F6C_DD1D;
    let mut zeros = 0usize;
    for i in 0..20_000 {
        let bound = 60;
        let pt = |s: &mut u64| [rng(s) % (2 * bound) - bound, rng(s) % (2 * bound) - bound];
        let (a, b, c) = (pt(&mut state), pt(&mut state), pt(&mut state));
        // Half the cases are placed exactly on a lattice circle so genuine
        // zeros occur; random points are essentially never cocircular.
        let d = if i % 2 == 0 {
            pt(&mut state)
        } else {
            // Reflect `a` across the perpendicular bisector midpoint of b, c:
            // integer construction, so it stays exact.
            [b[0] + c[0] - a[0], b[1] + c[1] - a[1]]
        };
        let got = incircle(
            Point2::new(a[0] as f64, a[1] as f64),
            Point2::new(b[0] as f64, b[1] as f64),
            Point2::new(c[0] as f64, c[1] as f64),
            Point2::new(d[0] as f64, d[1] as f64),
        )
        .sign()
        .expect("incircle always certifies");
        let want = incircle_oracle(a, b, c, d);
        assert_eq!(got, want, "disagreement on {a:?} {b:?} {c:?} {d:?}");
        if want == Sign::Zero {
            zeros += 1;
        }
    }
    assert!(zeros > 0, "no exactly cocircular case was generated");
}

/// `insphere` must match the independent oracle, including cospherical cases.
///
/// The bound is tighter than `incircle`'s: the 4x4 lifted determinant is
/// degree 6 in the coordinates, so 2^6 * bound^6 must stay inside i128.
#[test]
fn insphere_agrees_with_an_independent_exact_oracle() {
    let mut state = 0x9E37_79B9_7F4A_7C15;
    let mut zeros = 0usize;
    let bound = 40i64;
    for i in 0..20_000 {
        let pt = |s: &mut u64| {
            [
                rng(s) % (2 * bound) - bound,
                rng(s) % (2 * bound) - bound,
                rng(s) % (2 * bound) - bound,
            ]
        };
        // Build a, b, c, d ON a lattice sphere of radius r about a random
        // centre, so a fifth point on that sphere is exactly cospherical by
        // construction. Random points are essentially never cospherical, so an
        // unbiased sweep would only ever exercise the easy path.
        let centre = pt(&mut state);
        let r = 5i64;
        let on_sphere = |k: usize| -> [i64; 3] {
            // (+-r,0,0), (0,+-r,0), (0,0,+-r): all exactly at distance r.
            let (dx, dy, dz) = match k % 6 {
                0 => (r, 0, 0),
                1 => (-r, 0, 0),
                2 => (0, r, 0),
                3 => (0, -r, 0),
                4 => (0, 0, r),
                _ => (0, 0, -r),
            };
            [centre[0] + dx, centre[1] + dy, centre[2] + dz]
        };
        let (a, b, c, d) = if i % 2 == 1 {
            (on_sphere(0), on_sphere(2), on_sphere(4), on_sphere(1))
        } else {
            (
                pt(&mut state),
                pt(&mut state),
                pt(&mut state),
                pt(&mut state),
            )
        };
        // Points on an axis-aligned lattice sphere of radius^2 = r2, which has
        // exact integer solutions, so genuine cospherical cases occur.
        // On the odd iterations a, b, c, d sit on the lattice sphere, so
        // picking a fifth lattice point on it gives an exact zero.
        let e = if i % 2 == 1 {
            on_sphere(3)
        } else {
            pt(&mut state)
        };
        let f = |v: [i64; 3]| Point3::new(v[0] as f64, v[1] as f64, v[2] as f64);
        let got = insphere(f(a), f(b), f(c), f(d), f(e))
            .sign()
            .expect("insphere always certifies");
        let want = insphere_oracle(a, b, c, d, e);
        assert_eq!(got, want, "disagreement on {a:?} {b:?} {c:?} {d:?} {e:?}");
        if want == Sign::Zero {
            zeros += 1;
        }
    }
    assert!(zeros > 0, "no exactly cospherical case was generated");
}

/// Both filters must actually defer on degenerate input, or the exact paths
/// above are never exercised and these tests prove less than they appear to.
#[test]
fn the_filters_defer_on_exactly_degenerate_input() {
    // Four corners of a unit square are exactly cocircular.
    let square = [
        Point2::new(0.0, 0.0),
        Point2::new(1.0, 0.0),
        Point2::new(1.0, 1.0),
        Point2::new(0.0, 1.0),
    ];
    assert!(
        !incircle_filter(square[0], square[1], square[2], square[3]).is_certain(),
        "an exactly cocircular case must not be settled by the filter"
    );

    // Six vertices of an octahedron are exactly cospherical.
    let o = |x: f64, y: f64, z: f64| Point3::new(x, y, z);
    assert!(
        !insphere_filter(
            o(1.0, 0.0, 0.0),
            o(-1.0, 0.0, 0.0),
            o(0.0, 1.0, 0.0),
            o(0.0, 0.0, 1.0),
            o(0.0, -1.0, 0.0)
        )
        .is_certain(),
        "an exactly cospherical case must not be settled by the filter"
    );
}

/// The exact paths must recover a definite sign just off the boundary.
///
/// Mirrors the `orient3d` gate: replacing `incircle_exact`/`insphere_exact`
/// with a constant `Sign::Zero` passed every other test in this file, because
/// all of their degenerate cases are exactly on the ball and answer zero too.
/// Perturbing by one ULP forces the filter to defer AND demands a real answer.
#[test]
fn near_degenerate_delaunay_cases_recover_a_definite_sign() {
    let ulp = |v: f64| {
        let bits = v.to_bits();
        f64::from_bits(if v >= 0.0 { bits + 1 } else { bits - 1 })
    };

    // Unit square: exactly cocircular. Nudge the fourth point off the circle.
    let mut deferred = 0usize;
    let mut definite = 0usize;
    for k in 1..2_000i64 {
        let r = k as f64;
        let (a, b, c) = (Point2::new(-r, -r), Point2::new(r, -r), Point2::new(r, r));
        let d = Point2::new(-r, ulp(r));
        if !incircle_filter(a, b, c, d).is_certain() {
            deferred += 1;
        }
        if incircle(a, b, c, d).sign().expect("certified") != Sign::Zero {
            definite += 1;
        }
    }
    assert!(deferred > 1_000, "incircle: only {deferred} deferred");
    assert!(definite > 1_000, "incircle: only {definite} definite");

    // Octahedron: exactly cospherical. Nudge the fifth vertex off the sphere.
    let mut deferred = 0usize;
    let mut definite = 0usize;
    for k in 1..2_000i64 {
        let r = k as f64;
        let o = |x: f64, y: f64, z: f64| Point3::new(x, y, z);
        let (a, b, c, d) = (
            o(r, 0.0, 0.0),
            o(-r, 0.0, 0.0),
            o(0.0, r, 0.0),
            o(0.0, 0.0, r),
        );
        let e = o(0.0, -ulp(r), 0.0);
        if !insphere_filter(a, b, c, d, e).is_certain() {
            deferred += 1;
        }
        if insphere(a, b, c, d, e).sign().expect("certified") != Sign::Zero {
            definite += 1;
        }
    }
    assert!(deferred > 1_000, "insphere: only {deferred} deferred");
    assert!(definite > 1_000, "insphere: only {definite} definite");
}

/// Lattice points of the sphere of radius `m^2 + 2` from the Pythagorean
/// quadruple `(m^2, 2, 2m)`: every sign and permutation, plus the six axis
/// points. All exactly cospherical, with coordinate differences spanning
/// many bits.
fn quadruple_sphere(m: f64) -> Vec<Point3> {
    let r = m * m + 2.0;
    let base = [m * m, 2.0, 2.0 * m];
    let mut points = Vec::new();
    for perm in [
        [0, 1, 2],
        [0, 2, 1],
        [1, 0, 2],
        [1, 2, 0],
        [2, 0, 1],
        [2, 1, 0],
    ] {
        for signs in 0..8u32 {
            let s = |k: u32| if (signs >> k) & 1 == 1 { -1.0 } else { 1.0 };
            points.push(Point3::new(
                s(0) * base[perm[0]],
                s(1) * base[perm[1]],
                s(2) * base[perm[2]],
            ));
        }
    }
    for axis in 0..3 {
        for sign in [-1.0, 1.0] {
            let mut c = [0.0; 3];
            c[axis] = sign * r;
            points.push(Point3::new(c[0], c[1], c[2]));
        }
    }
    points
}

/// The `insphere` filter must never certify a sign for exactly cospherical
/// points (#126). It used to bound its error by the rounded 3x3 minors
/// instead of by the absolute values of the elementary products, and
/// certified "outside" for five points of one sphere whose tetrahedron was
/// thin.
#[test]
fn insphere_never_certifies_a_sign_for_cospherical_points() {
    let (m, k) = (16_785_409.0, 8194.0);
    let thin = [
        Point3::new(-2.0, -k, -m),
        Point3::new(-k, 2.0, -m),
        Point3::new(2.0, -k, -m),
        Point3::new(-k, -m, 2.0),
        Point3::new(2.0, -m, -k),
    ];
    assert!(!insphere_filter(thin[0], thin[1], thin[2], thin[3], thin[4]).is_certain());
    assert_eq!(
        insphere(thin[0], thin[1], thin[2], thin[3], thin[4]).sign(),
        Some(Sign::Zero)
    );

    for m in [65.0, 1025.0, 4097.0] {
        let points = quadruple_sphere(m);
        let mut state = 0x5DEE_CE66_D1CE_4E5Bu64;
        for _ in 0..50_000 {
            let mut pick = || points[(rng(&mut state) as usize) % points.len()];
            let (a, b, c, d, e) = (pick(), pick(), pick(), pick(), pick());
            assert_eq!(
                insphere(a, b, c, d, e).sign(),
                Some(Sign::Zero),
                "{a:?} {b:?} {c:?} {d:?} {e:?}"
            );
        }
    }
}
