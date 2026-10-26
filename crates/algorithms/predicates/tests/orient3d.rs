//! Differential gate for `orient3d`.
//!
//! The oracle is an exact i128 determinant over integer coordinates, which is
//! independent of the expansion arithmetic under test. Integer inputs make the
//! coordinate differences exact, so any disagreement is a real predicate bug.

use axiolid_core::Point3;
use axiolid_guarantees::Sign;
use axiolid_predicates::{orient3d, orient3d_filter, two_sum};

fn p(x: i64, y: i64, z: i64) -> Point3 {
    Point3::new(x as f64, y as f64, z as f64)
}

/// Exact 3x3 determinant sign over integers.
///
/// Coordinates are bounded by the caller so the products stay inside i128.
fn oracle(a: [i64; 3], b: [i64; 3], c: [i64; 3], d: [i64; 3]) -> Sign {
    let m = |v: [i64; 3]| {
        [
            i128::from(v[0] - d[0]),
            i128::from(v[1] - d[1]),
            i128::from(v[2] - d[2]),
        ]
    };
    let (a, b, c) = (m(a), m(b), m(c));
    let det = a[2] * (b[0] * c[1] - c[0] * b[1])
        + b[2] * (c[0] * a[1] - a[0] * c[1])
        + c[2] * (a[0] * b[1] - b[0] * a[1]);
    match det.signum() {
        1 => Sign::Positive,
        -1 => Sign::Negative,
        _ => Sign::Zero,
    }
}

#[test]
fn obvious_orientations_are_correct() {
    // Unit tetrahedron: d below the xy-plane triangle.
    let sign = orient3d(p(0, 0, 0), p(1, 0, 0), p(0, 1, 0), p(0, 0, -1))
        .sign()
        .expect("certified");
    assert_eq!(sign, Sign::Positive);

    // Mirror the apex: the sign must flip.
    let flipped = orient3d(p(0, 0, 0), p(1, 0, 0), p(0, 1, 0), p(0, 0, 1))
        .sign()
        .expect("certified");
    assert_eq!(flipped, Sign::Negative);
}

#[test]
fn finite_extreme_coordinates_do_not_receive_a_false_exact_sign() {
    let magnitude = 1.0e200;
    let origin = Point3::new(0.0, 0.0, 0.0);
    let x = Point3::new(magnitude, 0.0, 0.0);
    let y = Point3::new(0.0, magnitude, 0.0);
    let z = Point3::new(0.0, 0.0, magnitude);

    assert!(!orient3d_filter(x, y, z, origin).is_certain());
    assert_eq!(orient3d(x, y, z, origin).sign(), Some(Sign::Positive));
}

#[test]
fn maximum_finite_cubic_determinant_fits_the_exact_accumulator() {
    let magnitude = f64::MAX;
    let origin = Point3::new(0.0, 0.0, 0.0);
    let x = Point3::new(magnitude, 0.0, 0.0);
    let y = Point3::new(0.0, magnitude, 0.0);
    let z = Point3::new(0.0, 0.0, magnitude);

    assert!(!orient3d_filter(x, y, z, origin).is_certain());
    assert_eq!(orient3d(x, y, z, origin).sign(), Some(Sign::Positive));
}

#[test]
fn minimum_positive_cubic_determinant_uses_exact_accumulator() {
    let magnitude = f64::from_bits(1);
    let origin = Point3::new(0.0, 0.0, 0.0);
    let x = Point3::new(magnitude, 0.0, 0.0);
    let y = Point3::new(0.0, magnitude, 0.0);
    let z = Point3::new(0.0, 0.0, magnitude);

    assert!(!orient3d_filter(x, y, z, origin).is_certain());
    assert_eq!(orient3d(x, y, z, origin).sign(), Some(Sign::Positive));
}

#[test]
fn finite_extreme_exponent_spread_uses_exact_dyadic_fallback() {
    let origin = Point3::new(0.0, 0.0, 0.0);
    let x = Point3::new(1.0e308, 0.0, 0.0);
    let y = Point3::new(0.0, 1.0, 0.0);
    let z = Point3::new(0.0, 0.0, 1.0e-300);

    assert!(!orient3d_filter(x, y, z, origin).is_certain());
    assert_eq!(orient3d(x, y, z, origin).sign(), Some(Sign::Positive));
}

#[test]
fn finite_coordinate_difference_overflow_uses_exact_dyadic_fallback() {
    let a = Point3::new(f64::MAX, 0.0, 0.0);
    let b = Point3::new(0.0, 1.0, 0.0);
    let c = Point3::new(0.0, 0.0, 1.0);
    let d = Point3::new(-f64::MAX, 0.0, 0.0);

    assert!(!orient3d_filter(a, b, c, d).is_certain());
    assert_eq!(orient3d(a, b, c, d).sign(), Some(Sign::Positive));
}

#[test]
fn huge_exact_coplanarity_stays_zero() {
    let magnitude = 1.0e200;
    let a = Point3::new(magnitude, 0.0, magnitude);
    let b = Point3::new(0.0, magnitude, magnitude);
    let c = Point3::new(magnitude, magnitude, magnitude);
    let d = Point3::new(0.5 * magnitude, 0.5 * magnitude, magnitude);

    assert!(!orient3d_filter(a, b, c, d).is_certain());
    assert_eq!(orient3d(a, b, c, d).sign(), Some(Sign::Zero));
}

#[test]
fn non_finite_coordinates_are_uncertain() {
    let origin = Point3::new(0.0, 0.0, 0.0);
    let x = Point3::new(f64::NAN, 0.0, 0.0);
    let y = Point3::new(0.0, 1.0, 0.0);
    let z = Point3::new(0.0, 0.0, 1.0);

    assert_eq!(orient3d(x, y, z, origin).sign(), None);
}

#[test]
fn exact_coplanarity_survives_inexact_coordinate_differences() {
    // Each coordinate satisfies d = a + b - c exactly as a binary rational,
    // but several plain `point - d` operations round. The exact path must carry
    // those subtraction tails instead of applying Sterbenz without its ratio
    // precondition.
    let a = Point3::new(
        -16_777_216.000_015_26,
        2.793_967_723_846_435_5e-9,
        33_554_431.999_999_996,
    );
    let b = Point3::new(
        -6_442_450_944.0,
        -1.907_348_632_814_235e-6,
        -0.000_122_070_312_5,
    );
    let c = Point3::new(
        68_719_476_736.062_48,
        3.337_860_107_421_875_4e-6,
        -9.536_743_164_062_5e-7,
    );
    let d = Point3::new(
        -75_178_704_896.062_5,
        -5.242_414_772_512_264e-6,
        33_554_431.999_878_88,
    );

    for (left, right) in [
        ((a.x, b.x), (c.x, d.x)),
        ((a.y, b.y), (c.y, d.y)),
        ((a.z, b.z), (c.z, d.z)),
    ] {
        assert_eq!(two_sum(left.0, left.1), two_sum(right.0, right.1));
    }
    assert!(!orient3d_filter(a, b, c, d).is_certain());
    assert_eq!(orient3d(a, b, c, d).sign(), Some(Sign::Zero));
}

#[test]
fn exact_coplanarity_is_reported_as_zero() {
    // Four points on z = 0. A predicate that merely approximates would return
    // a tiny non-zero determinant here and corrupt every downstream decision.
    let sign = orient3d(p(0, 0, 0), p(5, 0, 0), p(0, 7, 0), p(3, 4, 0))
        .sign()
        .expect("certified");
    assert_eq!(sign, Sign::Zero);
}

/// Deterministic small integers, so the i128 oracle stays exact.
fn coords(state: &mut u64, bound: i64) -> [i64; 3] {
    let mut next = || {
        *state ^= *state << 13;
        *state ^= *state >> 7;
        *state ^= *state << 17;
        ((*state >> 33) as i64) % (2 * bound) - bound
    };
    [next(), next(), next()]
}

/// A point exactly on the plane through `a`, `b`, `c`.
///
/// Built as an integer affine combination `a + i*(b-a) + j*(c-a)`, which lies
/// on the plane by construction for any integers i, j -- no division, so no
/// rounding can push it off.
fn coplanar_with(a: [i64; 3], b: [i64; 3], c: [i64; 3], state: &mut u64) -> [i64; 3] {
    *state ^= *state << 13;
    *state ^= *state >> 7;
    *state ^= *state << 17;
    let i = ((*state >> 33) as i64) % 5 - 2;
    let j = ((*state >> 11) as i64) % 5 - 2;
    [
        a[0] + i * (b[0] - a[0]) + j * (c[0] - a[0]),
        a[1] + i * (b[1] - a[1]) + j * (c[1] - a[1]),
        a[2] + i * (b[2] - a[2]) + j * (c[2] - a[2]),
    ]
}

/// The predicate must agree with the independent oracle on every input.
#[test]
fn certified_signs_agree_with_an_independent_exact_oracle() {
    let mut state = 0x9E37_79B9_7F4A_7C15;
    let mut compared = 0usize;
    let mut zeros = 0usize;
    for _ in 0..20_000 {
        let (a, b, c) = (
            coords(&mut state, 50),
            coords(&mut state, 50),
            coords(&mut state, 50),
        );
        // Half the cases are free; half are forced onto the plane through
        // a, b, c so exact zeros actually occur. Random points in a 100^3 box
        // are essentially never coplanar, so an unbiased sweep would only ever
        // exercise the easy path.
        let d = if compared.is_multiple_of(2) {
            coords(&mut state, 50)
        } else {
            coplanar_with(a, b, c, &mut state)
        };
        let got = orient3d(
            p(a[0], a[1], a[2]),
            p(b[0], b[1], b[2]),
            p(c[0], c[1], c[2]),
            p(d[0], d[1], d[2]),
        )
        .sign()
        .expect("orient3d always certifies");
        let want = oracle(a, b, c, d);
        assert_eq!(got, want, "disagreement on {a:?} {b:?} {c:?} {d:?}");
        if want == Sign::Zero {
            zeros += 1;
        }
        compared += 1;
    }
    assert_eq!(compared, 20_000);
    // Small coordinates produce genuine coplanarities; if none appeared the
    // sweep would only be testing the easy path.
    assert!(zeros > 0, "no exactly coplanar case was generated");
}

/// Nearly-coplanar inputs must still be decided correctly.
///
/// This is the case the filter cannot settle, so it exercises the exact path
/// specifically -- and it asserts the escalation actually happened, or the
/// test would silently degrade into another easy-path check.
#[test]
fn near_coplanar_inputs_escalate_and_stay_correct() {
    let mut escalated = 0usize;
    let mut checked = 0usize;
    for k in 1..400i64 {
        // Three points spanning z = 0, and a fourth a hair off the plane.
        let (a, b, c) = ([0, 0, 0], [k, 0, 0], [0, k, 0]);
        for offset in [-1i64, 0, 1] {
            let d = [k / 3, k / 4, offset];
            let filtered = orient3d_filter(
                p(a[0], a[1], a[2]),
                p(b[0], b[1], b[2]),
                p(c[0], c[1], c[2]),
                p(d[0], d[1], d[2]),
            );
            if !filtered.is_certain() {
                escalated += 1;
            }
            let got = orient3d(
                p(a[0], a[1], a[2]),
                p(b[0], b[1], b[2]),
                p(c[0], c[1], c[2]),
                p(d[0], d[1], d[2]),
            )
            .sign()
            .expect("certified");
            assert_eq!(got, oracle(a, b, c, d), "k={k} offset={offset}");
            checked += 1;
        }
    }
    assert!(checked > 1_000);
    assert!(
        escalated > 0,
        "no input reached the exact path; this test proves nothing"
    );
}
