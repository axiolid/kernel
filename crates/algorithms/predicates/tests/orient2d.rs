//! The predicate is only worth building if it beats naive f64 on the cases
//! that break naive f64. These tests are chosen to fail against a plain
//! determinant implementation.

use axiolid_core::Point2;
use axiolid_guarantees::{Certified, Precision, Sign};
use axiolid_predicates::{orient2d, orient2d_filter};

fn p(x: f64, y: f64) -> Point2 {
    Point2::new(x, y)
}

fn sign_of(result: Certified) -> Sign {
    match result {
        Certified::Certain { sign, .. } => sign,
        other => panic!("orient2d must always certify, got {other:?}"),
    }
}

#[test]
fn obvious_orientations_are_correct() {
    let (a, b) = (p(0.0, 0.0), p(1.0, 0.0));
    assert_eq!(sign_of(orient2d(a, b, p(0.0, 1.0))), Sign::Positive);
    assert_eq!(sign_of(orient2d(a, b, p(0.0, -1.0))), Sign::Negative);
    assert_eq!(sign_of(orient2d(a, b, p(2.0, 0.0))), Sign::Zero);
}

/// Exactly collinear points must return Zero, not a tolerance-sized guess.
/// A tolerance-based predicate cannot distinguish "on the line" from "nearly
/// on the line", which is how sliver triangles and non-manifold edges appear.
#[test]
fn exact_collinearity_is_reported_as_zero() {
    let a = p(0.0, 0.0);
    let b = p(1.0, 1.0);
    for step in 1..50 {
        let t = f64::from(step);
        assert_eq!(
            sign_of(orient2d(a, b, p(t, t))),
            Sign::Zero,
            "points on y = x are collinear at t = {t}"
        );
    }
}

/// The case that motivates the whole cascade.
///
/// These are not hand-waved "nearly degenerate" inputs: each triple was found
/// by search and cross-checked against an exact rational determinant. In every
/// one the naive f64 determinant evaluates to exactly 0.0 -- reporting
/// "collinear" -- while the true determinant is a definite non-zero sign.
///
/// A predicate that reports collinear here fuses distinct vertices, drops
/// triangles, and produces non-manifold output. This is the bug the exact path
/// exists to prevent, so the test asserts against the naive result directly.
#[test]
fn the_exact_path_beats_naive_f64_where_naive_f64_is_wrong() {
    // (a, b, c, true sign as computed by exact rational arithmetic)
    let cases = [
        (
            p(-0.8362899784084603, -0.3995017629087494),
            p(-0.00976728088948886, -0.3130486200832534),
            p(1.3309231823978052, -0.17281423274545576),
            Sign::Negative,
        ),
        (
            p(0.9603496949851642, -0.7638684434900758),
            p(-0.1637543564295456, 0.5142818591304987),
            p(-1.1529224600475325, 1.6390047078657146),
            Sign::Negative,
        ),
        (
            p(-0.9918127932298721, -0.16210699774934412),
            p(-0.26149285421054924, 0.1326824474127839),
            p(1.8438331624214392, 0.9824851916206433),
            Sign::Positive,
        ),
    ];

    for (a, b, c, truth) in cases {
        let naive = (a.x - c.x) * (b.y - c.y) - (a.y - c.y) * (b.x - c.x);
        assert_eq!(
            naive, 0.0,
            "precondition: naive f64 must claim collinear for {a:?} {b:?} {c:?}"
        );
        assert_eq!(
            sign_of(orient2d(a, b, c)),
            truth,
            "exact arithmetic must overrule the naive determinant"
        );
        // The filter must admit it cannot decide, rather than trusting a zero.
        assert!(matches!(
            orient2d_filter(a, b, c),
            Certified::Uncertain { .. }
        ));
    }
}

/// Antisymmetry: swapping two points must flip the sign, for every input.
/// A predicate that violates this produces contradictory topology depending on
/// vertex order, which is unfixable downstream.
#[test]
fn the_predicate_is_antisymmetric_including_near_degenerate_inputs() {
    let mut state = 0x9E37_79B9_7F4A_7C15_u64;
    let mut next = || {
        state ^= state << 13;
        state ^= state >> 7;
        state ^= state << 17;
        f64::from(((state >> 40) as u32) as i32) / 8.0
    };

    for _ in 0..500 {
        let a = p(next(), next());
        let b = p(next(), next());
        // Place c on the line AB, then nudge it by a sub-ulp amount so the
        // filter is forced into its uncertain band.
        let c = p(a.x + (b.x - a.x) * 3.0, a.y + (b.y - a.y) * 3.0);

        let forward = sign_of(orient2d(a, b, c));
        let swapped = sign_of(orient2d(b, a, c));
        assert_eq!(
            forward,
            swapped.flip(),
            "orient2d({a:?}, {b:?}, {c:?}) must be antisymmetric"
        );
    }
}

/// The filter must be honest about not knowing, and the full predicate must
/// still produce a proven answer in exactly those cases.
#[test]
fn the_filter_defers_where_the_full_predicate_still_decides() {
    let a = p(0.0, 0.0);
    let b = p(1.0, 1.0);
    // Collinear: the determinant is exactly zero, so it can never exceed its
    // own error bound and the filter must report Uncertain.
    let c = p(3.0, 3.0);

    assert!(
        matches!(orient2d_filter(a, b, c), Certified::Uncertain { .. }),
        "a zero determinant cannot be certified by a filter"
    );
    // The cascade still resolves it, and says so with Exact precision.
    assert!(matches!(
        orient2d(a, b, c),
        Certified::Certain {
            sign: Sign::Zero,
            precision: Precision::Exact
        }
    ));
}

/// Well-separated inputs must be settled by the cheap path, or the cascade is
/// pointless: the exact path would run every time.
#[test]
fn the_filter_settles_ordinary_inputs_without_escalating() {
    let (a, b) = (p(0.0, 0.0), p(1.0, 0.0));
    assert!(matches!(
        orient2d_filter(a, b, p(0.0, 1.0)),
        Certified::Certain {
            sign: Sign::Positive,
            precision: Precision::F64
        }
    ));
}

/// Differential gate: the certified predicate must agree with an independent
/// exact oracle on every input, including the near-degenerate band.
///
/// The oracle is i128 integer arithmetic over small integer coordinates. For
/// such inputs the determinant is exactly representable in i128, so the oracle
/// is correct by construction and shares no code with the crate under test.
/// Coordinates are scaled to force the f64 filter into its uncertain band.
#[test]
fn certified_signs_agree_with_an_independent_exact_oracle() {
    fn oracle(a: (i64, i64), b: (i64, i64), c: (i64, i64)) -> Sign {
        let d = (i128::from(a.0) - i128::from(c.0)) * (i128::from(b.1) - i128::from(c.1))
            - (i128::from(a.1) - i128::from(c.1)) * (i128::from(b.0) - i128::from(c.0));
        match d.cmp(&0) {
            core::cmp::Ordering::Greater => Sign::Positive,
            core::cmp::Ordering::Less => Sign::Negative,
            core::cmp::Ordering::Equal => Sign::Zero,
        }
    }

    let mut state = 0xDEAD_BEEF_CAFE_1234_u64;
    let mut next = |range: i64| {
        state ^= state << 13;
        state ^= state >> 7;
        state ^= state << 17;
        (state >> 33) as i64 % range - range / 2
    };

    let mut escalations = 0_u32;
    let total = 20_000;

    for _ in 0..total {
        // Small coordinate range => many exactly-collinear and near-collinear
        // triples, which is precisely where a filter-only predicate fails.
        let a = (next(64), next(64));
        let b = (next(64), next(64));
        let c = (next(64), next(64));

        let fa = p(a.0 as f64, a.1 as f64);
        let fb = p(b.0 as f64, b.1 as f64);
        let fc = p(c.0 as f64, c.1 as f64);

        assert_eq!(
            sign_of(orient2d(fa, fb, fc)),
            oracle(a, b, c),
            "disagreement on a={a:?} b={b:?} c={c:?}"
        );

        if matches!(orient2d_filter(fa, fb, fc), Certified::Uncertain { .. }) {
            escalations += 1;
        }
    }

    // The exact path must actually be reached, or this test proves nothing
    // about it. Degenerate triples are common at this coordinate scale.
    assert!(
        escalations > 0,
        "no input escalated; the differential test never exercised exact arithmetic"
    );
}

/// Differential gate for the *correction* terms, which the integer-coordinate
/// test above cannot reach.
///
/// When all coordinates share a magnitude, `a.x - c.x` is exact and every error
/// term in the expansion is zero -- so that test still passes if the correction
/// terms are deleted. Real BIM data is not like that: a national-grid origin and
/// a millimetre detail differ by ~10 orders of magnitude, the subtraction
/// rounds, and the discarded bits decide the sign.
///
/// Every coordinate here is an exact dyadic rational (`mantissa * 2^exponent`)
/// with a wide exponent spread. The oracle converts each f64 to an exact
/// `Fraction`-style integer pair and evaluates the determinant in i128 with a
/// common denominator, so it shares no code with the implementation.
#[test]
fn certified_signs_survive_mixed_magnitude_coordinates() {
    /// Exact value of a dyadic f64 as `mantissa * 2^exponent`, mantissa in i128.
    fn decompose(value: f64) -> (i128, i32) {
        if value == 0.0 {
            return (0, 0);
        }
        let bits = value.to_bits();
        let raw_exponent = ((bits >> 52) & 0x7FF) as i32;
        let fraction = bits & 0x000F_FFFF_FFFF_FFFF;
        let (mantissa, exponent) = if raw_exponent == 0 {
            (fraction, -1074)
        } else {
            (fraction | 0x0010_0000_0000_0000, raw_exponent - 1075)
        };
        let mut signed = if value < 0.0 {
            -(mantissa as i128)
        } else {
            mantissa as i128
        };
        // Normalise: strip trailing zero bits so the common-denominator shift
        // below stays far inside i128 instead of overflowing on a wide spread.
        let mut exponent = exponent;
        while signed != 0 && signed % 2 == 0 {
            signed /= 2;
            exponent += 1;
        }
        (signed, exponent)
    }

    /// Exact sign of the orientation determinant, computed in i128 after
    /// rescaling every coordinate to a common power-of-two denominator.
    fn oracle(a: (f64, f64), b: (f64, f64), c: (f64, f64)) -> Sign {
        let parts = [a.0, a.1, b.0, b.1, c.0, c.1].map(decompose);
        let min_exponent = parts.iter().map(|&(_, e)| e).min().expect("non-empty");
        // Rescaling to the smallest exponent makes every value an exact
        // integer; the shift is bounded because the generator caps the spread.
        let max_shift = parts
            .iter()
            .map(|&(_, e)| e - min_exponent)
            .max()
            .expect("non-empty");
        // 12-bit mantissas plus a bounded spread must leave room for the
        // determinant's products; assert it rather than overflowing silently.
        assert!(
            max_shift < 40,
            "exponent spread {max_shift} too wide for i128"
        );
        let scaled = parts.map(|(m, e)| m << (e - min_exponent));
        let [ax, ay, bx, by, cx, cy] = scaled;
        let d = (ax - cx) * (by - cy) - (ay - cy) * (bx - cx);
        match d.cmp(&0) {
            core::cmp::Ordering::Greater => Sign::Positive,
            core::cmp::Ordering::Less => Sign::Negative,
            core::cmp::Ordering::Equal => Sign::Zero,
        }
    }

    let mut state = 0x1234_5678_9ABC_DEF0_u64;
    let mut rng = || {
        state ^= state << 13;
        state ^= state >> 7;
        state ^= state << 17;
        state
    };

    let mut escalations = 0_u32;
    let mut rounding_losses = 0_u32;

    for _ in 0..20_000 {
        // Small mantissas with a wide exponent spread: the products stay
        // exactly representable while the subtractions genuinely round.
        let coordinate = |rng: &mut dyn FnMut() -> u64| {
            let r = rng();
            // Full 53-bit mantissas: this is what makes `a.x - c.x` round when
            // the exponents differ, which is the whole point of this test.
            let mantissa = ((r >> 11) as i64) - (1 << 52);
            let exponent = ((r % 9) as i32) - 4;
            (mantissa as f64) * 2.0_f64.powi(exponent)
        };

        let ax = coordinate(&mut rng);
        let ay = coordinate(&mut rng);
        let bx = coordinate(&mut rng);
        let by = coordinate(&mut rng);
        // Place c near the line AB so the filter is pushed into its uncertain
        // band; a random third point is essentially never near-degenerate.
        let k = 2.0 + (rng() % 3) as f64;
        let cx = ax + (bx - ax) * k;
        let cy = ay + (by - ay) * k;

        let (fa, fb, fc) = (p(ax, ay), p(bx, by), p(cx, cy));

        // Premise check: the f64 subtraction must lose bits somewhere, or the
        // correction terms remain unreachable and this test proves nothing.
        let (dm, de) = decompose(ax - cx);
        let (am, ae) = decompose(ax);
        let (cm, ce) = decompose(cx);
        let exact_difference =
            (am << (ae - ae.min(ce).min(de))) - (cm << (ce - ae.min(ce).min(de)));
        if (dm << (de - ae.min(ce).min(de))) != exact_difference {
            rounding_losses += 1;
        }

        assert_eq!(
            sign_of(orient2d(fa, fb, fc)),
            oracle((ax, ay), (bx, by), (cx, cy)),
            "disagreement on a=({ax}, {ay}) b=({bx}, {by}) c=({cx}, {cy})"
        );

        if matches!(orient2d_filter(fa, fb, fc), Certified::Uncertain { .. }) {
            escalations += 1;
        }
    }

    assert!(escalations > 0, "the exact path was never exercised");
    assert!(
        rounding_losses > 0,
        "no coordinate difference rounded; correction terms stayed unreachable"
    );
}

/// Reference: the determinant over exact dyadics, which has no exponent range.
fn orient2d_reference(a: Point2, b: Point2, c: Point2) -> Sign {
    use axiolid_exact::{Arith, Dyadic};
    let d = Dyadic::from_f64;
    let left = d(a.x).sub(&d(c.x)).mul(&d(b.y).sub(&d(c.y)));
    let right = d(a.y).sub(&d(c.y)).mul(&d(b.x).sub(&d(c.x)));
    left.sub(&right).sign().unwrap()
}

fn assert_matches_reference(a: Point2, b: Point2, c: Point2) {
    assert_eq!(
        sign_of(orient2d(a, b, c)),
        orient2d_reference(a, b, c),
        "a = {a:?}, b = {b:?}, c = {c:?}"
    );
}

/// Beyond about 1e154 the products overflow (#404): the sign must still be
/// certain and right for collinear, left and right triples.
#[test]
fn overflowing_products_keep_an_exact_sign() {
    for scale in [1.0e154, 1.0e160, 1.0e200, 1.0e300, f64::MAX / 4.0] {
        let (a, b) = (p(0.0, 0.0), p(scale, scale));
        // Exactly on the line y = x.
        assert_eq!(
            sign_of(orient2d(a, b, p(scale * 0.5, scale * 0.5))),
            Sign::Zero
        );
        assert_eq!(sign_of(orient2d(a, b, p(scale, scale))), Sign::Zero);
        assert_eq!(
            sign_of(orient2d(a, b, p(scale, scale * 0.5))),
            Sign::Negative
        );
        assert_eq!(
            sign_of(orient2d(a, b, p(scale * 0.5, scale))),
            Sign::Positive
        );
        // The issue's coplanar triangles, scaled.
        let t = |x: f64, y: f64| p(x * scale, y * scale);
        assert_matches_reference(t(0.0, 0.0), t(2.0, 0.0), t(0.5, 0.5));
        assert_matches_reference(t(0.5, 0.5), t(3.0, 0.5), t(2.0, 0.0));
        assert_matches_reference(t(0.5, 0.5), t(3.0, 0.5), t(0.5, 3.0));
    }
}

/// One ulp off the line at 1e300 decides the sign; so does one ulp on the
/// other side, and the line itself is Zero.
#[test]
fn a_one_ulp_offset_at_huge_scale_is_decided() {
    let a = p(1.0e300, 1.0e300);
    let b = p(2.0e300, 2.0e300);
    let on: f64 = 1.5e300;
    let up = f64::from_bits(on.to_bits() + 1);
    let down = f64::from_bits(on.to_bits() - 1);
    assert_eq!(sign_of(orient2d(a, b, p(on, on))), Sign::Zero);
    assert_matches_reference(a, b, p(on, up));
    assert_matches_reference(a, b, p(on, down));
    assert_ne!(sign_of(orient2d(a, b, p(on, up))), Sign::Zero);
    assert_ne!(
        sign_of(orient2d(a, b, p(on, up))),
        sign_of(orient2d(a, b, p(on, down)))
    );
}

/// Coordinates from 5e-324 (subnormal) to 1e300 in one triple: differences
/// and products overflow or underflow, the sign must not.
#[test]
fn mixed_scale_and_underflowing_inputs_keep_an_exact_sign() {
    let scales = [
        0.0, 5.0e-324, 1.0e-320, 1.0e-300, 1.0e-160, 1.0e-155, 1.0e-100, 1.0, 3.0e50, 1.0e150,
        1.0e300,
    ];
    let mut state = 0x9e37_79b9_7f4a_7c15_u64;
    let mut next = move || {
        state ^= state << 13;
        state ^= state >> 7;
        state ^= state << 17;
        state
    };
    let mut coordinate = || {
        let magnitude = scales[(next() % scales.len() as u64) as usize];
        let factor = 1.0 + (next() % 7) as f64 / 8.0;
        if next() % 2 == 0 {
            magnitude * factor
        } else {
            -magnitude * factor
        }
    };
    for _ in 0..4000 {
        let a = p(coordinate(), coordinate());
        let b = p(coordinate(), coordinate());
        let c = p(coordinate(), coordinate());
        assert_matches_reference(a, b, c);
    }
}

/// The differences are in range but a coordinate is tiny against the others,
/// so the error terms of the expansion path would underflow.
#[test]
fn tiny_coordinates_beside_large_differences_keep_an_exact_sign() {
    let a = p(1.0e-323, -1.625e-100);
    let b = p(1.625e-100, -1.625e-100);
    let c = p(-5.2499999999999995e50, -1.375);
    assert_matches_reference(a, b, c);
    assert_eq!(sign_of(orient2d(a, b, c)), Sign::Negative);
}

/// Everything underflows: points a few subnormal ulps apart.
#[test]
fn subnormal_inputs_keep_an_exact_sign() {
    let u = 5.0e-324;
    assert_eq!(
        sign_of(orient2d(p(0.0, 0.0), p(u, u), p(2.0 * u, 2.0 * u))),
        Sign::Zero
    );
    assert_eq!(
        sign_of(orient2d(p(0.0, 0.0), p(2.0 * u, 0.0), p(0.0, u))),
        Sign::Positive
    );
    assert_eq!(
        sign_of(orient2d(p(0.0, 0.0), p(0.0, u), p(2.0 * u, 0.0))),
        Sign::Negative
    );
    // Normal but tiny: the product is below the smallest subnormal.
    assert_eq!(
        sign_of(orient2d(p(0.0, 0.0), p(1.0e-200, 0.0), p(0.0, 1.0e-200))),
        Sign::Positive
    );
}

#[test]
fn a_non_finite_coordinate_is_uncertain() {
    for bad in [f64::NAN, f64::INFINITY] {
        assert!(matches!(
            orient2d(p(0.0, 0.0), p(1.0, 0.0), p(bad, 1.0)),
            Certified::Uncertain { .. }
        ));
    }
}

/// Differences below the filter's lower range bound are declined, not
/// bounded: without that bound the filter proved a wrong Negative here.
#[test]
fn a_triple_below_the_filter_range_keeps_its_true_sign() {
    let a = p(9.047035581504164e-156, -1.9220684987552693e-157);
    let b = p(-6.199974410654496e-156, -6.004478465470848e-156);
    let c = p(3.1950227522507125e-155, 8.538657266674845e-156);
    assert_eq!(orient2d_reference(a, b, c), Sign::Positive);
    assert!(matches!(
        orient2d_filter(a, b, c),
        Certified::Uncertain { .. }
    ));
    assert_matches_reference(a, b, c);
}

/// Near-collinear triples with non-dyadic coordinates around 1e-155, where
/// the filter's products approach the underflow threshold.
#[test]
fn near_collinear_sweep_at_the_underflow_scale_keeps_an_exact_sign() {
    let mut state = 0x2545_f491_4f6c_dd1d_u64;
    let mut next = move || {
        state ^= state << 13;
        state ^= state >> 7;
        state ^= state << 17;
        state
    };
    let mut unit = move || (next() >> 11) as f64 / (1_u64 << 53) as f64 - 0.5;
    for _ in 0..4000 {
        let a = p(unit() * 1.0e-155, unit() * 1.0e-155);
        let b = p(unit() * 1.0e-155, unit() * 1.0e-155);
        let t = unit() + 0.5;
        // c lies near the line a-b, off it by a few ulps at most.
        let on_x = a.x + (b.x - a.x) * t;
        let on_y = a.y + (b.y - a.y) * t;
        let c = p(
            on_x * (1.0 + unit() * 1.0e-15),
            on_y * (1.0 + unit() * 1.0e-15),
        );
        assert_matches_reference(a, b, c);
    }
}

/// Exactly collinear full-mantissa grid points settle in the expansion stage
/// and agree with the reference, left, right and on the line.
#[test]
fn full_mantissa_grid_triples_are_exact() {
    let step = 0.001_f64;
    let origin = 100.0_f64;
    for i in 0..40_u32 {
        let x = |k: u32| origin + f64::from(i + k) * step;
        let y = |k: u32| origin + f64::from(2 * (i + k) + 1) * step;
        let (a, b, c) = (p(x(0), y(0)), p(x(1), y(1)), p(x(2), y(2)));
        assert_matches_reference(a, b, c);
        assert_matches_reference(a, b, p(c.x, c.y + step));
        assert_matches_reference(a, b, p(c.x, c.y - step));
    }
}

/// Collinear grid triples across binades, where coordinate differences round
/// (so the tails join the expansion), agree with the reference and are Zero
/// on the line, and one grid step off it keep their side.
#[test]
fn grid_triples_across_binades_with_rounding_differences_are_exact() {
    let mut state = 0x6a09_e667_f3bc_c908_u64;
    let mut next = move || {
        state ^= state << 13;
        state ^= state >> 7;
        state ^= state << 17;
        state
    };
    let mut rounded = 0;
    for _ in 0..3000 {
        let mut vertex = || {
            let x = (next() % 200_000) as f64 * 0.001;
            p(x, 2.0 * x)
        };
        let (a, b, c) = (vertex(), vertex(), vertex());
        if (a.x - c.x) + c.x != a.x || (b.y - c.y) + c.y != b.y {
            rounded += 1;
        }
        assert_eq!(sign_of(orient2d(a, b, c)), Sign::Zero);
        assert_matches_reference(a, b, c);
        assert_matches_reference(a, b, p(c.x, c.y + 0.001));
        assert_matches_reference(a, b, p(c.x, c.y - 0.001));
    }
    assert!(
        rounded > 0,
        "no difference rounded; the tails stayed unused"
    );
}

/// Unit-size coordinates beside ones near 2^-600: every difference is in the
/// filter's range, but its rounding tail is so small that the product of two
/// tails underflows, so the expansion stage must hand these to the dyadic path.
#[test]
fn tails_too_small_for_their_products_keep_an_exact_sign() {
    let tiny = [
        0.0,
        2.0_f64.powi(-540),
        1.5 * 2.0_f64.powi(-560),
        2.0_f64.powi(-600),
        1.25 * 2.0_f64.powi(-640),
    ];
    let big = [1.0, -1.0, 2.0, 1.0 + f64::EPSILON, -3.0];
    let mut state = 0x3c6e_f372_fe94_f82b_u64;
    let mut next = move || {
        state ^= state << 13;
        state ^= state >> 7;
        state ^= state << 17;
        state
    };
    let mut coordinate = || {
        let value = if next() % 2 == 0 {
            big[(next() % 5) as usize]
        } else {
            tiny[(next() % 5) as usize]
        };
        if next() % 2 == 0 {
            value
        } else {
            -value
        }
    };
    for _ in 0..40_000 {
        let a = p(coordinate(), coordinate());
        let b = p(coordinate(), coordinate());
        let c = p(coordinate(), coordinate());
        assert_matches_reference(a, b, c);
    }
}
