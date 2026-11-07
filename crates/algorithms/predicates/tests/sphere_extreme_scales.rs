//! `incircle` and `insphere` at extreme coordinate scales (#430).
//!
//! The exact path used to overflow or underflow beyond about 1e60 and report a
//! certain Zero for points clearly inside or outside the ball. Every answer
//! here is compared with an exact `Dyadic` determinant, which has no exponent
//! range.

use axiolid_core::{Point2, Point3};
use axiolid_exact::{Arith, Dyadic};
use axiolid_guarantees::{Certified, Precision, Sign};
use axiolid_predicates::{incircle, incircle_filter, insphere, insphere_filter};

fn d(x: f64) -> Dyadic {
    Dyadic::from_f64(x)
}

/// Determinant by Laplace expansion along the first row.
fn det(m: &[Vec<Dyadic>]) -> Dyadic {
    if m.len() == 1 {
        return m[0][0].clone();
    }
    let mut total = d(0.0);
    for column in 0..m.len() {
        let minor: Vec<Vec<Dyadic>> = m[1..]
            .iter()
            .map(|row| {
                row.iter()
                    .enumerate()
                    .filter(|(k, _)| *k != column)
                    .map(|(_, v)| v.clone())
                    .collect()
            })
            .collect();
        let term = m[0][column].mul(&det(&minor));
        total = if column % 2 == 0 {
            total.add(&term)
        } else {
            total.sub(&term)
        };
    }
    total
}

/// Rows `(p - origin, |p - origin|^2)`, the matrix both predicates take the
/// determinant of.
fn lifted(points: &[&[f64]], origin: &[f64]) -> Vec<Vec<Dyadic>> {
    points
        .iter()
        .map(|p| {
            let diffs: Vec<Dyadic> = p
                .iter()
                .zip(origin)
                .map(|(x, o)| d(*x).sub(&d(*o)))
                .collect();
            let lift = diffs.iter().fold(d(0.0), |sum, v| sum.add(&v.mul(v)));
            let mut row = diffs;
            row.push(lift);
            row
        })
        .collect()
}

fn incircle_reference(a: Point2, b: Point2, c: Point2, p: Point2) -> Sign {
    det(&lifted(
        &[&[a.x, a.y], &[b.x, b.y], &[c.x, c.y]],
        &[p.x, p.y],
    ))
    .sign()
    .unwrap()
}

fn insphere_reference(a: Point3, b: Point3, c: Point3, e: Point3, p: Point3) -> Sign {
    det(&lifted(
        &[
            &[a.x, a.y, a.z],
            &[b.x, b.y, b.z],
            &[c.x, c.y, c.z],
            &[e.x, e.y, e.z],
        ],
        &[p.x, p.y, p.z],
    ))
    .sign()
    .unwrap()
}

fn certain(c: Certified) -> Sign {
    match c {
        Certified::Certain { sign, .. } => sign,
        other => panic!("{other:?}"),
    }
}

fn p2(x: f64, y: f64) -> Point2 {
    Point2::new(x, y)
}

fn p3(x: f64, y: f64, z: f64) -> Point3 {
    Point3::new(x, y, z)
}

/// The issue's inputs: a circle through (s, 0), (0, s), (-s, 0).
#[test]
fn a_point_clearly_inside_or_outside_is_not_reported_on_the_circle() {
    for s in [
        1.0, 1e-80, 1e80, 1e100, 1e160, 1e-100, 1e-160, 1e300, 1e-300,
    ] {
        let (a, b, c) = (p2(s, 0.0), p2(0.0, s), p2(-s, 0.0));
        assert_eq!(
            certain(incircle(a, b, c, p2(0.0, 0.0))),
            Sign::Positive,
            "s = {s:e}"
        );
        assert_eq!(
            certain(incircle(a, b, c, p2(2.0 * s, 2.0 * s))),
            Sign::Negative,
            "s = {s:e}"
        );
    }
}

/// The 3D analogue: a sphere through the six points (+-s, 0, 0), (0, +-s, 0),
/// (0, 0, +-s), taken four at a time.
#[test]
fn a_point_clearly_inside_or_outside_is_not_reported_on_the_sphere() {
    for s in [
        1.0, 1e-60, 1e60, 1e80, 1e100, 1e160, 1e-100, 1e-160, 1e300, 1e-300,
    ] {
        let (a, b, c, e) = (
            p3(s, 0.0, 0.0),
            p3(0.0, s, 0.0),
            p3(-s, 0.0, 0.0),
            p3(0.0, 0.0, s),
        );
        let inside = insphere(a, b, c, e, p3(0.0, 0.0, 0.0));
        let outside = insphere(a, b, c, e, p3(2.0 * s, 2.0 * s, 2.0 * s));
        let (inside, outside) = (certain(inside), certain(outside));
        assert_eq!(
            inside,
            insphere_reference(a, b, c, e, p3(0.0, 0.0, 0.0)),
            "s = {s:e}"
        );
        assert_ne!(inside, Sign::Zero, "s = {s:e}");
        assert_eq!(
            outside,
            insphere_reference(a, b, c, e, p3(2.0 * s, 2.0 * s, 2.0 * s))
        );
        assert_eq!(outside, inside.flip(), "s = {s:e}");
    }
}

/// Exactly cocircular and cospherical points keep a certain Zero at any scale.
#[test]
fn exactly_cocircular_points_stay_zero_at_every_scale() {
    for k in [-1000, -300, -100, -10, 0, 10, 100, 266, 1000] {
        let s = 2f64.powi(k);
        let (a, b, c) = (
            p2(5.0 * s, 0.0),
            p2(3.0 * s, 4.0 * s),
            p2(-3.0 * s, 4.0 * s),
        );
        assert_eq!(
            certain(incircle(a, b, c, p2(4.0 * s, 3.0 * s))),
            Sign::Zero,
            "k = {k}"
        );
        assert_eq!(
            certain(incircle(a, b, c, p2(-4.0 * s, -3.0 * s))),
            Sign::Zero,
            "k = {k}"
        );
        let (a, b, c, e) = (
            p3(3.0 * s, 0.0, 0.0),
            p3(0.0, 3.0 * s, 0.0),
            p3(0.0, 0.0, 3.0 * s),
            p3(-3.0 * s, 0.0, 0.0),
        );
        assert_eq!(
            certain(insphere(a, b, c, e, p3(s, 2.0 * s, 2.0 * s))),
            Sign::Zero,
            "k = {k}"
        );
    }
}

#[test]
fn a_non_finite_coordinate_is_uncertain() {
    let uncertain = |c: Certified| {
        assert!(
            matches!(
                c,
                Certified::Uncertain {
                    attempted: Precision::Exact
                }
            ),
            "{c:?}"
        );
    };
    for bad in [f64::NAN, f64::INFINITY, f64::NEG_INFINITY] {
        for slot in 0..8 {
            let mut v = [1.0, 0.0, 0.0, 1.0, -1.0, 0.0, 0.25, 0.25];
            v[slot] = bad;
            let q = |i: usize| p2(v[2 * i], v[2 * i + 1]);
            assert!(!incircle_filter(q(0), q(1), q(2), q(3)).is_certain());
            uncertain(incircle(q(0), q(1), q(2), q(3)));
        }
        for slot in 0..15 {
            let mut v = [
                1.0, 0.0, 0.0, 0.0, 1.0, 0.0, -1.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.1, 0.2, 0.3,
            ];
            v[slot] = bad;
            let q = |i: usize| p3(v[3 * i], v[3 * i + 1], v[3 * i + 2]);
            assert!(!insphere_filter(q(0), q(1), q(2), q(3), q(4)).is_certain());
            uncertain(insphere(q(0), q(1), q(2), q(3), q(4)));
        }
    }
}

/// Differences beyond the filter's range are declined, not evaluated.
#[test]
fn the_filter_declines_outside_its_range() {
    let (a, b, c) = (p2(1e100, 0.0), p2(0.0, 1e100), p2(-1e100, 0.0));
    assert!(!incircle_filter(a, b, c, p2(0.0, 0.0)).is_certain());
    let (a, b, c) = (p2(1e-100, 0.0), p2(0.0, 1e-100), p2(-1e-100, 0.0));
    assert!(!incircle_filter(a, b, c, p2(0.0, 0.0)).is_certain());
    let (a, b, c, e) = (
        p3(1e70, 0.0, 0.0),
        p3(0.0, 1e70, 0.0),
        p3(-1e70, 0.0, 0.0),
        p3(0.0, 0.0, 1e70),
    );
    assert!(!insphere_filter(a, b, c, e, p3(0.0, 0.0, 0.0)).is_certain());
    let (a, b, c, e) = (
        p3(1e-70, 0.0, 0.0),
        p3(0.0, 1e-70, 0.0),
        p3(-1e-70, 0.0, 0.0),
        p3(0.0, 0.0, 1e-70),
    );
    assert!(!insphere_filter(a, b, c, e, p3(0.0, 0.0, 0.0)).is_certain());
}

/// A point a tiny offset `t` off the circle: the differences round, with a
/// tail far below the expansion stage's range, and the sign hinges on it.
#[test]
fn a_rounding_tail_far_below_its_difference_is_kept() {
    for k in [-200, -300, -700, -1000] {
        let t = 2f64.powi(k);
        let (a, b, c) = (p2(1.0, 0.0), p2(-1.0, 0.0), p2(0.0, 1.0));
        for p in [
            p2(t, 1.0),
            p2(-t, 1.0),
            p2(t, -1.0),
            p2(1.0, t),
            p2(1.0 - 2f64.powi(-53), t),
        ] {
            assert_eq!(
                certain(incircle(a, b, c, p)),
                incircle_reference(a, b, c, p),
                "k = {k}: {p:?}"
            );
        }
        let (a, b, c, e) = (
            p3(1.0, 0.0, 0.0),
            p3(-1.0, 0.0, 0.0),
            p3(0.0, 1.0, 0.0),
            p3(0.0, 0.0, 1.0),
        );
        for p in [
            p3(t, 0.0, -1.0),
            p3(0.0, t, -1.0),
            p3(0.0, 0.0, -1.0 + t),
            p3(t, t, 1.0),
        ] {
            assert_eq!(
                certain(insphere(a, b, c, e, p)),
                insphere_reference(a, b, c, e, p),
                "k = {k}: {p:?}"
            );
        }
    }
}

struct Rng(u64);

impl Rng {
    fn next(&mut self) -> u64 {
        self.0 ^= self.0 << 13;
        self.0 ^= self.0 >> 7;
        self.0 ^= self.0 << 17;
        self.0
    }

    fn below(&mut self, n: usize) -> usize {
        (self.next() % n as u64) as usize
    }

    fn unit(&mut self) -> f64 {
        (self.next() >> 11) as f64 / (1u64 << 53) as f64 * 2.0 - 1.0
    }
}

/// One ulp towards zero or away from it: a point just off a lattice circle.
fn nudge(v: f64, rng: &mut Rng) -> f64 {
    if v == 0.0 {
        return v;
    }
    let bits = v.to_bits();
    f64::from_bits(if rng.next().is_multiple_of(2) {
        bits + 1
    } else {
        bits - 1
    })
}

/// Zero, or a signed offset of 2^-20 to 2^-400.
fn tiny(rng: &mut Rng) -> f64 {
    if rng.below(3) == 0 {
        return 0.0;
    }
    let magnitude = 2f64.powi(-(20 + rng.below(381) as i32));
    magnitude * rng.unit().signum()
}

fn lattice2(radius: i64) -> Vec<(f64, f64)> {
    (-radius..=radius)
        .flat_map(|x| (-radius..=radius).map(move |y| (x, y)))
        .filter(|(x, y)| x * x + y * y == radius * radius)
        .map(|(x, y)| (x as f64, y as f64))
        .collect()
}

fn lattice3(radius: i64) -> Vec<(f64, f64, f64)> {
    let mut out = Vec::new();
    for x in -radius..=radius {
        for y in -radius..=radius {
            for z in -radius..=radius {
                if x * x + y * y + z * z == radius * radius {
                    out.push((x as f64, y as f64, z as f64));
                }
            }
        }
    }
    out
}

/// A scale 2^k: half the time anywhere in the range, half the time within
/// the windows where the degree-4 (`degree` 4) or degree-5 products cross
/// the overflow or underflow threshold, where partial underflow is possible.
fn scale_exponent(rng: &mut Rng, degree: i32) -> i32 {
    if rng.below(2) == 0 {
        return rng.below(2071) as i32 - 1060;
    }
    let edge = 1024 / degree;
    let offset = rng.below(81) as i32 - 40;
    if rng.below(2) == 0 {
        edge + offset
    } else {
        -edge + offset
    }
}

/// Seeded fuzz of `incircle` over scales 2^-1060 to 2^1010: lattice circle
/// points (exactly cocircular), with and without a one-ulp nudge, and random
/// points, every answer compared with the exact determinant.
#[test]
fn incircle_matches_the_dyadic_reference_across_scales() {
    let circle = lattice2(65);
    let mut rng = Rng(0x9E37_79B9_7F4A_7C15);
    let (mut cases, mut zeros, mut nonzero_out_of_filter) = (0usize, 0usize, 0usize);
    for _ in 0..30_000 {
        let k = scale_exponent(&mut rng, 4);
        let s = 2f64.powi(k);
        let pick = |rng: &mut Rng| circle[rng.below(circle.len())];
        let mut points = [(0.0, 0.0); 4];
        let mode = rng.below(4);
        for q in &mut points {
            *q = if mode == 2 {
                (rng.unit() * 65.0, rng.unit() * 65.0)
            } else {
                pick(&mut rng)
            };
        }
        let mut pts = points.map(|(x, y)| p2(x * s, y * s));
        if mode == 3 {
            // Microscopic offsets: the differences round, with tails of every size.
            for q in &mut pts {
                q.x += s * tiny(&mut rng);
                q.y += s * tiny(&mut rng);
            }
        }
        if mode == 1 {
            let which = rng.below(4);
            let (x, y) = (pts[which].x, pts[which].y);
            pts[which] = if rng.next().is_multiple_of(2) {
                p2(nudge(x, &mut rng), y)
            } else {
                p2(x, nudge(y, &mut rng))
            };
        }
        let expected = incircle_reference(pts[0], pts[1], pts[2], pts[3]);
        let got = certain(incircle(pts[0], pts[1], pts[2], pts[3]));
        assert_eq!(got, expected, "k = {k}: {pts:?}");
        cases += 1;
        zeros += usize::from(expected == Sign::Zero);
        nonzero_out_of_filter += usize::from(
            expected != Sign::Zero && !incircle_filter(pts[0], pts[1], pts[2], pts[3]).is_certain(),
        );
    }
    println!(
        "incircle fuzz: {cases} cases, {zeros} exactly zero, {nonzero_out_of_filter} nonzero beyond the filter"
    );
    assert!(zeros > 2_000 && nonzero_out_of_filter > 2_000);
}

/// The same for `insphere`, over scales 2^-1060 to 2^1010.
#[test]
fn insphere_matches_the_dyadic_reference_across_scales() {
    let sphere = lattice3(21);
    let mut rng = Rng(0xD1B5_4A32_D192_ED03);
    let (mut cases, mut zeros, mut nonzero_out_of_filter) = (0usize, 0usize, 0usize);
    for _ in 0..4_000 {
        let k = scale_exponent(&mut rng, 5);
        let s = 2f64.powi(k);
        let mode = rng.below(4);
        let mut points = [(0.0, 0.0, 0.0); 5];
        for q in &mut points {
            *q = if mode == 2 {
                (rng.unit() * 21.0, rng.unit() * 21.0, rng.unit() * 21.0)
            } else {
                sphere[rng.below(sphere.len())]
            };
        }
        let mut pts = points.map(|(x, y, z)| p3(x * s, y * s, z * s));
        if mode == 3 {
            for q in &mut pts {
                q.x += s * tiny(&mut rng);
                q.y += s * tiny(&mut rng);
                q.z += s * tiny(&mut rng);
            }
        }
        if mode == 1 {
            let which = rng.below(5);
            let (x, y, z) = (pts[which].x, pts[which].y, pts[which].z);
            pts[which] = match rng.below(3) {
                0 => p3(nudge(x, &mut rng), y, z),
                1 => p3(x, nudge(y, &mut rng), z),
                _ => p3(x, y, nudge(z, &mut rng)),
            };
        }
        let expected = insphere_reference(pts[0], pts[1], pts[2], pts[3], pts[4]);
        let got = certain(insphere(pts[0], pts[1], pts[2], pts[3], pts[4]));
        assert_eq!(got, expected, "k = {k}: {pts:?}");
        cases += 1;
        zeros += usize::from(expected == Sign::Zero);
        nonzero_out_of_filter += usize::from(
            expected != Sign::Zero
                && !insphere_filter(pts[0], pts[1], pts[2], pts[3], pts[4]).is_certain(),
        );
    }
    println!(
        "insphere fuzz: {cases} cases, {zeros} exactly zero, {nonzero_out_of_filter} nonzero beyond the filter"
    );
    assert!(zeros > 300 && nonzero_out_of_filter > 300);
}

/// Differences of very different magnitudes: a degree-2 sub-product underflows
/// and is then multiplied by a large lift, which a gate on the permanent alone
/// did not exclude (#430).
#[test]
fn mixed_scale_inputs_are_not_certified_wrong() {
    let a = p2(1e100, 0.0);
    let b = p2(0.0, -5e-324);
    let c = p2(0.3, 0.0);
    let origin = p2(0.0, 0.0);
    assert_eq!(
        certain(incircle(a, b, c, origin)),
        incircle_reference(a, b, c, origin)
    );
    assert_eq!(incircle_reference(a, b, c, origin), Sign::Positive);

    let (a, b, c) = (p2(0.0, 0.0), p2(0.0, 3.0000000000000002e100), p2(0.0, 0.3));
    let p = p2(-5e-324, 0.0);
    assert_eq!(
        certain(incircle(a, b, c, p)),
        incircle_reference(a, b, c, p)
    );
    assert_eq!(incircle_reference(a, b, c, p), Sign::Negative);

    let a = p3(-0.0625, 1e100, -7.888609052210118e-31);
    let b = p3(-1.5e-323, -0.0, 1.5e-323);
    let c = p3(-5e-324, 5e-324, -0.6);
    let e = p3(0.0, 1.88079096131566e-37, 7.427640235712282e27);
    let p = p3(-0.0, 1e-323, -0.0);
    assert_eq!(
        certain(insphere(a, b, c, e, p)),
        insphere_reference(a, b, c, e, p)
    );
    assert_eq!(insphere_reference(a, b, c, e, p), Sign::Positive);
}

/// A coordinate of random sign with a random binary exponent from the
/// subnormals to about 2^340, or zero.
fn mixed(rng: &mut Rng) -> f64 {
    if rng.below(8) == 0 {
        return 0.0;
    }
    let k = rng.below(1500) as i32 - 1070;
    rng.unit() * 2f64.powi(k.min(340))
}

/// Seeded fuzz with an independent scale per coordinate, so that differences
/// of every magnitude meet in one determinant.
#[test]
fn mixed_scale_fuzz_matches_the_dyadic_reference() {
    let mut rng = Rng(0xA5A5_1234_DEAD_BEEF);
    for _ in 0..40_000 {
        let pts = [(); 4].map(|()| p2(mixed(&mut rng), mixed(&mut rng)));
        assert_eq!(
            certain(incircle(pts[0], pts[1], pts[2], pts[3])),
            incircle_reference(pts[0], pts[1], pts[2], pts[3]),
            "{pts:?}"
        );
    }
    for _ in 0..4_000 {
        let pts = [(); 5].map(|()| p3(mixed(&mut rng), mixed(&mut rng), mixed(&mut rng)));
        assert_eq!(
            certain(insphere(pts[0], pts[1], pts[2], pts[3], pts[4])),
            insphere_reference(pts[0], pts[1], pts[2], pts[3], pts[4]),
            "{pts:?}"
        );
    }
}

/// The mixed-scale inputs under every permutation of their points and
/// coordinate axes, so that each of the filters' differences is the one that
/// carries the tiny value.
#[test]
fn mixed_scale_inputs_hold_under_every_permutation() {
    fn permutations(n: usize) -> Vec<Vec<usize>> {
        if n == 1 {
            return vec![vec![0]];
        }
        let mut out = Vec::new();
        for rest in permutations(n - 1) {
            for at in 0..n {
                let mut p = rest.clone();
                p.insert(at, n - 1);
                out.push(p);
            }
        }
        out
    }
    let flat2: [[[f64; 2]; 4]; 2] = [
        [[1e100, 0.0], [0.0, -5e-324], [0.3, 0.0], [0.0, 0.0]],
        [
            [0.0, 0.0],
            [0.0, 3.0000000000000002e100],
            [0.0, 0.3],
            [-5e-324, 0.0],
        ],
    ];
    for case in flat2 {
        for order in permutations(4) {
            // Each point swaps its own axes or not, so the large and the tiny
            // value meet in every combination of difference slots.
            for swaps in 0..16_u32 {
                let q = |i: usize| {
                    let v = case[order[i]];
                    if swaps & (1 << i) != 0 {
                        p2(v[1], v[0])
                    } else {
                        p2(v[0], v[1])
                    }
                };
                assert_eq!(
                    certain(incircle(q(0), q(1), q(2), q(3))),
                    incircle_reference(q(0), q(1), q(2), q(3)),
                    "{order:?} {swaps:b} {case:?}"
                );
            }
        }
    }
    let case: [[f64; 3]; 5] = [
        [-0.0625, 1e100, -7.888609052210118e-31],
        [-1.5e-323, -0.0, 1.5e-323],
        [-5e-324, 5e-324, -0.6],
        [0.0, 1.88079096131566e-37, 7.427640235712282e27],
        [-0.0, 1e-323, -0.0],
    ];
    let axes = permutations(3);
    let mut rng = Rng(0x0DDB_1A5E_5BAD_5EED);
    for order in permutations(5) {
        for _ in 0..60 {
            let pick: Vec<&Vec<usize>> = (0..5).map(|_| &axes[rng.below(axes.len())]).collect();
            let q = |i: usize| {
                let (v, a) = (case[order[i]], pick[i]);
                p3(v[a[0]], v[a[1]], v[a[2]])
            };
            assert_eq!(
                certain(insphere(q(0), q(1), q(2), q(3), q(4))),
                insphere_reference(q(0), q(1), q(2), q(3), q(4)),
                "{order:?} {pick:?}"
            );
        }
    }
}

/// Seeded fuzz over a small palette of awkward values (zeros, subnormals,
/// large and mid-range ones) so that sub-products underflow beside large
/// lifts in every slot of the determinant.
#[test]
fn palette_fuzz_matches_the_dyadic_reference() {
    const PALETTE: [f64; 14] = [
        0.0,
        -0.0,
        5e-324,
        -5e-324,
        1.5e-323,
        1e-323,
        0.3,
        -0.6,
        1e100,
        -3.0000000000000002e100,
        1.88079096131566e-37,
        -7.888609052210118e-31,
        7.427640235712282e27,
        1.0,
    ];
    let mut rng = Rng(0xFEED_FACE_CAFE_F00D);
    let pick = |rng: &mut Rng| PALETTE[rng.below(PALETTE.len())];
    for _ in 0..150_000 {
        let pts = [(); 4].map(|()| p2(pick(&mut rng), pick(&mut rng)));
        assert_eq!(
            certain(incircle(pts[0], pts[1], pts[2], pts[3])),
            incircle_reference(pts[0], pts[1], pts[2], pts[3]),
            "{pts:?}"
        );
    }
    for _ in 0..120_000 {
        let pts = [(); 5].map(|()| p3(pick(&mut rng), pick(&mut rng), pick(&mut rng)));
        assert_eq!(
            certain(insphere(pts[0], pts[1], pts[2], pts[3], pts[4])),
            insphere_reference(pts[0], pts[1], pts[2], pts[3], pts[4]),
            "{pts:?}"
        );
    }
}
