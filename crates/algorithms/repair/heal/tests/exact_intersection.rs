//! The triangle-triangle test is exact, with no float direction (#399).
//!
//! A rotation rounds the vertices of a flat face slightly off their common
//! plane. Two disjoint triangles of that face are then nearly coplanar, and
//! a test that projects onto a float `n_p x n_q` reads noise. Every case
//! here either pins such a configuration or checks the narrow phase against
//! an independent exact rational oracle.

use axiolid_core::Point3;
use axiolid_heal::{self_intersections, self_intersections_brute_force};
use axiolid_mesh::TriMesh;
use num_bigint::BigInt;
use std::cmp::Ordering;

type P = [f64; 3];

fn mesh(points: &[P], indices: Vec<u32>) -> TriMesh {
    TriMesh::new(
        points
            .iter()
            .map(|&[x, y, z]| Point3::new(x, y, z))
            .collect(),
        indices,
    )
}

fn rotate_z(points: &[P], degrees: f64) -> Vec<P> {
    let (s, c) = degrees.to_radians().sin_cos();
    points
        .iter()
        .map(|[x, y, z]| [x * c - y * s, x * s + y * c, *z])
        .collect()
}

fn pairs(points: &[P], indices: &[u32]) -> usize {
    let m = mesh(points, indices.to_vec());
    let found = self_intersections(&m);
    assert_eq!(found, self_intersections_brute_force(&m));
    found.len()
}

/// The issue's repro: two disjoint triangles in the plane `x = 0.025`.
const REPRO: [P; 6] = [
    [0.025, -0.96, 0.0],
    [0.025, -0.81, 0.0],
    [0.025, -0.81, 2.26],
    [0.025, 0.0, 2.7],
    [0.025, -1.77, 2.7],
    [0.025, -0.96, 2.26],
];

#[test]
fn disjoint_near_coplanar_triangles_stay_disjoint_at_every_whole_degree() {
    for degrees in 0..=90 {
        let rotated = rotate_z(&REPRO, f64::from(degrees));
        assert_eq!(pairs(&rotated, &[0, 1, 2, 3, 4, 5]), 0, "{degrees} degrees");
    }
}

/// The bottom face of a block less a through-slot, a U in the plane
/// `x = 0.025`, triangulated into six triangles that share their vertices.
const SLOT_FACE: [P; 8] = [
    [0.025, -1.77, 0.0],
    [0.025, -0.96, 0.0],
    [0.025, -0.96, 2.26],
    [0.025, -0.81, 2.26],
    [0.025, -0.81, 0.0],
    [0.025, 0.0, 0.0],
    [0.025, 0.0, 2.7],
    [0.025, -1.77, 2.7],
];
const SLOT_TRIANGLES: [u32; 18] = [0, 1, 2, 0, 2, 7, 2, 3, 6, 2, 6, 7, 3, 4, 5, 3, 5, 6];

#[test]
fn a_slotted_face_has_no_self_intersection_at_any_whole_degree() {
    for degrees in 0..=90 {
        let rotated = rotate_z(&SLOT_FACE, f64::from(degrees));
        assert_eq!(pairs(&rotated, &SLOT_TRIANGLES), 0, "{degrees} degrees");
    }
}

/// A flat blade through a flat triangle: `crossing`, edge through face,
/// vertex on face and coplanar overlap, each as `(points, name)`.
fn intersecting_cases() -> Vec<(Vec<P>, &'static str)> {
    vec![
        (
            vec![
                [0.0, 0.0, 0.0],
                [2.0, 0.0, 0.0],
                [0.0, 2.0, 0.0],
                [0.5, 0.5, -1.0],
                [0.5, 0.5, 1.0],
                [1.5, 1.5, 0.0],
            ],
            "crossing",
        ),
        (
            vec![
                [0.0, 0.0, 0.0],
                [2.0, 0.0, 0.0],
                [0.0, 2.0, 0.0],
                [0.5, 0.5, -1.0],
                [0.5, 0.5, 1.0],
                [3.0, 0.5, 3.0],
            ],
            "an edge through the face",
        ),
        (
            vec![
                [0.0, 0.0, 0.0],
                [2.0, 0.0, 0.0],
                [0.0, 2.0, 0.0],
                [0.5, 0.5, 0.0],
                [3.0, 0.5, 1.0],
                [0.5, 3.0, 1.0],
            ],
            "a vertex on the face",
        ),
        (
            vec![
                [0.0, 0.0, 0.0],
                [2.0, 0.0, 0.0],
                [0.0, 2.0, 0.0],
                [0.5, 0.5, 0.0],
                [3.0, 0.5, 0.0],
                [0.5, 3.0, 0.0],
            ],
            "coplanar overlap",
        ),
        (
            vec![
                [0.0, 0.0, 0.0],
                [2.0, 0.0, 0.0],
                [0.0, 2.0, 0.0],
                [0.0, 0.0, 0.0],
                [2.0, 0.0, 0.0],
                [0.0, 2.0, 0.0],
            ],
            "the same triangle twice",
        ),
    ]
}

#[test]
fn intersecting_pairs_are_reported_at_every_rotation() {
    for (points, name) in intersecting_cases() {
        // Tilt the frame so that the rotation about z is not axis aligned.
        let tilted: Vec<P> = points.iter().map(|[x, y, z]| [*z, *y, *x]).collect();
        for frame in [&points, &tilted] {
            for degrees in 0..=90 {
                let rotated = rotate_z(frame, f64::from(degrees));
                // Coplanar overlap is exactly planar only until it is
                // rounded; its overlap area is large, so it stays reported.
                assert_eq!(
                    pairs(&rotated, &[0, 1, 2, 3, 4, 5]),
                    1,
                    "{name} at {degrees} degrees"
                );
            }
        }
    }
}

#[test]
fn coplanar_touching_is_not_an_intersection() {
    let touching: [(&str, [P; 6]); 3] = [
        (
            "a shared edge",
            [
                [0.0, 0.0, 0.0],
                [1.0, 0.0, 0.0],
                [0.0, 1.0, 0.0],
                [1.0, 0.0, 0.0],
                [1.0, 1.0, 0.0],
                [0.0, 1.0, 0.0],
            ],
        ),
        (
            "a shared corner",
            [
                [0.0, 0.0, 0.0],
                [1.0, 0.0, 0.0],
                [0.0, 1.0, 0.0],
                [0.0, 0.0, 0.0],
                [-1.0, 0.0, 0.0],
                [0.0, -1.0, 0.0],
            ],
        ),
        (
            "a vertex on an edge",
            [
                [0.0, 0.0, 0.0],
                [2.0, 0.0, 0.0],
                [0.0, 2.0, 0.0],
                [1.0, 0.0, 0.0],
                [2.0, -1.0, 0.0],
                [0.0, -1.0, 0.0],
            ],
        ),
    ];
    for (name, points) in touching {
        assert_eq!(pairs(&points, &[0, 1, 2, 3, 4, 5]), 0, "{name}");
    }
}

// ---------------------------------------------------------------------------
// Independent exact oracle: BigInt rationals, plane-distance intervals on the
// cross-product direction in 3-D and a clipped-polygon area test in 2-D. It
// shares no code or formulation with the crate's orientation tests.

#[derive(Clone, Debug)]
struct Rat {
    n: BigInt,
    d: BigInt,
}

impl Rat {
    fn new(n: BigInt, d: BigInt) -> Self {
        if d.sign() == num_bigint::Sign::Minus {
            Self { n: -n, d: -d }
        } else {
            Self { n, d }
        }
    }
    fn int(n: &BigInt) -> Self {
        Self::new(n.clone(), BigInt::from(1))
    }
    fn add(&self, o: &Self) -> Self {
        Self::new(&self.n * &o.d + &o.n * &self.d, &self.d * &o.d)
    }
    fn sub(&self, o: &Self) -> Self {
        Self::new(&self.n * &o.d - &o.n * &self.d, &self.d * &o.d)
    }
    fn mul(&self, o: &Self) -> Self {
        Self::new(&self.n * &o.n, &self.d * &o.d)
    }
    fn cmp(&self, o: &Self) -> Ordering {
        (&self.n * &o.d).cmp(&(&o.n * &self.d))
    }
    fn sign(&self) -> i32 {
        match self.n.sign() {
            num_bigint::Sign::Minus => -1,
            num_bigint::Sign::NoSign => 0,
            num_bigint::Sign::Plus => 1,
        }
    }
}

type V = [BigInt; 3];

fn sub3(a: &V, b: &V) -> V {
    [&a[0] - &b[0], &a[1] - &b[1], &a[2] - &b[2]]
}
fn dot3(a: &V, b: &V) -> BigInt {
    &a[0] * &b[0] + &a[1] * &b[1] + &a[2] * &b[2]
}
fn cross3(a: &V, b: &V) -> V {
    [
        &a[1] * &b[2] - &a[2] * &b[1],
        &a[2] * &b[0] - &a[0] * &b[2],
        &a[0] * &b[1] - &a[1] * &b[0],
    ]
}
fn sign_of(n: &BigInt) -> i32 {
    Rat::int(n).sign()
}

/// The exact integer image of finite doubles: every coordinate times one
/// common power of two, which is exact because a double is `m * 2^e`.
fn exact(points: &[P]) -> Vec<V> {
    let decode = |x: f64| -> (BigInt, i32) {
        if x == 0.0 {
            return (BigInt::from(0), 0);
        }
        let bits = x.to_bits();
        let exponent = ((bits >> 52) & 0x7ff) as i32;
        let fraction = (bits & ((1u64 << 52) - 1)) as i64;
        let (mantissa, e) = if exponent == 0 {
            (fraction, -1074)
        } else {
            (fraction | (1i64 << 52), exponent - 1075)
        };
        let m = BigInt::from(mantissa);
        (if x < 0.0 { -m } else { m }, e)
    };
    let parts: Vec<[(BigInt, i32); 3]> = points
        .iter()
        .map(|p| [decode(p[0]), decode(p[1]), decode(p[2])])
        .collect();
    let low = parts.iter().flatten().map(|(_, e)| *e).min().unwrap();
    parts
        .into_iter()
        .map(|p| p.map(|(m, e)| m << (e - low) as usize))
        .collect()
}

fn oracle(points: &[P]) -> bool {
    let v = exact(points);
    let (p, q) = (&v[0..3], &v[3..6]);
    let np = cross3(&sub3(&p[1], &p[0]), &sub3(&p[2], &p[0]));
    let nq = cross3(&sub3(&q[1], &q[0]), &sub3(&q[2], &q[0]));
    let zero = |n: &V| n.iter().all(|c| sign_of(c) == 0);
    if zero(&np) || zero(&nq) {
        return false;
    }
    let dq: Vec<BigInt> = q.iter().map(|x| dot3(&np, &sub3(x, &p[0]))).collect();
    let dp: Vec<BigInt> = p.iter().map(|x| dot3(&nq, &sub3(x, &q[0]))).collect();
    if dq.iter().all(|d| sign_of(d) == 0) {
        return coplanar_oracle(p, q, &np);
    }
    let direction = cross3(&np, &nq);
    if zero(&direction) {
        return false; // parallel, distinct planes
    }
    let (Some(a), Some(b)) = (interval(p, &dp, &direction), interval(q, &dq, &direction)) else {
        return false;
    };
    a.0.cmp(&b.1) != Ordering::Greater && b.0.cmp(&a.1) != Ordering::Greater
}

/// The triangle's extent along `direction` on the other plane's line, where
/// `d` are its vertices' signed distances to that plane.
fn interval(t: &[V], d: &[BigInt], direction: &V) -> Option<(Rat, Rat)> {
    let mut hits: Vec<Rat> = Vec::new();
    for i in 0..3 {
        let j = (i + 1) % 3;
        if sign_of(&d[i]) == 0 {
            hits.push(Rat::int(&dot3(&t[i], direction)));
        }
        if sign_of(&d[i]) * sign_of(&d[j]) < 0 {
            let denominator = &d[i] - &d[j];
            let numerator = dot3(&t[i], direction) * &denominator
                + dot3(&sub3(&t[j], &t[i]), direction) * &d[i];
            hits.push(Rat::new(numerator, denominator));
        }
    }
    let low = hits
        .iter()
        .cloned()
        .reduce(|a, b| if a.cmp(&b).is_le() { a } else { b })?;
    let high = hits
        .iter()
        .cloned()
        .reduce(|a, b| if a.cmp(&b).is_ge() { a } else { b })?;
    Some((low, high))
}

type R2 = [Rat; 2];

fn coplanar_oracle(p: &[V], q: &[V], normal: &V) -> bool {
    // Drop the axis of the largest normal component (exact integers).
    let axis = (0..3)
        .max_by_key(|&i| normal[i].magnitude().clone())
        .unwrap();
    let (u, w) = ((axis + 1) % 3, (axis + 2) % 3);
    let flat = |t: &[V]| -> Vec<R2> {
        t.iter()
            .map(|x| [Rat::int(&x[u]), Rat::int(&x[w])])
            .collect()
    };
    let (p, q) = (flat(p), flat(q));
    let turn = |t: &[R2]| cross2(&t[0], &t[1], &t[2]).sign();
    let (tp, tq) = (turn(&p), turn(&q));
    // Closed membership and strict interior, by signed areas.
    let area =
        |t: &[R2], x: &R2, i: usize, turn: i32| cross2(&t[i], &t[(i + 1) % 3], x).sign() * turn;
    let inside = |t: &[R2], turn: i32, x: &R2, strict: bool| {
        (0..3).all(|i| {
            let s = area(t, x, i, turn);
            if strict {
                s > 0
            } else {
                s >= 0
            }
        })
    };
    let mut candidates: Vec<R2> = Vec::new();
    for x in &p {
        if inside(&q, tq, x, false) {
            candidates.push(x.clone());
        }
    }
    for x in &q {
        if inside(&p, tp, x, false) {
            candidates.push(x.clone());
        }
    }
    for i in 0..3 {
        for j in 0..3 {
            let (a, b) = (&p[i], &p[(i + 1) % 3]);
            let (c, d) = (&q[j], &q[(j + 1) % 3]);
            let r = [b[0].sub(&a[0]), b[1].sub(&a[1])];
            let s = [d[0].sub(&c[0]), d[1].sub(&c[1])];
            let den = r[0].mul(&s[1]).sub(&r[1].mul(&s[0]));
            if den.sign() == 0 {
                continue;
            }
            let ca = [c[0].sub(&a[0]), c[1].sub(&a[1])];
            let t = div(&ca[0].mul(&s[1]).sub(&ca[1].mul(&s[0])), &den);
            let u = div(&ca[0].mul(&r[1]).sub(&ca[1].mul(&r[0])), &den);
            let one = Rat::int(&BigInt::from(1));
            let zero = Rat::int(&BigInt::from(0));
            if t.cmp(&zero).is_ge()
                && t.cmp(&one).is_le()
                && u.cmp(&zero).is_ge()
                && u.cmp(&one).is_le()
            {
                candidates.push([a[0].add(&r[0].mul(&t)), a[1].add(&r[1].mul(&t))]);
            }
        }
    }
    if candidates.is_empty() {
        return false;
    }
    // The mean of the vertices of a convex region is interior to it exactly
    // when the region has area.
    let count = Rat::int(&BigInt::from(candidates.len()));
    let mut mean = [Rat::int(&BigInt::from(0)), Rat::int(&BigInt::from(0))];
    for c in &candidates {
        mean = [mean[0].add(&c[0]), mean[1].add(&c[1])];
    }
    let mean = [div(&mean[0], &count), div(&mean[1], &count)];
    inside(&p, tp, &mean, true) && inside(&q, tq, &mean, true)
}

fn cross2(a: &R2, b: &R2, c: &R2) -> Rat {
    b[0].sub(&a[0])
        .mul(&c[1].sub(&a[1]))
        .sub(&b[1].sub(&a[1]).mul(&c[0].sub(&a[0])))
}

fn div(a: &Rat, b: &Rat) -> Rat {
    Rat::new(&a.n * &b.d, &a.d * &b.n)
}

// ---------------------------------------------------------------------------

struct Rng(u64);

impl Rng {
    fn next(&mut self) -> u64 {
        self.0 ^= self.0 << 13;
        self.0 ^= self.0 >> 7;
        self.0 ^= self.0 << 17;
        self.0
    }
    fn range(&mut self, low: i64, high: i64) -> i64 {
        low + (self.next() % (high - low + 1) as u64) as i64
    }
}

fn check_against_oracle(points: &[P], label: &str) -> bool {
    let expected = oracle(points);
    let got = pairs(points, &[0, 1, 2, 3, 4, 5]) == 1;
    assert_eq!(got, expected, "{label}: {points:?}");
    expected
}

#[test]
fn the_narrow_phase_agrees_with_an_exact_rational_oracle() {
    let mut rng = Rng(0x9e37_79b9_7f4a_7c15);
    let (mut hits, mut misses) = (0u32, 0u32);
    let mut tally = |hit: bool| {
        if hit {
            hits += 1;
        } else {
            misses += 1;
        }
    };
    for round in 0..3000 {
        let unit = 1.0 / 16.0;
        let point = |rng: &mut Rng, span: i64, depth: i64| -> P {
            [
                rng.range(-span, span) as f64 * unit,
                rng.range(-span, span) as f64 * unit,
                rng.range(-depth, depth) as f64 * unit,
            ]
        };
        // Small grids give many exact contacts and coplanar pairs.
        let tiny: Vec<P> = (0..6).map(|_| point(&mut rng, 3, 3)).collect();
        tally(check_against_oracle(&tiny, "grid"));
        // Flat in z, or one unit off: exactly or nearly coplanar pairs.
        let flat: Vec<P> = (0..6)
            .map(|i| point(&mut rng, 40, if i < 3 { 0 } else { 1 }))
            .collect();
        tally(check_against_oracle(&flat, "flat"));
        // The same, rotated: rounding takes them off their common plane.
        let degrees = rng.range(0, 359) as f64 + round as f64 * 1e-3;
        let rotated = rotate_z(&flat, degrees);
        tally(check_against_oracle(&rotated, "rotated flat"));
        // A tilted second plane crossing the first near its rim.
        let tilted: Vec<P> = (0..6)
            .map(|i| {
                let [x, y, z] = point(&mut rng, 30, 0);
                if i < 3 {
                    [x, y, z]
                } else {
                    [x, y * 0.25, y * 0.125 + z]
                }
            })
            .collect();
        tally(check_against_oracle(&tilted, "tilted"));
    }
    assert!(hits > 500 && misses > 500, "{hits} hits, {misses} misses");
}

/// Coplanar triangles at coordinates whose `orient2d` products overflow
/// (1e154 and above) are decided exactly, not dropped or panicking (#404).
#[test]
fn coplanar_triangles_at_overflowing_scale_are_decided_exactly() {
    let overlapping: [P; 6] = [
        [0.0, 0.0, 0.0],
        [2.0, 0.0, 0.0],
        [0.0, 2.0, 0.0],
        [0.5, 0.5, 0.0],
        [3.0, 0.5, 0.0],
        [0.5, 3.0, 0.0],
    ];
    let disjoint: [P; 6] = [
        [0.0, 0.0, 0.0],
        [1.0, 0.0, 0.0],
        [0.0, 1.0, 0.0],
        [2.0, 2.0, 0.0],
        [3.0, 2.0, 0.0],
        [2.0, 3.0, 0.0],
    ];
    for scale in [1.0, 1.0e160, 1.0e200, 1.0e300] {
        let scaled = |points: [P; 6]| points.map(|[x, y, z]| [x * scale, y * scale, z]);
        let indices = [0, 1, 2, 3, 4, 5];
        assert_eq!(pairs(&scaled(overlapping), &indices), 1, "scale {scale}");
        assert_eq!(pairs(&scaled(disjoint), &indices), 0, "scale {scale}");
    }
}
