//! The separation and the piece certificate of [`super`] against the sampled
//! system: `mu` must stay below the smallest singular value of the scaled
//! matrix at every sampled point of the boxes, and the speed and bend of a
//! piece must hold the sampled rates of the section.

use super::*;
use axiolid_core::Point2;
use axiolid_curve::{BSplineSurface, Carrier, KnotSpec, PairNode};

const SAMPLES: usize = 400;

fn patch(rows: Vec<Vec<Point3>>, u_degree: u16) -> BSplineSurface {
    BSplineSurface {
        u_degree,
        v_degree: 1,
        control_points: rows,
        u_knots: vec![0.0, 1.0],
        u_multiplicities: vec![u32::from(u_degree) + 1; 2],
        v_knots: vec![0.0, 1.0],
        v_multiplicities: vec![2; 2],
        weights: None,
        u_closed: false,
        v_closed: false,
        knot_spec: KnotSpec::Unspecified,
        self_intersect: None,
    }
}

/// `z = x^2` and `z = 1/4 + y/2 + k (x^2 - 1/2)`, `x = 2u - 1`, `y = v`,
/// meeting where `(1 - k) x^2 = 1/4 - k/2 + y/2`.
fn section(k: Scalar) -> PairSection3 {
    let z = [1.0, -1.0, 1.0];
    let x = [-1.0, 0.0, 1.0];
    let cylinder: Vec<Vec<Point3>> = (0..3)
        .map(|i| vec![Point3::new(x[i], 0.0, z[i]), Point3::new(x[i], 1.0, z[i])])
        .collect();
    let tilted: Vec<Vec<Point3>> = (0..3)
        .map(|i| {
            let c = 0.25 + k * (z[i] - 0.5);
            vec![Point3::new(x[i], 0.0, c), Point3::new(x[i], 1.0, c + 0.5)]
        })
        .collect();
    let node = |v: Scalar| {
        let x = ((0.25 - 0.5 * k + 0.5 * v) / (1.0 - k)).sqrt();
        let u = (1.0 + x) / 2.0;
        PairNode {
            point: Point3::new(x, v, x * x),
            first: Point2::new(u, v),
            second: Point2::new(u, v),
        }
    };
    PairSection3 {
        first: Carrier::Spline(Box::new(patch(cylinder, 2))),
        second: Carrier::Spline(Box::new(patch(tilted, 2))),
        nodes: [0.0, 0.25, 0.5, 0.75, 1.0].map(node).to_vec(),
    }
}

/// The section with its surfaces, and the nodes' parameters on them,
/// exchanged.
fn swapped(section: &PairSection3) -> PairSection3 {
    PairSection3 {
        first: section.second.clone(),
        second: section.first.clone(),
        nodes: section
            .nodes
            .iter()
            .map(|n| PairNode {
                first: n.second,
                second: n.first,
                ..*n
            })
            .collect(),
    }
}

type Matrix = [[Scalar; 4]; 4];

fn product(a: &Matrix, b: &Matrix) -> Matrix {
    a.map(|row| core::array::from_fn(|j| row.iter().zip(b).map(|(x, r)| x * r[j]).sum()))
}

fn transpose(a: &Matrix) -> Matrix {
    core::array::from_fn(|i| core::array::from_fn(|j| a[j][i]))
}

/// The smallest singular value of a 4x4 matrix: the square root of the
/// smallest eigenvalue of `M^T M`, by Jacobi rotations.
fn smallest_singular(m: Matrix) -> Scalar {
    let mut a = product(&transpose(&m), &m);
    for _ in 0..60 {
        for p in 0..3 {
            for q in p + 1..4 {
                if a[p][q].abs() < 1e-300 {
                    continue;
                }
                let theta = (a[q][q] - a[p][p]) / (2.0 * a[p][q]);
                let t = theta.signum() / (theta.abs() + (theta * theta + 1.0).sqrt());
                let (c, s) = (1.0 / (t * t + 1.0).sqrt(), t / (t * t + 1.0).sqrt());
                let mut g: Matrix = core::array::from_fn(|i| {
                    core::array::from_fn(|j| if i == j { 1.0 } else { 0.0 })
                });
                (g[p][p], g[q][q], g[p][q], g[q][p]) = (c, c, s, -s);
                a = product(&transpose(&g), &product(&a, &g));
            }
        }
    }
    (0..4)
        .map(|i| a[i][i])
        .fold(Scalar::INFINITY, Scalar::min)
        .max(0.0)
        .sqrt()
}

/// The scaled matrix at parameters `(a, b)`, scaled by the boxes' centres.
fn scaled(chord: &Chord<'_>, scales: [[Scalar; 2]; 2], a: Point2, b: Point2) -> Matrix {
    let (_, a_u, a_v) = chord.oracles[0].jet(a.x, a.y).unwrap();
    let (_, b_u, b_v) = chord.oracles[1].jet(b.x, b.y).unwrap();
    let n = chord.d / chord.length.hi;
    let columns = [
        a_u / scales[0][0],
        a_v / scales[0][1],
        -b_u / scales[1][0],
        -b_v / scales[1][1],
    ];
    let mut m = [[0.0; 4]; 4];
    for (j, c) in columns.iter().enumerate() {
        m[0][j] = c.x;
        m[1][j] = c.y;
        m[2][j] = c.z;
        m[3][j] = if j < 2 { n.dot(*c) } else { 0.0 };
    }
    m
}

/// A deterministic pseudo-random number in `[0, 1)`.
fn next(state: &mut u64) -> Scalar {
    *state = state
        .wrapping_mul(6364136223846793005)
        .wrapping_add(1442695040888963407);
    (*state >> 11) as Scalar / (1_u64 << 53) as Scalar
}

#[test]
fn mu_is_below_the_scaled_singular_values_over_the_boxes() {
    let mut rng = 11_u64;
    for k in [0.0, 0.3, -0.4] {
        let section = section(k);
        let chord = Chord::new(&section, true, 1.2, 1.4).unwrap();
        let (n0, n1) = (section.nodes[1], section.nodes[2]);
        for half in [0.004, 0.01, 0.02] {
            let around = |p: Point2, q: Point2| {
                let (lo, hi) = (p.min(q), p.max(q));
                ((lo.x - half, hi.x + half), (lo.y - half, hi.y + half))
            };
            let boxes = [around(n0.first, n1.first), around(n0.second, n1.second)];
            let jets = [
                chord.oracles[0].box_jet(boxes[0].0, boxes[0].1).unwrap(),
                chord.oracles[1].box_jet(boxes[1].0, boxes[1].1).unwrap(),
            ];
            let mu = chord.separation(&jets).expect("separated");
            let scales = jets.map(|j| j.scale);
            let mut lowest = Scalar::INFINITY;
            for _ in 0..SAMPLES {
                let mut at = |b: &((Scalar, Scalar), (Scalar, Scalar))| {
                    Point2::new(
                        b.0 .0 + (b.0 .1 - b.0 .0) * next(&mut rng),
                        b.1 .0 + (b.1 .1 - b.1 .0) * next(&mut rng),
                    )
                };
                let (a, b) = (at(&boxes[0]), at(&boxes[1]));
                lowest = lowest.min(smallest_singular(scaled(&chord, scales, a, b)));
            }
            assert!(mu <= lowest, "k {k} half {half}: mu {mu} over {lowest}");
            // Not slack beyond the boxes' spread.
            assert!(
                mu >= 0.3 * lowest || half > 0.01,
                "k {k} half {half}: {mu} {lowest}"
            );
        }
    }
}

#[test]
fn a_piece_holds_the_sampled_rates_of_the_section() {
    for k in [0.0, 0.3, -0.4] {
        // Either surface may be the curved one: the plane's own curvature
        // is zero, so the bend must come from the other's.
        for section in [section(k), swapped(&section(k))] {
            for first in [true, false] {
                let (lo, hi) = (1.2, 1.23);
                let chord = Chord::new(&section, first, lo, hi).unwrap();
                let p = chord
                    .piece(lo, hi)
                    .unwrap_or_else(|| panic!("k {k} not certified"));
                let (mut slope, mut bend): (Scalar, Scalar) = (0.0, 0.0);
                for i in 0..=SAMPLES {
                    let t = lo + (hi - lo) * i as Scalar / SAMPLES as Scalar;
                    let (da, db, _) = section.rates(t).unwrap();
                    let (dda, ddb, _) = section.second_rates(t).unwrap();
                    let (d, dd) = if first { (da, dda) } else { (db, ddb) };
                    slope = slope.max(d.length());
                    bend = bend.max(dd.length());
                }
                assert!(slope <= p.slope, "k {k} first {first}: {slope} over {p:?}");
                assert!(
                    bend <= p.bend.unwrap(),
                    "k {k} first {first}: {bend} over {p:?}"
                );
                // The nodes are exact: the ends carry rounding only.
                assert!(p.start < 1e-8 && p.end < 1e-8, "{p:?}");
            }
        }
    }
}

#[test]
fn a_wider_box_about_the_same_centre_has_a_smaller_mu() {
    // The same matrix at the centre; only the perturbation differs.
    let section = section(0.3);
    let chord = Chord::new(&section, true, 1.2, 1.4).unwrap();
    let node = section.nodes[1];
    let mu = |half: Scalar| {
        let around = |p: Point2| ((p.x - half, p.x + half), (p.y - half, p.y + half));
        let jets = [
            chord.oracles[0]
                .box_jet(around(node.first).0, around(node.first).1)
                .unwrap(),
            chord.oracles[1]
                .box_jet(around(node.second).0, around(node.second).1)
                .unwrap(),
        ];
        chord.separation(&jets).unwrap()
    };
    assert!(mu(0.001) > mu(0.01) + 1e-3, "{} {}", mu(0.001), mu(0.01));
}
