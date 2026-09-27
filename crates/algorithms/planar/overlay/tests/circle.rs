//! Minimum enclosing circle (#118). Oracles are closed forms, exact
//! containment checks in dyadics, and a brute force over every circle on
//! two or three of the points.

use axiolid_core::Point2;
use axiolid_exact::{Arith, Dyadic};
use axiolid_guarantees::Sign;
use axiolid_overlay::{minimum_enclosing_circle, CircleError, MinimumCircle};

fn p(x: f64, y: f64) -> Point2 {
    Point2::new(x, y)
}

fn exact(x: f64) -> Dyadic {
    Dyadic::from_f64(x)
}

/// Every point lies in the returned circle, decided exactly.
fn contains_all(points: &[Point2], m: &MinimumCircle) {
    let c = m.circle.centre;
    let r = exact(m.circle.radius);
    for q in points {
        let dx = exact(q.x).sub(&exact(c.x));
        let dy = exact(q.y).sub(&exact(c.y));
        let slack = r.mul(&r).sub(&dx.mul(&dx).add(&dy.mul(&dy)));
        assert_ne!(slack.sign(), Some(Sign::Negative), "{q:?} outside {m:?}");
    }
}

/// The smallest circle on two or three of the points that holds them all,
/// in plain floating point with a relative slack: an independent oracle.
fn brute_force(points: &[Point2]) -> f64 {
    let holds = |c: Point2, r: f64| {
        points
            .iter()
            .all(|q| (*q - c).length() <= r * (1.0 + 1e-12) + 1e-12)
    };
    let mut best = f64::INFINITY;
    let n = points.len();
    for i in 0..n {
        for j in i + 1..n {
            let c = (points[i] + points[j]) * 0.5;
            let r = (points[i] - c).length();
            if r < best && holds(c, r) {
                best = r;
            }
            for k in j + 1..n {
                let (a, b, cc) = (points[i], points[j], points[k]);
                let (u, v) = (b - a, cc - a);
                let d = 2.0 * (u.x * v.y - u.y * v.x);
                if d == 0.0 {
                    continue;
                }
                let (uu, vv) = (u.dot(u), v.dot(v));
                let centre = a + Point2::new(uu * v.y - vv * u.y, vv * u.x - uu * v.x) / d;
                let r = (a - centre).length();
                if r < best && holds(centre, r) {
                    best = r;
                }
            }
        }
    }
    best
}

/// A deterministic cloud.
fn cloud(n: usize, seed: u64, scale: f64) -> Vec<Point2> {
    let mut s = seed;
    let mut next = move || {
        s = s
            .wrapping_mul(6_364_136_223_846_793_005)
            .wrapping_add(1_442_695_040_888_963_407);
        ((s >> 11) as f64 / (1u64 << 53) as f64 - 0.5) * scale
    };
    (0..n).map(|_| p(next(), next())).collect()
}

#[test]
fn two_far_points_are_a_diameter() {
    let points = [
        p(0.0, 0.0),
        p(1.0, 0.5),
        p(4.0, 0.0),
        p(2.0, -1.0),
        p(3.0, 1.5),
    ];
    let m = minimum_enclosing_circle(&points).unwrap();
    assert_eq!(m.evidence.support, vec![0, 2]);
    assert_eq!(m.circle.centre, p(2.0, 0.0));
    assert!(m.circle.radius >= 2.0 && m.circle.radius - 2.0 <= m.evidence.error);
    assert!(m.evidence.error < 1e-14, "{}", m.evidence.error);
    contains_all(&points, &m);
}

#[test]
fn an_acute_triangle_gives_its_circumcircle() {
    // Centre (2, 5/6), radius 13/6.
    let points = [p(0.0, 0.0), p(4.0, 0.0), p(2.0, 3.0), p(2.0, 1.0)];
    let m = minimum_enclosing_circle(&points).unwrap();
    assert_eq!(m.evidence.support, vec![0, 1, 2]);
    let e = m.evidence.error;
    assert!(e < 1e-14, "{e}");
    assert!((m.circle.centre - p(2.0, 5.0 / 6.0)).length() <= e + 1e-16);
    assert!(m.circle.radius >= 13.0 / 6.0 - 1e-16 && m.circle.radius - 13.0 / 6.0 <= e);
    contains_all(&points, &m);
}

#[test]
fn an_obtuse_triangle_gives_its_longest_side_as_diameter() {
    let points = [p(0.0, 0.0), p(4.0, 0.0), p(2.0, 0.5)];
    let m = minimum_enclosing_circle(&points).unwrap();
    assert_eq!(m.evidence.support, vec![0, 1]);
    assert_eq!(m.circle.centre, p(2.0, 0.0));
    assert!(m.circle.radius - 2.0 <= m.evidence.error);
}

#[test]
fn points_on_a_known_circle() {
    let (centre, radius) = (p(3.0, -1.0), 5.0);
    let mut points: Vec<Point2> = (0..97)
        .map(|k| {
            let t = k as f64 * std::f64::consts::TAU / 97.0;
            centre + p(t.cos(), t.sin()) * radius
        })
        .collect();
    points.extend(cloud(200, 7, 6.0).into_iter().map(|q| centre + q));
    let m = minimum_enclosing_circle(&points).unwrap();
    assert_eq!(m.evidence.support.len(), 3);
    assert!((m.circle.radius - radius).abs() < 1e-12, "{m:?}");
    assert!((m.circle.centre - centre).length() < 1e-12);
    contains_all(&points, &m);
}

#[test]
fn a_cocircular_square_needs_two_or_three_supports() {
    let points = [p(1.0, 1.0), p(-1.0, 1.0), p(-1.0, -1.0), p(1.0, -1.0)];
    let m = minimum_enclosing_circle(&points).unwrap();
    assert!(m.evidence.support.len() >= 2);
    assert!((m.circle.centre - p(0.0, 0.0)).length() <= m.evidence.error);
    assert!(m.circle.radius >= 2f64.sqrt() - 1e-16);
    assert!(m.circle.radius - 2f64.sqrt() <= m.evidence.error + 1e-16);
    contains_all(&points, &m);
}

#[test]
fn duplicates_and_collinear_points() {
    let same = [p(2.5, -7.0); 5];
    let m = minimum_enclosing_circle(&same).unwrap();
    assert_eq!(m.circle.centre, p(2.5, -7.0));
    assert_eq!(m.circle.radius, 0.0);
    assert_eq!(m.evidence.error, 0.0);
    assert_eq!(m.evidence.support.len(), 1);

    let line = [
        p(0.0, 0.0),
        p(1.0, 0.0),
        p(5.0, 0.0),
        p(3.0, 0.0),
        p(5.0, 0.0),
    ];
    let m = minimum_enclosing_circle(&line).unwrap();
    assert_eq!(m.evidence.support.len(), 2);
    assert_eq!(m.circle.centre, p(2.5, 0.0));
    assert!(m.circle.radius >= 2.5 && m.circle.radius - 2.5 <= m.evidence.error);
    contains_all(&line, &m);
}

#[test]
fn empty_and_non_finite_input_are_refused() {
    assert_eq!(minimum_enclosing_circle(&[]), Err(CircleError::Empty));
    for bad in [f64::NAN, f64::INFINITY, f64::NEG_INFINITY] {
        assert_eq!(
            minimum_enclosing_circle(&[p(0.0, 0.0), p(bad, 1.0)]),
            Err(CircleError::NonFinite)
        );
    }
}

#[test]
fn random_clouds_match_the_brute_force_and_contain_every_point() {
    for seed in 1..=12u64 {
        let points = cloud(
            12 + seed as usize * 3,
            seed,
            10.0_f64.powi(seed as i32 % 5 - 2),
        );
        let m = minimum_enclosing_circle(&points).unwrap();
        contains_all(&points, &m);
        let oracle = brute_force(&points);
        let scale = oracle.max(1e-300);
        assert!(
            (m.circle.radius - oracle).abs() <= 1e-10 * scale,
            "seed {seed}: {} vs brute force {oracle}",
            m.circle.radius
        );
        assert!(m.evidence.error <= 1e-12 * scale, "seed {seed}: {m:?}");
        // The support points lie on the returned circle, to within the
        // error bound.
        for &i in &m.evidence.support {
            let d = (points[i] - m.circle.centre).length();
            assert!((d - m.circle.radius).abs() <= 2.0 * m.evidence.error + 1e-15 * scale);
        }
    }
}

#[test]
fn the_circle_does_not_depend_on_the_input_order() {
    let points = cloud(60, 99, 3.0);
    let forward = minimum_enclosing_circle(&points).unwrap();
    let mut reversed = points.clone();
    reversed.reverse();
    let backward = minimum_enclosing_circle(&reversed).unwrap();
    let e = forward.evidence.error + backward.evidence.error;
    assert!((forward.circle.centre - backward.circle.centre).length() <= e);
    assert!((forward.circle.radius - backward.circle.radius).abs() <= e);
}
