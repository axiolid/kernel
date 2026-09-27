//! Minimum enclosing sphere and oriented bounding box (#118).
//!
//! Oracles: closed forms, exact containment checks in dyadics, and a brute
//! force over every sphere on two, three or four of the points.

use axiolid_construct::bounding::{
    minimum_enclosing_sphere, oriented_bounding_box, BoundingError, MinimumSphere, OrientedBox,
};
use axiolid_core::{Point2, Point3, Vec3};
use axiolid_exact::{Arith, Dyadic};
use axiolid_guarantees::Sign;

fn p(x: f64, y: f64, z: f64) -> Point3 {
    Point3::new(x, y, z)
}

fn exact(x: f64) -> Dyadic {
    Dyadic::from_f64(x)
}

fn cloud(n: usize, seed: u64, scale: f64) -> Vec<Point3> {
    let mut s = seed;
    let mut next = move || {
        s = s
            .wrapping_mul(6_364_136_223_846_793_005)
            .wrapping_add(1_442_695_040_888_963_407);
        ((s >> 11) as f64 / (1u64 << 53) as f64 - 0.5) * scale
    };
    (0..n).map(|_| p(next(), next(), next())).collect()
}

/// A rotation taking the world axes to `(e0, e1, e2)`.
fn rotation() -> [Vec3; 3] {
    let e0 = Vec3::new(1.0, 2.0, 2.0) / 3.0;
    let e1 = Vec3::new(2.0, 1.0, -2.0) / 3.0;
    let tilt = 0.3f64;
    let (a, b) = (
        e0 * tilt.cos() + e1 * tilt.sin(),
        e1 * tilt.cos() - e0 * tilt.sin(),
    );
    [a, b, a.cross(b)]
}

fn rotated(q: Point3, offset: Point3) -> Point3 {
    let [a, b, c] = rotation();
    offset + a * q.x + b * q.y + c * q.z
}

// ---------------------------------------------------------------------------
// Sphere
// ---------------------------------------------------------------------------

fn sphere_contains_all(points: &[Point3], m: &MinimumSphere) {
    let c = m.sphere.centre;
    let r = exact(m.sphere.radius);
    for q in points {
        let dx = exact(q.x).sub(&exact(c.x));
        let dy = exact(q.y).sub(&exact(c.y));
        let dz = exact(q.z).sub(&exact(c.z));
        let slack = r
            .mul(&r)
            .sub(&dx.mul(&dx).add(&dy.mul(&dy)).add(&dz.mul(&dz)));
        assert_ne!(slack.sign(), Some(Sign::Negative), "{q:?} outside {m:?}");
    }
}

/// Every sphere on two, three or four of the points that holds them all,
/// in floating point: an independent oracle for small sets.
fn brute_force(points: &[Point3]) -> f64 {
    let holds = |c: Point3, r: f64| {
        points
            .iter()
            .all(|q| (*q - c).length() <= r * (1.0 + 1e-10) + 1e-12)
    };
    let n = points.len();
    let mut best = f64::INFINITY;
    let mut consider = |c: Point3, r: f64| {
        if r.is_finite() && r < best && holds(c, r) {
            best = r;
        }
    };
    for i in 0..n {
        for j in i + 1..n {
            let (a, b) = (points[i], points[j]);
            consider((a + b) * 0.5, (b - a).length() * 0.5);
            for k in j + 1..n {
                let c = points[k];
                let (u, v) = (b - a, c - a);
                let w = u.cross(v);
                let ww = w.length_squared();
                if ww > 0.0 {
                    let centre = a
                        + (v.cross(w) * u.length_squared() + w.cross(u) * v.length_squared())
                            / (2.0 * ww);
                    consider(centre, (a - centre).length());
                }
                for &d in &points[k + 1..] {
                    let t = d - a;
                    let det = 2.0 * u.dot(v.cross(t));
                    if det != 0.0 {
                        let centre = a
                            + (v.cross(t) * u.length_squared()
                                + t.cross(u) * v.length_squared()
                                + u.cross(v) * t.length_squared())
                                / det;
                        consider(centre, (a - centre).length());
                    }
                }
            }
        }
    }
    best
}

#[test]
fn a_diameter_a_triangle_and_a_tetrahedron() {
    let pole = [
        p(0.0, 0.0, 0.0),
        p(0.0, 0.0, 6.0),
        p(1.0, 1.0, 3.0),
        p(-2.0, 0.0, 2.5),
    ];
    let m = minimum_enclosing_sphere(&pole).unwrap();
    assert_eq!(m.evidence.support, vec![0, 1]);
    assert_eq!(m.sphere.centre, p(0.0, 0.0, 3.0));
    assert!(m.sphere.radius >= 3.0 && m.sphere.radius - 3.0 <= m.evidence.error);
    sphere_contains_all(&pole, &m);

    // Centre (1/3, 1/3, 1/3), radius sqrt(2/3); the origin is inside.
    let corner = [
        p(1.0, 0.0, 0.0),
        p(0.0, 1.0, 0.0),
        p(0.0, 0.0, 1.0),
        p(0.0, 0.0, 0.0),
    ];
    let m = minimum_enclosing_sphere(&corner).unwrap();
    assert_eq!(m.evidence.support, vec![0, 1, 2]);
    let e = m.evidence.error;
    assert!(e < 1e-14, "{e}");
    assert!((m.sphere.centre - p(1.0, 1.0, 1.0) / 3.0).length() <= e);
    let r = (2.0f64 / 3.0).sqrt();
    assert!(m.sphere.radius >= r - 1e-16 && m.sphere.radius - r <= e + 1e-16);
    sphere_contains_all(&corner, &m);

    // A regular tetrahedron: centre 0, radius sqrt 3.
    let tetra = [
        p(1.0, 1.0, 1.0),
        p(1.0, -1.0, -1.0),
        p(-1.0, 1.0, -1.0),
        p(-1.0, -1.0, 1.0),
        p(0.2, -0.1, 0.3),
    ];
    let m = minimum_enclosing_sphere(&tetra).unwrap();
    assert_eq!(m.evidence.support, vec![0, 1, 2, 3]);
    assert!(m.sphere.centre.length() <= m.evidence.error);
    assert!(m.sphere.radius >= 3f64.sqrt() - 1e-16);
    assert!(m.sphere.radius - 3f64.sqrt() <= m.evidence.error + 1e-16);
    sphere_contains_all(&tetra, &m);
}

#[test]
fn points_on_a_known_sphere() {
    let (centre, radius) = (p(1.0, -2.0, 3.0), 4.0);
    let n = 400;
    let golden = std::f64::consts::PI * (3.0 - 5f64.sqrt());
    let mut points: Vec<Point3> = (0..n)
        .map(|i| {
            let z = 1.0 - 2.0 * (i as f64 + 0.5) / n as f64;
            let rho = (1.0 - z * z).sqrt();
            let t = golden * i as f64;
            centre + p(rho * t.cos(), rho * t.sin(), z) * radius
        })
        .collect();
    points.extend(cloud(300, 5, 4.0).into_iter().map(|q| centre + q));
    let m = minimum_enclosing_sphere(&points).unwrap();
    assert!((m.sphere.radius - radius).abs() < 1e-3, "{m:?}");
    assert!(m.sphere.radius <= radius + 1e-12);
    assert!(m.evidence.error < 1e-12);
    sphere_contains_all(&points, &m);
}

#[test]
fn duplicates_collinear_and_coplanar_points() {
    let same = [p(1.0, 2.0, 3.0); 4];
    let m = minimum_enclosing_sphere(&same).unwrap();
    assert_eq!(m.sphere.radius, 0.0);
    assert_eq!(m.evidence.error, 0.0);
    assert_eq!(m.sphere.centre, p(1.0, 2.0, 3.0));

    let line: Vec<Point3> = [0.0, 3.0, 1.0, 8.0, 8.0, 2.0]
        .iter()
        .map(|&t| p(t, 2.0 * t, -t))
        .collect();
    let m = minimum_enclosing_sphere(&line).unwrap();
    assert_eq!(m.evidence.support.len(), 2);
    assert_eq!(m.sphere.centre, p(4.0, 8.0, -4.0));
    let r = 4.0 * 6f64.sqrt();
    assert!(m.sphere.radius >= r - 1e-15 && m.sphere.radius - r <= m.evidence.error + 1e-15);

    // Coplanar points: the sphere is the plane's minimum circle.
    let flat: Vec<Point2> = cloud(40, 3, 5.0)
        .into_iter()
        .map(|q| Point2::new(q.x, q.y))
        .collect();
    let circle = axiolid_overlay::minimum_enclosing_circle(&flat).unwrap();
    let lifted: Vec<Point3> = flat.iter().map(|q| p(q.x, q.y, 0.0)).collect();
    let m = minimum_enclosing_sphere(&lifted).unwrap();
    assert!(m.evidence.support.len() <= 3);
    assert!(
        (m.sphere.radius - circle.circle.radius).abs() <= m.evidence.error + circle.evidence.error
    );
    sphere_contains_all(&lifted, &m);
}

#[test]
fn random_clouds_match_the_brute_force() {
    for seed in 1..=10u64 {
        let scale = 10f64.powi(seed as i32 % 5 - 2);
        let points = cloud(9 + seed as usize % 4, seed, scale);
        let m = minimum_enclosing_sphere(&points).unwrap();
        sphere_contains_all(&points, &m);
        let oracle = brute_force(&points);
        assert!(
            (m.sphere.radius - oracle).abs() <= 1e-9 * oracle,
            "seed {seed}: {} vs {oracle}",
            m.sphere.radius
        );
        assert!(m.evidence.error <= 1e-12 * oracle, "seed {seed}: {m:?}");
        for &i in &m.evidence.support {
            let d = (points[i] - m.sphere.centre).length();
            assert!((d - m.sphere.radius).abs() <= 2.0 * m.evidence.error + 1e-14 * oracle);
        }
    }
}

#[test]
fn larger_clouds_are_contained_and_order_free() {
    let points = cloud(2000, 17, 50.0);
    let m = minimum_enclosing_sphere(&points).unwrap();
    sphere_contains_all(&points, &m);
    let mut reversed = points.clone();
    reversed.reverse();
    let r = minimum_enclosing_sphere(&reversed).unwrap();
    let e = m.evidence.error + r.evidence.error;
    assert!((m.sphere.centre - r.sphere.centre).length() <= e);
    assert!((m.sphere.radius - r.sphere.radius).abs() <= e);
}

#[test]
fn empty_and_non_finite_input_are_refused() {
    assert_eq!(minimum_enclosing_sphere(&[]), Err(BoundingError::Empty));
    assert_eq!(oriented_bounding_box(&[]), Err(BoundingError::Empty));
    for bad in [f64::NAN, f64::INFINITY, f64::NEG_INFINITY] {
        let points = [p(0.0, 0.0, 0.0), p(1.0, bad, 0.0)];
        assert_eq!(
            minimum_enclosing_sphere(&points),
            Err(BoundingError::NonFinite)
        );
        assert_eq!(
            oriented_bounding_box(&points),
            Err(BoundingError::NonFinite)
        );
    }
}

// ---------------------------------------------------------------------------
// Oriented box
// ---------------------------------------------------------------------------

/// `|(p - centre) . axes[i]| <= half_extents[i]` for every point, exactly.
fn box_contains_all(points: &[Point3], b: &OrientedBox) {
    let c = b.centre;
    for q in points {
        for k in 0..3 {
            let a = b.axes[k];
            let t = exact(q.x)
                .sub(&exact(c.x))
                .mul(&exact(a.x))
                .add(&exact(q.y).sub(&exact(c.y)).mul(&exact(a.y)))
                .add(&exact(q.z).sub(&exact(c.z)).mul(&exact(a.z)));
            let h = exact(b.half_extents[k]);
            assert_ne!(
                h.sub(&t).sign(),
                Some(Sign::Negative),
                "{q:?} past axis {k}"
            );
            assert_ne!(
                h.add(&t).sign(),
                Some(Sign::Negative),
                "{q:?} past axis {k}"
            );
        }
    }
}

fn sorted(mut h: [f64; 3]) -> [f64; 3] {
    h.sort_by(f64::total_cmp);
    h
}

fn box_points(sides: [f64; 3], interior: usize) -> Vec<Point3> {
    let offset = p(4.0, -3.0, 7.0);
    let mut points: Vec<Point3> = (0..8)
        .map(|i| {
            let s = |k: usize| if i & (1 << k) == 0 { -0.5 } else { 0.5 };
            rotated(p(s(0) * sides[0], s(1) * sides[1], s(2) * sides[2]), offset)
        })
        .collect();
    points.extend(
        cloud(interior, 11, 1.0)
            .into_iter()
            .map(|q| rotated(p(q.x * sides[0], q.y * sides[1], q.z * sides[2]), offset)),
    );
    points
}

#[test]
fn a_rotated_box_is_recovered() {
    for sides in [[1.0, 2.0, 3.0], [1.0, 1.0, 1.0], [0.5, 4.0, 4.0]] {
        let points = box_points(sides, 50);
        let obb = oriented_bounding_box(&points).unwrap();
        let b = obb.bounding_box;
        let volume = sides[0] * sides[1] * sides[2];
        assert!(
            (b.volume() - volume).abs() <= 1e-12 * volume,
            "{sides:?}: {} vs {volume}",
            b.volume()
        );
        let h = sorted(b.half_extents);
        let expected = sorted(sides.map(|s| s * 0.5));
        for k in 0..3 {
            assert!((h[k] - expected[k]).abs() <= 1e-12, "{h:?} vs {expected:?}");
        }
        assert!((b.centre - p(4.0, -3.0, 7.0)).length() <= 1e-12);
        assert!(b.volume() <= obb.evidence.axis_aligned_volume);
        assert!(obb.evidence.axis_aligned_volume > 1.5 * volume);
        assert!(obb.evidence.orthogonality <= 1e-15);
        box_contains_all(&points, &b);
    }
}

#[test]
fn an_axis_aligned_box_keeps_the_world_axes() {
    let points: Vec<Point3> = (0..8)
        .map(|i| {
            let s = |k: usize, lo: f64, hi: f64| if i & (1 << k) == 0 { lo } else { hi };
            p(s(0, -1.0, 3.0), s(1, 0.0, 2.0), s(2, 5.0, 6.0))
        })
        .collect();
    let obb = oriented_bounding_box(&points).unwrap();
    let b = obb.bounding_box;
    assert!((b.volume() - 8.0).abs() <= 1e-12);
    assert_eq!(obb.evidence.orthogonality, 0.0);
    box_contains_all(&points, &b);
}

#[test]
fn flat_and_degenerate_input() {
    // A 2 x 3 rectangle in a tilted plane: one extent vanishes.
    let flat: Vec<Point3> = [(0.0, 0.0), (2.0, 0.0), (2.0, 3.0), (0.0, 3.0), (1.0, 1.0)]
        .iter()
        .map(|&(x, y)| rotated(p(x, y, 0.0), p(1.0, 1.0, 1.0)))
        .collect();
    let b = oriented_bounding_box(&flat).unwrap().bounding_box;
    let h = sorted(b.half_extents);
    assert!(h[0] <= 1e-14, "{h:?}");
    assert!(
        (h[1] - 1.0).abs() <= 1e-12 && (h[2] - 1.5).abs() <= 1e-12,
        "{h:?}"
    );
    box_contains_all(&flat, &b);

    let line: Vec<Point3> = [0.0, 1.0, 5.0, 2.0]
        .iter()
        .map(|&t| rotated(p(t, 0.0, 0.0), p(0.0, 0.0, 0.0)))
        .collect();
    let b = oriented_bounding_box(&line).unwrap().bounding_box;
    let h = sorted(b.half_extents);
    assert!(h[0] <= 1e-14 && h[1] <= 1e-14, "{h:?}");
    assert!((h[2] - 2.5).abs() <= 1e-12, "{h:?}");
    box_contains_all(&line, &b);

    let one = [p(2.0, 3.0, 4.0); 3];
    let b = oriented_bounding_box(&one).unwrap().bounding_box;
    assert_eq!(b.half_extents, [0.0; 3]);
    assert_eq!(b.centre, p(2.0, 3.0, 4.0));
    box_contains_all(&one, &b);
}

#[test]
fn random_clouds_are_contained_and_no_worse_than_axis_aligned() {
    for seed in 1..=8u64 {
        let raw = cloud(300, seed, 10.0);
        // Squash and turn the cloud so the axis-aligned box is poor.
        let points: Vec<Point3> = raw
            .iter()
            .map(|q| rotated(p(q.x, 0.2 * q.y, 0.05 * q.z), p(1.0, 2.0, 3.0)))
            .collect();
        let obb = oriented_bounding_box(&points).unwrap();
        let b = obb.bounding_box;
        box_contains_all(&points, &b);
        assert!(b.volume() <= obb.evidence.axis_aligned_volume);
        assert!(
            b.volume() < 0.5 * obb.evidence.axis_aligned_volume,
            "seed {seed}"
        );
        assert!(obb.evidence.orthogonality <= 1e-15);
        assert!(obb.evidence.candidates > 6);
        // Every corner of the box is a point of it.
        for corner in b.corners() {
            assert!(
                (corner - b.centre).length()
                    <= b.half_extents.iter().map(|h| h * h).sum::<f64>().sqrt() * (1.0 + 1e-12)
            );
        }
    }
}
