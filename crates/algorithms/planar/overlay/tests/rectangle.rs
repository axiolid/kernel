//! Minimum-area rectangles (#182), against a brute-force oracle.

use axiolid_core::{Point2, Vec2};
use axiolid_overlay::{minimum_area_rectangle, RectangleError};

fn p(x: f64, y: f64) -> Point2 {
    Point2::new(x, y)
}

/// The least area over many directions, a lower bound on no rectangle but
/// an upper bound on the minimum, tight to the sampling.
fn sampled_minimum(points: &[Point2], samples: usize) -> f64 {
    (0..samples)
        .map(|k| {
            let a = std::f64::consts::FRAC_PI_2 * k as f64 / samples as f64;
            let (u, v) = (Vec2::new(a.cos(), a.sin()), Vec2::new(-a.sin(), a.cos()));
            let span = |axis: Vec2| {
                let t = points.iter().map(|q| q.dot(axis));
                t.clone().fold(f64::NEG_INFINITY, f64::max) - t.fold(f64::INFINITY, f64::min)
            };
            span(u) * span(v)
        })
        .fold(f64::INFINITY, f64::min)
}

/// Every point is inside the rectangle, to its stated error.
fn encloses(points: &[Point2], r: &axiolid_overlay::OrientedRectangle, error: f64) {
    for q in points {
        let d = *q - r.centre;
        for k in 0..2 {
            assert!(
                d.dot(r.axes[k]).abs() <= r.half_extents[k] + 2.0 * error,
                "{q:?} outside {r:?}"
            );
        }
    }
}

#[test]
fn a_rotated_rectangle_is_its_own_minimum() {
    let (a, b) = (30f64.to_radians(), 120f64.to_radians());
    let (u, v) = (Vec2::new(a.cos(), a.sin()), Vec2::new(b.cos(), b.sin()));
    let c = p(5.0, -2.0);
    let corners = [
        c - 3.0 * u - v,
        c + 3.0 * u - v,
        c + 3.0 * u + v,
        c - 3.0 * u + v,
    ];
    let mut points = corners.to_vec();
    points.push(c);
    points.push(c + 0.5 * u);
    let m = minimum_area_rectangle(&points).unwrap();
    let r = m.rectangle;
    let e = m.evidence.error;
    assert_eq!(m.evidence.hull_vertices, 4);
    assert_eq!(m.evidence.minimal_orientations, 1);
    assert!((r.axes[0] - u).length() < 1e-15);
    assert!((r.half_extents[0] - 3.0).abs() <= 1e-14 && (r.half_extents[1] - 1.0).abs() <= 1e-14);
    assert!((r.centre - c).length() <= 1e-14);
    assert!(e > 0.0 && e < 1e-12);
    // Corners within the stated error of the input's own (which were
    // themselves rounded, by less than an ulp of 8).
    for (x, y) in r.corners().iter().zip(corners) {
        assert!(
            (*x - y).length() <= e + 8.0 * f64::EPSILON * 8.0,
            "{x:?} vs {y:?}"
        );
    }
}

#[test]
fn a_convex_polygon_matches_the_sampled_minimum() {
    let points: Vec<Point2> = (0..11)
        .map(|k| {
            let a = k as f64 * 0.57 + 0.1;
            p(
                (2.0 + (k % 3) as f64) * a.cos(),
                (1.0 + 0.3 * k as f64) * a.sin(),
            )
        })
        .collect();
    let m = minimum_area_rectangle(&points).unwrap();
    encloses(&points, &m.rectangle, m.evidence.error);
    let sampled = sampled_minimum(&points, 200_000);
    assert!(
        m.rectangle.area() <= sampled + 1e-9,
        "{} > {sampled}",
        m.rectangle.area()
    );
    assert!(m.rectangle.area() >= sampled - 1e-4 * sampled);
}

#[test]
fn every_hull_edge_is_tried_not_just_the_first_few() {
    // A long thin diamond whose minimum is on the last hull edges walked:
    // skipping a caliper step or an edge gives a larger area.
    let points = [
        p(0.0, 0.0),
        p(10.0, 1.0),
        p(11.0, 3.0),
        p(1.0, 2.0),
        p(4.0, 1.5),
    ];
    let m = minimum_area_rectangle(&points).unwrap();
    encloses(&points, &m.rectangle, m.evidence.error);
    let sampled = sampled_minimum(&points, 400_000);
    assert!(
        (m.rectangle.area() - sampled).abs() <= 1e-4 * sampled,
        "{} vs {sampled}",
        m.rectangle.area()
    );
}

#[test]
fn random_sets_match_the_sampled_minimum_and_are_enclosed() {
    let mut s = 0x2545_f491_4f6c_dd1du64;
    let mut next = || {
        s ^= s << 13;
        s ^= s >> 7;
        s ^= s << 17;
        (s >> 11) as f64 / (1u64 << 53) as f64
    };
    for _ in 0..40 {
        let n = 3 + (next() * 30.0) as usize;
        let points: Vec<Point2> = (0..n)
            .map(|_| p(next() * 10.0 - 5.0, next() * 4.0 - 1.0))
            .collect();
        let m = minimum_area_rectangle(&points).unwrap();
        encloses(&points, &m.rectangle, m.evidence.error);
        let sampled = sampled_minimum(&points, 20_000);
        let area = m.rectangle.area();
        assert!(area <= sampled + 1e-9, "{area} > {sampled}");
        assert!(area >= sampled * (1.0 - 2e-3), "{area} < {sampled}");
        // Far from the origin the output rounds more, within its bound.
        let far: Vec<Point2> = points.iter().map(|q| *q + Vec2::new(1e6, -2e6)).collect();
        let m = minimum_area_rectangle(&far).unwrap();
        assert!(m.evidence.error > 1e-9 && m.evidence.error < 1e-7);
        encloses(&far, &m.rectangle, m.evidence.error);
        assert!((m.rectangle.area() - area).abs() <= 1e-6 * area.max(1.0));
    }
}

#[test]
fn collinear_points_give_a_zero_width_rectangle() {
    let points = [p(1.0, 1.0), p(3.0, 2.0), p(-1.0, 0.0), p(5.0, 3.0)];
    let m = minimum_area_rectangle(&points).unwrap();
    let r = m.rectangle;
    assert_eq!(m.evidence.hull_vertices, 2);
    assert_eq!(r.area(), 0.0);
    assert!((r.centre - p(2.0, 1.5)).length() <= 1e-15);
    assert!((r.half_extents[0] - 0.5 * 45f64.sqrt()).abs() < 1e-14);
    encloses(&points, &r, m.evidence.error);
    // A vertical segment's first axis is turned into [0, 90) degrees, so
    // the segment runs along the second.
    let v = minimum_area_rectangle(&[p(0.0, 2.0), p(0.0, 0.0)])
        .unwrap()
        .rectangle;
    assert_eq!(v.axes[0], Vec2::X);
    assert_eq!(v.half_extents, [0.0, 1.0]);
    assert_eq!(v.centre, p(0.0, 1.0));
}

#[test]
fn a_single_point_and_no_points() {
    let m = minimum_area_rectangle(&[p(2.0, 3.0), p(2.0, 3.0)]).unwrap();
    assert_eq!(m.rectangle.centre, p(2.0, 3.0));
    assert_eq!(m.rectangle.half_extents, [0.0, 0.0]);
    assert_eq!(m.evidence.error, 0.0);
    assert_eq!(minimum_area_rectangle(&[]), Err(RectangleError::Empty));
    assert_eq!(
        minimum_area_rectangle(&[p(f64::NAN, 0.0)]),
        Err(RectangleError::NonFinite)
    );
}

#[test]
fn ties_resolve_the_same_way_for_any_input_order() {
    // A square's four edges give one orientation.
    let square = [p(0.0, 0.0), p(2.0, 0.0), p(2.0, 2.0), p(0.0, 2.0)];
    let m = minimum_area_rectangle(&square).unwrap();
    assert_eq!(m.evidence.minimal_orientations, 1);
    assert_eq!(m.rectangle.axes[0], Vec2::X);
    // This triangle's boxes at 0 and 45 degrees both have area 2, exactly.
    // The hull walk meets 45 degrees first; 0 degrees must still win.
    let triangle = [
        p(0.0, 0.0),
        p(1.0, -1.0),
        p(1.0, 1.0),
        p(0.5, 0.0),
        p(0.9, 0.2),
    ];
    let first = minimum_area_rectangle(&triangle).unwrap();
    assert_eq!(first.evidence.minimal_orientations, 2, "{first:?}");
    assert_eq!(first.rectangle.axes[0], Vec2::X);
    assert_eq!(first.rectangle.half_extents, [0.5, 1.0]);
    for shift in 1..5 {
        let mut turned = triangle.to_vec();
        turned.rotate_left(shift);
        if shift % 2 == 1 {
            turned.reverse();
        }
        assert_eq!(minimum_area_rectangle(&turned).unwrap(), first);
    }
}

#[test]
fn an_axis_aligned_box_reports_the_rounding_actually_done() {
    // A 2.5 m x 5 m bay: nothing is rounded, so a 5 m side is exactly 5 m.
    let bay = [
        p(10.0, 20.0),
        p(12.5, 20.0),
        p(12.5, 25.0),
        p(10.0, 25.0),
        p(11.0, 22.0),
    ];
    let m = minimum_area_rectangle(&bay).unwrap();
    assert_eq!(m.rectangle.axes, [Vec2::X, Vec2::Y]);
    assert_eq!(m.rectangle.half_extents, [1.25, 2.5]);
    assert_eq!(m.rectangle.centre, p(11.25, 22.5));
    assert_eq!(m.evidence.error, 0.0);
    // Decimal coordinates round: the error is what was rounded, no more,
    // and it covers the corners.
    let decimal = [p(0.1, 0.3), p(2.6, 0.3), p(2.6, 5.3), p(0.1, 5.3)];
    let m = minimum_area_rectangle(&decimal).unwrap();
    let e = m.evidence.error;
    assert!(e < 4.0 * f64::EPSILON * 5.3, "{e}");
    let [c0, _, c2, _] = m.rectangle.corners();
    assert!(
        (c0.x - 0.1).abs() <= e
            && (c0.y - 0.3).abs() <= e
            && (c2.x - 2.6).abs() <= e
            && (c2.y - 5.3).abs() <= e
    );
    // Along an axis, collinear: exact too.
    let line = [p(0.0, 1.0), p(3.0, 1.0), p(1.0, 1.0)];
    assert_eq!(minimum_area_rectangle(&line).unwrap().evidence.error, 0.0);
}

#[test]
fn the_measured_error_covers_every_corner() {
    // From x = 1.3 to 4.22 the low corner, centre minus half extent,
    // rounds twice as far as the centre or the half extent do.
    let m =
        minimum_area_rectangle(&[p(1.3, 0.0), p(4.22, 0.0), p(4.22, 1.0), p(1.3, 1.0)]).unwrap();
    let e = m.evidence.error;
    assert!(e > 0.0);
    let exact = [p(1.3, 0.0), p(4.22, 0.0), p(4.22, 1.0), p(1.3, 1.0)];
    for (c, x) in m.rectangle.corners().iter().zip(exact) {
        assert!(
            (c.x - x.x).abs() <= e && (c.y - x.y).abs() <= e,
            "{c:?} vs {x:?}, {e}"
        );
    }
}
