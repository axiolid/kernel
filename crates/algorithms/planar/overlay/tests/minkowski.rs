//! Minkowski sums and erosions by a convex polygon (#145), and disc
//! morphology on a stated side (#163).
//!
//! Oracles written apart from the construction:
//! - a point is in `R + K` exactly when the polygon `x - K` meets `R`, and
//!   in `R - K` exactly when `x + K` lies within `R`, decided here by
//!   point-in-polygon and segment crossings on random points kept away
//!   from the result's boundary;
//! - for a disc, a point is in the exact dilation when its distance to `R`
//!   is at most `r`, and in the exact erosion when it lies in `R` at least
//!   `r` from its boundary; for convex polygons the areas are
//!   `A + L r + pi r^2` and `A - L r + r^2 sum cot(theta_i / 2)`.

use axiolid_core::{Point2, Tolerance};
use axiolid_overlay::{BoundSide, MinkowskiError, Polygon, Region, Ring};

fn tol() -> Tolerance {
    Tolerance::METRE
}

fn ring(points: &[(f64, f64)]) -> Ring {
    Ring {
        points: points.iter().map(|&(x, y)| Point2::new(x, y)).collect(),
    }
}

fn region(polygons: Vec<Polygon>) -> Region {
    Region::new(polygons, tol()).expect("a valid region")
}

/// A U-shaped block with a square hole in its base: non-convex, with a hole.
fn u_with_hole() -> Region {
    region(vec![Polygon {
        outer: ring(&[
            (0.0, 0.0),
            (6.0, 0.0),
            (6.0, 5.0),
            (4.0, 5.0),
            (4.0, 2.0),
            (2.0, 2.0),
            (2.0, 5.0),
            (0.0, 5.0),
        ]),
        holes: vec![ring(&[(0.5, 0.5), (1.5, 0.5), (1.5, 1.5), (0.5, 1.5)])],
    }])
}

/// A rectangle `w x h` about the origin, turned by `angle`.
fn rectangle(w: f64, h: f64, angle: f64) -> Vec<Point2> {
    let (s, c) = angle.sin_cos();
    [(-w, -h), (w, -h), (w, h), (-w, h)]
        .iter()
        .map(|&(x, y)| Point2::new(0.5 * (c * x - s * y), 0.5 * (s * x + c * y)))
        .collect()
}

fn in_ring(points: &[Point2], p: Point2) -> bool {
    let mut inside = false;
    let n = points.len();
    for i in 0..n {
        let (a, b) = (points[i], points[(i + 1) % n]);
        if (a.y > p.y) != (b.y > p.y) && p.x < (b.x - a.x) * (p.y - a.y) / (b.y - a.y) + a.x {
            inside = !inside;
        }
    }
    inside
}

fn in_region(region: &Region, p: Point2) -> bool {
    region.polygons().iter().any(|poly| {
        in_ring(&poly.outer.points, p) && !poly.holes.iter().any(|h| in_ring(&h.points, p))
    })
}

/// Edges of a region.
fn edges(region: &Region) -> Vec<(Point2, Point2)> {
    region
        .boundary_rings()
        .iter()
        .flat_map(|r| {
            let n = r.points.len();
            (0..n).map(move |i| (r.points[i], r.points[(i + 1) % n]))
        })
        .collect()
}

fn cross(a: Point2, b: Point2, c: Point2) -> f64 {
    (b - a).perp_dot(c - a)
}

fn segments_cross(a: Point2, b: Point2, c: Point2, d: Point2) -> bool {
    let (d1, d2) = (cross(a, b, c), cross(a, b, d));
    let (d3, d4) = (cross(c, d, a), cross(c, d, b));
    (d1 > 0.0) != (d2 > 0.0) && (d3 > 0.0) != (d4 > 0.0)
}

/// Distance from `p` to a segment.
fn to_segment(p: Point2, a: Point2, b: Point2) -> f64 {
    let d = b - a;
    let t = ((p - a).dot(d) / d.dot(d)).clamp(0.0, 1.0);
    (a + d * t - p).length()
}

fn to_boundary(region: &Region, p: Point2) -> f64 {
    edges(region)
        .iter()
        .map(|&(a, b)| to_segment(p, a, b))
        .fold(f64::INFINITY, f64::min)
}

/// Whether the convex polygon `k` meets `region` (boundaries or either
/// inside the other).
fn meets(k: &[Point2], region: &Region) -> bool {
    let n = k.len();
    for &(a, b) in &edges(region) {
        for i in 0..n {
            if segments_cross(a, b, k[i], k[(i + 1) % n]) {
                return true;
            }
        }
        if in_ring(k, a) {
            return true;
        }
    }
    k.iter().any(|&q| in_region(region, q))
}

/// Whether the convex polygon `k` lies within `region`.
fn within(k: &[Point2], region: &Region) -> bool {
    let n = k.len();
    k.iter().all(|&q| in_region(region, q))
        && !edges(region).iter().any(|&(a, b)| {
            (0..n).any(|i| segments_cross(a, b, k[i], k[(i + 1) % n])) || in_ring(k, a)
        })
}

/// Deterministic points over a box.
fn samples(lo: Point2, hi: Point2, count: usize) -> Vec<Point2> {
    let mut state = 0x9e37_79b9_7f4a_7c15_u64;
    let mut next = || {
        state = state
            .wrapping_mul(6_364_136_223_846_793_005)
            .wrapping_add(1_442_695_040_888_963_407);
        (state >> 11) as f64 / (1u64 << 53) as f64
    };
    (0..count)
        .map(|_| Point2::new(lo.x + (hi.x - lo.x) * next(), lo.y + (hi.y - lo.y) * next()))
        .collect()
}

#[test]
fn an_axis_rectangle_widens_a_frame_by_its_half_sides() {
    // [0, 4]^2 with the hole [1, 3]^2, plus [-1/4, 1/4] x [-1/10, 1/10]:
    // the outside grows to [-1/4, 17/4] x [-1/10, 41/10] and the hole
    // shrinks to [5/4, 11/4] x [11/10, 29/10].
    let frame = region(vec![Polygon {
        outer: ring(&[(0.0, 0.0), (4.0, 0.0), (4.0, 4.0), (0.0, 4.0)]),
        holes: vec![ring(&[(1.0, 1.0), (3.0, 1.0), (3.0, 3.0), (1.0, 3.0)])],
    }]);
    let k = ring(&[(-0.25, -0.1), (0.25, -0.1), (0.25, 0.1), (-0.25, 0.1)]);
    let sum = frame.minkowski_sum(&k, tol()).unwrap();
    assert_eq!(sum.component_count(), 1);
    assert_eq!(sum.polygons()[0].holes.len(), 1);
    let want = 4.5 * 4.2 - 1.5 * 1.8;
    assert!(
        (sum.area() - want).abs() < 1e-12,
        "{} vs {want}",
        sum.area()
    );
    // And eroding it back by the same rectangle recovers the frame.
    let back = sum.minkowski_erosion(&k, tol()).unwrap();
    assert!(
        (back.area() - frame.area()).abs() < 1e-12,
        "{}",
        back.area()
    );
}

#[test]
fn a_turned_rectangle_summed_with_a_non_convex_region_with_a_hole() {
    let r = u_with_hole();
    for angle in [0.0, 0.3, 1.1] {
        let k = rectangle(0.8, 0.3, angle);
        let sum = r.minkowski_sum(&Ring { points: k.clone() }, tol()).unwrap();
        let mut checked = 0;
        for p in samples(Point2::new(-1.5, -1.5), Point2::new(7.5, 6.5), 3000) {
            if to_boundary(&sum, p) < 1e-6 {
                continue;
            }
            let moved: Vec<Point2> = k
                .iter()
                .map(|q| Point2::new(p.x - q.x, p.y - q.y))
                .collect();
            assert_eq!(
                in_region(&sum, p),
                meets(&moved, &r),
                "angle {angle}, point {p:?}"
            );
            checked += 1;
        }
        assert!(checked > 2500);
    }
}

#[test]
fn a_turned_rectangle_eroding_a_non_convex_region_with_a_hole() {
    let r = u_with_hole();
    for angle in [0.0, 0.4] {
        let k = rectangle(0.8, 0.3, angle);
        let eroded = r
            .minkowski_erosion(&Ring { points: k.clone() }, tol())
            .unwrap();
        assert!(!eroded.is_empty());
        for p in samples(Point2::new(-0.5, -0.5), Point2::new(6.5, 5.5), 3000) {
            if to_boundary(&eroded, p) < 1e-6 {
                continue;
            }
            let moved: Vec<Point2> = k
                .iter()
                .map(|q| Point2::new(p.x + q.x, p.y + q.y))
                .collect();
            assert_eq!(
                in_region(&eroded, p),
                within(&moved, &r),
                "angle {angle}, point {p:?}"
            );
        }
    }
}

#[test]
fn a_non_convex_structuring_polygon_is_refused_by_name() {
    let dart = ring(&[(0.0, 0.0), (2.0, 1.0), (0.0, 2.0), (0.5, 1.0)]);
    assert_eq!(
        u_with_hole().minkowski_sum(&dart, tol()),
        Err(MinkowskiError::NotConvex)
    );
}

/// The exact disc dilation and erosion, pointwise.
fn in_exact_dilation(r: &Region, p: Point2, radius: f64) -> bool {
    in_region(r, p) || to_boundary(r, p) <= radius
}

fn in_exact_erosion(r: &Region, p: Point2, radius: f64) -> bool {
    in_region(r, p) && to_boundary(r, p) >= radius
}

#[test]
fn inner_and_outer_discs_bracket_the_exact_morphology() {
    let convex = region(vec![Polygon {
        outer: ring(&[(0.0, 0.0), (5.0, 0.0), (6.0, 3.0), (2.0, 4.0)]),
        holes: vec![],
    }]);
    for r in [convex, u_with_hole()] {
        for radius in [0.2, 0.45] {
            let di = r.dilate_inner(radius, tol()).unwrap();
            let dout = r.dilate_outer(radius, tol()).unwrap();
            let ei = r.erode_inner(radius, tol()).unwrap();
            let eo = r.erode_outer(radius, tol()).unwrap();
            assert_eq!(di.bound().unwrap().side, BoundSide::Inner);
            assert_eq!(dout.bound().unwrap().side, BoundSide::Outer);
            assert_eq!(ei.bound().unwrap().side, BoundSide::Inner);
            assert_eq!(eo.bound().unwrap().side, BoundSide::Outer);
            for b in [di.bound(), dout.bound(), ei.bound(), eo.bound()] {
                // The 64-gon's deviation, about 0.12% of the radius, plus
                // the margin.
                let d = b.unwrap().deviation;
                assert!(d > 0.0 && d < 2e-3 * radius + 1e-5, "{d}");
            }
            assert!(r.dilate(radius, tol()).unwrap().bound().is_none());
            for p in samples(Point2::new(-1.0, -1.0), Point2::new(7.0, 6.0), 4000) {
                // Inner results lie within the exact one, the exact one
                // within the outer results.
                if in_region(&di, p) {
                    assert!(in_exact_dilation(&r, p, radius), "dilate_inner {p:?}");
                }
                if in_exact_dilation(&r, p, radius) {
                    assert!(in_region(&dout, p), "dilate_outer {p:?}");
                }
                if in_region(&ei, p) {
                    assert!(in_exact_erosion(&r, p, radius), "erode_inner {p:?}");
                }
                if in_exact_erosion(&r, p, radius) {
                    assert!(in_region(&eo, p), "erode_outer {p:?}");
                }
            }
            assert!(di.area() <= dout.area() && ei.area() <= eo.area());
        }
    }
}

#[test]
fn inner_and_outer_areas_bracket_the_exact_areas_on_convex_input() {
    // A convex quadrilateral: exact dilation area A + L r + pi r^2, exact
    // erosion area A - L r + r^2 sum cot(theta_i / 2) (for r below where
    // an edge vanishes).
    let pts = [(0.0, 0.0), (5.0, 0.0), (6.0, 3.0), (2.0, 4.0)];
    let convex = region(vec![Polygon {
        outer: ring(&pts),
        holes: vec![],
    }]);
    let p: Vec<Point2> = pts.iter().map(|&(x, y)| Point2::new(x, y)).collect();
    let n = p.len();
    let perimeter: f64 = (0..n).map(|i| (p[(i + 1) % n] - p[i]).length()).sum();
    let cot_sum: f64 = (0..n)
        .map(|i| {
            let (a, b, c) = (p[(i + n - 1) % n], p[i], p[(i + 1) % n]);
            let (u, v) = ((a - b).normalize(), (c - b).normalize());
            let theta = u.dot(v).clamp(-1.0, 1.0).acos();
            1.0 / (0.5 * theta).tan()
        })
        .sum();
    let area = convex.area();
    for radius in [0.1, 0.3, 0.6] {
        let dilation = area + perimeter * radius + core::f64::consts::PI * radius * radius;
        let erosion = area - perimeter * radius + radius * radius * cot_sum;
        let (di, dout) = (
            convex.dilate_inner(radius, tol()).unwrap(),
            convex.dilate_outer(radius, tol()).unwrap(),
        );
        let (ei, eo) = (
            convex.erode_inner(radius, tol()).unwrap(),
            convex.erode_outer(radius, tol()).unwrap(),
        );
        assert!(
            di.area() < dilation && dilation < dout.area(),
            "{} {dilation} {}",
            di.area(),
            dout.area()
        );
        assert!(
            ei.area() < erosion && erosion < eo.area(),
            "{} {erosion} {}",
            ei.area(),
            eo.area()
        );
        // And close: within the deviation times the boundary's length.
        let slack = 2.0
            * di.bound()
                .unwrap()
                .deviation
                .max(dout.bound().unwrap().deviation)
            * (perimeter + 7.0 * radius);
        assert!(
            dout.area() - di.area() < slack,
            "{} {}",
            dout.area() - di.area(),
            slack
        );
    }
}
