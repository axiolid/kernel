//! Minkowski sums and erosions by a polygon, convex or not, and by a
//! region (#145), and disc morphology on a stated side (#163).
//!
//! Oracles written apart from the construction:
//! - a point is in `R + K` exactly when the polygon `x - K` meets `R`, and
//!   in `R - K` exactly when `x + K` lies within `R`, decided here by
//!   point-in-polygon and segment crossings on random points kept away
//!   from the result's boundary (for non-convex operands: the boundaries
//!   cross or a vertex of one lies in the other);
//! - for a disc, a point is in the exact dilation when its distance to `R`
//!   is at most `r`, and in the exact erosion when it lies in `R` at least
//!   `r` from its boundary; for convex polygons the areas are
//!   `A + L r + pi r^2` and `A - L r + r^2 sum cot(theta_i / 2)`.

use axiolid_core::Vec2;
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

// ---------------------------------------------------------------------
// Non-convex structuring polygons and regions (#145).
// ---------------------------------------------------------------------

/// The L-shape `[0, 2s] x [0, s]  union  [0, s] x [0, 2s]`.
fn l_shape(s: f64) -> Ring {
    ring(&[
        (0.0, 0.0),
        (2.0 * s, 0.0),
        (2.0 * s, s),
        (s, s),
        (s, 2.0 * s),
        (0.0, 2.0 * s),
    ])
}

/// An arrowhead pointing right, reflex at its back vertex, turned by
/// `angle` about its middle and moved to `centre`.
fn dart(angle: f64, centre: (f64, f64)) -> Ring {
    let (s, c) = angle.sin_cos();
    Ring {
        points: [(0.0, 0.0), (2.0, 1.0), (0.0, 2.0), (0.5, 1.0)]
            .iter()
            .map(|&(x, y)| {
                let (x, y) = (x - 1.0, y - 1.0);
                Point2::new(c * x - s * y + centre.0, s * x + c * y + centre.1)
            })
            .collect(),
    }
}

fn single(r: &Ring) -> Region {
    region(vec![Polygon {
        outer: r.clone(),
        holes: vec![],
    }])
}

/// The region reflected and moved to `x` (`x - B`), or moved by `x`
/// (`x + B`).
fn placed(b: &Region, x: Point2, reflect: bool) -> Region {
    let s = if reflect { -1.0 } else { 1.0 };
    let map = |r: &Ring| Ring {
        points: r
            .points
            .iter()
            .map(|q| Point2::new(x.x + s * q.x, x.y + s * q.y))
            .collect(),
    };
    Region::new(
        b.polygons()
            .iter()
            .map(|p| Polygon {
                outer: map(&p.outer),
                holes: p.holes.iter().map(map).collect(),
            })
            .collect(),
        tol(),
    )
    .unwrap()
}

/// Whether two regions meet: their boundaries cross, or a vertex of one
/// lies in the other (if the boundaries do not cross, an outer ring of
/// one lies inside the other). Exact enough away from touching.
fn regions_meet(a: &Region, b: &Region) -> bool {
    let (ea, eb) = (edges(a), edges(b));
    ea.iter()
        .any(|&(p, q)| eb.iter().any(|&(r, s)| segments_cross(p, q, r, s)))
        || ea.iter().any(|&(p, _)| in_region(b, p))
        || eb.iter().any(|&(p, _)| in_region(a, p))
}

/// Whether `inner` lies within `outer`: no crossing, every vertex of
/// `inner` in `outer`, no vertex of `outer` in `inner`.
fn region_within(inner: &Region, outer: &Region) -> bool {
    let (ei, eo) = (edges(inner), edges(outer));
    !ei.iter()
        .any(|&(p, q)| eo.iter().any(|&(r, s)| segments_cross(p, q, r, s)))
        && ei.iter().all(|&(p, _)| in_region(outer, p))
        && !eo.iter().any(|&(p, _)| in_region(inner, p))
}

/// Check a sum and an erosion against the definitions on sampled points
/// clear of the results' boundaries.
fn check_by_definition(
    r: &Region,
    k: &Region,
    sum: &Region,
    eroded: &Region,
    lo: Point2,
    hi: Point2,
) {
    let mut checked = (0, 0);
    for p in samples(lo, hi, 1500) {
        if to_boundary(sum, p) > 1e-6 {
            assert_eq!(
                in_region(sum, p),
                regions_meet(&placed(k, p, true), r),
                "sum at {p:?}"
            );
            checked.0 += 1;
        }
        if to_boundary(eroded, p) > 1e-6 {
            assert_eq!(
                in_region(eroded, p),
                region_within(&placed(k, p, false), r),
                "erosion at {p:?}"
            );
            checked.1 += 1;
        }
    }
    assert!(checked.0 > 1400 && checked.1 > 1400, "{checked:?}");
}

/// The area of the symmetric difference.
fn apart(a: &Region, b: &Region) -> f64 {
    a.difference(b, tol()).unwrap().area() + b.difference(a, tol()).unwrap().area()
}

/// The corners of a region, sorted: vertices where the boundary turns
/// (the result keeps vertices where pieces met on a straight edge). The
/// turn is exact for the small dyadic coordinates used here.
fn vertices(region: &Region) -> Vec<(f64, f64)> {
    let mut v: Vec<(f64, f64)> = region
        .boundary_rings()
        .iter()
        .flat_map(|r| {
            let n = r.points.len();
            (0..n)
                .filter(move |&i| {
                    cross(
                        r.points[(i + n - 1) % n],
                        r.points[i],
                        r.points[(i + 1) % n],
                    ) != 0.0
                })
                .map(move |i| (r.points[i].x, r.points[i].y))
        })
        .collect();
    v.sort_by(|a, b| a.0.total_cmp(&b.0).then(a.1.total_cmp(&b.1)));
    v
}

#[test]
fn a_rectangle_plus_an_l_shape_is_an_l_shape() {
    // [0, 4] x [0, 1] + L(1) = [0, 6] x [0, 2]  union  [0, 5] x [0, 3].
    let r = single(&ring(&[(0.0, 0.0), (4.0, 0.0), (4.0, 1.0), (0.0, 1.0)]));
    let sum = r.minkowski_sum(&l_shape(1.0), tol()).unwrap();
    assert_eq!(sum.area(), 17.0);
    assert_eq!(
        vertices(&sum),
        vec![
            (0.0, 0.0),
            (0.0, 3.0),
            (5.0, 2.0),
            (5.0, 3.0),
            (6.0, 0.0),
            (6.0, 2.0)
        ]
    );
}

#[test]
fn an_l_shape_plus_itself_is_a_staircase() {
    // The pieces' sums are [0, 4] x [0, 2], [0, 3]^2 and [0, 2] x [0, 4].
    let l = l_shape(1.0);
    let sum = single(&l).minkowski_sum(&l, tol()).unwrap();
    assert_eq!(sum.area(), 13.0);
    assert_eq!(
        vertices(&sum),
        vec![
            (0.0, 0.0),
            (0.0, 4.0),
            (2.0, 3.0),
            (2.0, 4.0),
            (3.0, 2.0),
            (3.0, 3.0),
            (4.0, 0.0),
            (4.0, 2.0)
        ]
    );
}

#[test]
fn a_hole_stays_open_or_closes_with_the_operand_size() {
    // [0, 10]^2 with the hole [3, 7]^2, plus L(s): the outside becomes
    // [0, 10 + 2s] x [0, 10 + s]  union  [0, 10 + s] x [0, 10 + 2s], and a
    // point x stays in the hole while the box of x - L(s), [x - 2s, x]^2,
    // fits inside it: the hole [3 + 2s, 7]^2 while 2s < 4.
    let frame = region(vec![Polygon {
        outer: ring(&[(0.0, 0.0), (10.0, 0.0), (10.0, 10.0), (0.0, 10.0)]),
        holes: vec![ring(&[(3.0, 3.0), (7.0, 3.0), (7.0, 7.0), (3.0, 7.0)])],
    }]);
    let outer = |s: f64| 2.0 * (10.0 + 2.0 * s) * (10.0 + s) - (10.0 + s) * (10.0 + s);
    for (s, hole) in [(0.5, 3.0), (1.0, 2.0), (1.5, 1.0), (2.5, 0.0), (4.0, 0.0)] {
        let sum = frame.minkowski_sum(&l_shape(s), tol()).unwrap();
        assert_eq!(sum.component_count(), 1);
        let holes = sum.polygons()[0].holes.len();
        assert_eq!(holes, usize::from(hole > 0.0), "s = {s}");
        assert_eq!(sum.area(), outer(s) - hole * hole, "s = {s}");
        if hole > 0.0 {
            let corners = vertices(&single(&sum.polygons()[0].holes[0]));
            let (a, b) = (3.0 + 2.0 * s, 7.0);
            assert_eq!(corners, vec![(a, a), (a, b), (b, a), (b, b)]);
        }
    }
}

#[test]
fn non_convex_sums_and_erosions_match_the_definition() {
    let r = u_with_hole();
    let (lo, hi) = (Point2::new(-2.5, -2.5), Point2::new(8.5, 7.5));
    for angle in [0.0, 0.7, 2.3] {
        // About the origin, and off it (an erosion must not assume the
        // origin lies in the polygon).
        for centre in [(0.0, 0.0), (0.7, -0.5)] {
            let k = dart(angle, centre);
            let k = Ring {
                points: k
                    .points
                    .iter()
                    .map(|p| Point2::new(p.x * 0.6, p.y * 0.6))
                    .collect(),
            };
            let sum = r.minkowski_sum(&k, tol()).unwrap();
            let eroded = r.minkowski_erosion(&k, tol()).unwrap();
            assert!(!eroded.is_empty());
            check_by_definition(&r, &single(&k), &sum, &eroded, lo, hi);
        }
    }
    let k = l_shape(0.4);
    check_by_definition(
        &r,
        &single(&k),
        &r.minkowski_sum(&k, tol()).unwrap(),
        &r.minkowski_erosion(&k, tol()).unwrap(),
        lo,
        hi,
    );
}

#[test]
fn a_region_with_a_hole_as_the_structuring_shape() {
    // Two components, one a square ring: the region's own hole shows
    // through where neither reaches.
    let k = region(vec![
        Polygon {
            outer: ring(&[(0.0, 0.0), (1.2, 0.0), (1.2, 1.2), (0.0, 1.2)]),
            holes: vec![ring(&[(0.3, 0.3), (0.9, 0.3), (0.9, 0.9), (0.3, 0.9)])],
        },
        Polygon {
            outer: ring(&[(1.5, -0.5), (2.0, -0.5), (1.75, 0.2)]),
            holes: vec![],
        },
    ]);
    let r = u_with_hole();
    let sum = r.minkowski_sum_region(&k, tol()).unwrap();
    let eroded = r.minkowski_erosion_region(&k, tol()).unwrap();
    assert!(!eroded.is_empty());
    check_by_definition(
        &r,
        &k,
        &sum,
        &eroded,
        Point2::new(-3.0, -3.0),
        Point2::new(9.0, 8.0),
    );
    // Erosion by nothing would be the whole plane.
    assert_eq!(
        r.minkowski_erosion_region(&Region::empty(), tol()),
        Err(MinkowskiError::EmptyStructuring)
    );
    assert!(r
        .minkowski_sum_region(&Region::empty(), tol())
        .unwrap()
        .is_empty());
}

#[test]
fn swapping_the_operands_gives_the_same_sum() {
    // Each call cuts the other operand into convex pieces.
    let l = l_shape(1.0);
    let d = dart(0.4, (0.3, 0.2));
    let ld = single(&l).minkowski_sum(&d, tol()).unwrap();
    let dl = single(&d).minkowski_sum(&l, tol()).unwrap();
    assert!(apart(&ld, &dl) < 1e-12, "{}", apart(&ld, &dl));
    assert!((ld.area() - dl.area()).abs() < 1e-12);
    // With holes and several components, through the region form.
    let r = u_with_hole();
    let k = single(&l_shape(0.3))
        .union(&single(&dart(0.4, (2.0, -1.0))), tol())
        .unwrap();
    let rk = r.minkowski_sum_region(&k, tol()).unwrap();
    let kr = k.minkowski_sum_region(&r, tol()).unwrap();
    assert_eq!(rk, kr);
    // The region form cuts the smaller operand, the ring form the ring.
    let ring_form = single(&l_shape(0.3).clone())
        .minkowski_sum_region(&r, tol())
        .unwrap();
    let cut_region = r.minkowski_sum(&l_shape(0.3), tol()).unwrap();
    assert!(apart(&ring_form, &cut_region) < 1e-12);
}

#[test]
fn translating_an_operand_translates_the_sum() {
    let r = u_with_hole();
    let k = l_shape(0.5);
    let t = Vec2::new(3.0, -2.0);
    let sum = r.minkowski_sum(&k, tol()).unwrap();
    // Coordinates are dyadic, so every vertex is exact and the corners
    // match bit for bit (straight-edge vertices may differ).
    let same = |a: &Region, b: &Region| {
        assert_eq!(vertices(a), vertices(b));
        assert_eq!(a.area(), b.area());
        assert_eq!(a.component_count(), b.component_count());
    };
    let moved = r.translate(t).unwrap().minkowski_sum(&k, tol()).unwrap();
    same(&moved, &sum.translate(t).unwrap());
    let k_moved = Ring {
        points: k.points.iter().map(|p| *p + t).collect(),
    };
    same(
        &r.minkowski_sum(&k_moved, tol()).unwrap(),
        &sum.translate(t).unwrap(),
    );
    let eroded = r.minkowski_erosion(&k, tol()).unwrap();
    assert!(!eroded.is_empty());
    same(
        &r.minkowski_erosion(&k_moved, tol()).unwrap(),
        &eroded.translate(-t).unwrap(),
    );
    same(
        &r.translate(t)
            .unwrap()
            .minkowski_erosion(&k, tol())
            .unwrap(),
        &eroded.translate(t).unwrap(),
    );
}

#[test]
fn eroding_by_a_polygon_off_the_origin() {
    // [0, 1] x [0, 10] eroded by [5, 5.5] x [0, 1]: x + K lies inside
    // for x in [-5, -4.5] x [0, 9], away from the region itself. Before
    // 0.3.10 the erosion was also cut to the region, and came out empty.
    let strip = single(&ring(&[(0.0, 0.0), (1.0, 0.0), (1.0, 10.0), (0.0, 10.0)]));
    let k = ring(&[(5.0, 0.0), (5.5, 0.0), (5.5, 1.0), (5.0, 1.0)]);
    let eroded = strip.minkowski_erosion(&k, tol()).unwrap();
    assert_eq!(
        vertices(&eroded),
        vec![(-5.0, 0.0), (-5.0, 9.0), (-4.5, 0.0), (-4.5, 9.0)]
    );
}

/// FNV-1a over the bits of every vertex.
fn digest(regions: &[Region]) -> u64 {
    let mut h = 0xcbf2_9ce4_8422_2325_u64;
    for region in regions {
        for p in region
            .boundary_rings()
            .iter()
            .flat_map(|r| r.points.clone())
        {
            for b in
                p.x.to_bits()
                    .to_le_bytes()
                    .into_iter()
                    .chain(p.y.to_bits().to_le_bytes())
            {
                h = (h ^ u64::from(b)).wrapping_mul(0x0100_0000_01b3);
            }
        }
        h = (h ^ 0xff).wrapping_mul(0x0100_0000_01b3);
    }
    h
}

#[test]
fn the_convex_path_is_unchanged_bit_for_bit() {
    // A digest of results as axiolid-overlay 0.3.9 gave them, before
    // non-convex polygons were accepted: convex polygons holding the
    // origin, and the disc morphology built on them.
    let r = u_with_hole();
    let mut out = Vec::new();
    for angle in [0.0, 0.3, 1.1] {
        let k = Ring {
            points: rectangle(0.8, 0.3, angle),
        };
        out.push(r.minkowski_sum(&k, tol()).unwrap());
        out.push(r.minkowski_erosion(&k, tol()).unwrap());
    }
    let off = ring(&[(0.1, 0.1), (0.9, 0.2), (0.3, 0.7)]);
    out.push(r.minkowski_sum(&off, tol()).unwrap());
    for radius in [0.2, 0.45] {
        out.push(r.dilate_inner(radius, tol()).unwrap());
        out.push(r.dilate_outer(radius, tol()).unwrap());
        out.push(r.erode_inner(radius, tol()).unwrap());
        out.push(r.erode_outer(radius, tol()).unwrap());
    }
    let about = ring(&[(-0.1, -0.1), (0.9, 0.2), (0.3, 0.7)]);
    out.push(r.minkowski_erosion(&about, tol()).unwrap());
    assert_eq!(digest(&out), 13_114_369_021_911_355_074);
}
