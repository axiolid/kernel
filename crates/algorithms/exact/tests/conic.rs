//! Exact conic intersection, checked against independent oracles:
//! hand-built configurations with exactly known coordinates, the existing
//! line/circle construction, and a dense numeric sign-change count for
//! random ellipse pairs.

use axiolid_core::Point2;
use axiolid_exact::{
    conic_intersections, line_circle_hits, line_conic_hits, Arith, Circle, Conic,
    ConicIntersection, ConicLineHits, ConicPoint, Dyadic, HitCount, Line,
};
use axiolid_guarantees::Sign;

fn p(x: f64, y: f64) -> Point2 {
    Point2::new(x, y)
}

fn circle(cx: f64, cy: f64, r: f64) -> Conic {
    Conic::circle(p(cx, cy), r).unwrap()
}

fn points(a: &Conic, b: &Conic) -> Vec<ConicPoint> {
    match conic_intersections(a, b).unwrap() {
        ConicIntersection::Points(points) => points,
        ConicIntersection::Overlapping => panic!("unexpected overlap"),
    }
}

/// Every point lies exactly on both conics.
fn on_both(points: &[ConicPoint], a: &Conic, b: &Conic) {
    for pt in points {
        assert_eq!(pt.sign_of_conic(a), Sign::Zero, "on first");
        assert_eq!(pt.sign_of_conic(b), Sign::Zero, "on second");
    }
}

fn dy(v: f64) -> Dyadic {
    Dyadic::from_f64(v)
}

#[test]
fn two_circles_meet_at_integer_points() {
    let (a, b) = (circle(0.0, 0.0, 5.0), circle(6.0, 0.0, 5.0));
    let pts = points(&a, &b);
    assert_eq!(pts.len(), 2);
    on_both(&pts, &a, &b);
    for q in &pts {
        assert_eq!(q.x().cmp_dyadic(&dy(3.0)), Sign::Zero, "x = 3 exactly");
        let y = q.y();
        let four = y.cmp_dyadic(&dy(4.0)) == Sign::Zero || y.cmp_dyadic(&dy(-4.0)) == Sign::Zero;
        assert!(four, "y = +-4 exactly");
        assert!(!q.is_tangent());
    }
}

#[test]
fn tangency_is_exact_and_one_ulp_breaks_it() {
    let (a, b) = (circle(0.0, 0.0, 1.0), circle(2.0, 0.0, 1.0));
    let pts = points(&a, &b);
    assert_eq!(pts.len(), 1);
    assert!(pts[0].is_tangent(), "externally tangent circles");
    assert_eq!(pts[0].x().cmp_dyadic(&dy(1.0)), Sign::Zero);
    // One ulp larger: they now cross at two points, neither tangent.
    let bigger = circle(2.0, 0.0, 1.0f64.next_up());
    let crossing = points(&a, &bigger);
    assert_eq!(crossing.len(), 2);
    assert!(crossing.iter().all(|q| !q.is_tangent()));
    on_both(&crossing, &a, &bigger);
    // One ulp smaller: they miss.
    let smaller = circle(2.0, 0.0, 1.0f64.next_down());
    assert!(points(&a, &smaller).is_empty());
}

#[test]
fn ellipse_and_circle_with_irrational_points() {
    // x^2/4 + y^2 = 1 and x^2 + y^2 = 9/4: x^2 = 5/3, y^2 = 7/12.
    let e = Conic::ellipse(p(0.0, 0.0), p(1.0, 0.0), 2.0, 1.0).unwrap();
    let c = circle(0.0, 0.0, 1.5);
    let pts = points(&e, &c);
    assert_eq!(pts.len(), 4);
    on_both(&pts, &e, &c);
    for q in &pts {
        // Compare x^2 with 5/3 without division: |x| vs sqrt(5/3) via
        // refinement brackets. x is a root of 3x^2 - 5 (up to factors).
        let x = q.x().approx();
        let y = q.y().approx();
        assert!((x * x - 5.0 / 3.0).abs() < 1e-12, "x^2 = 5/3, got {x}");
        assert!((y * y - 7.0 / 12.0).abs() < 1e-12, "y^2 = 7/12, got {y}");
        assert!(!q.is_tangent());
    }
    // Tangent pair: circle of radius 1 touches x^2/4 + y^2 = 1 at (0, +-1).
    let unit = circle(0.0, 0.0, 1.0);
    let touch = points(&e, &unit);
    assert_eq!(touch.len(), 2);
    assert!(touch.iter().all(ConicPoint::is_tangent));
    for q in &touch {
        assert_eq!(q.x().cmp_dyadic(&dy(0.0)), Sign::Zero);
    }
}

#[test]
fn a_hyperbola_needs_a_shear() {
    // xy = 1 has A = C = 0: the unsheared y-coefficient vanishes.
    let h = Conic::from_coefficients([0.0, 1.0, 0.0, 0.0, 0.0, -1.0]).unwrap();
    let c = circle(0.0, 0.0, 2f64.sqrt());
    // sqrt(2) as a double is not sqrt(2): the circle is a hair off
    // tangency, so compute with the exact radius-squared form instead.
    let c_exact = Conic::from_coefficients([1.0, 0.0, 1.0, 0.0, 0.0, -2.0]).unwrap();
    let touch = points(&h, &c_exact);
    assert_eq!(touch.len(), 2, "tangent at (1,1) and (-1,-1)");
    assert!(touch.iter().all(ConicPoint::is_tangent));
    on_both(&touch, &h, &c_exact);
    assert!(touch.iter().all(|q| q.shear() != 0));
    // With the rounded radius the answer is exact too: rounded down, it
    // misses; the circle is just inside the hyperbola's branches.
    let rounded = points(&h, &c);
    let r2 = dy(2f64.sqrt()).square().sub(&dy(2.0)).sign().unwrap();
    match r2 {
        Sign::Negative => assert!(rounded.is_empty()),
        Sign::Positive => assert_eq!(rounded.len(), 4),
        _ => unreachable!(),
    }
}

#[test]
fn line_hits_are_ordered_when_the_leading_coefficient_is_negative() {
    // x^2 - y^2 = 1 along the vertical line x = 2, from (2,-10) to (2,10):
    // alpha = dx^2 - dy^2 = -400 < 0, so the root formula's "minus" branch
    // is the LATER hit and the order must be swapped. Hits at y = -sqrt(3)
    // (t = (10 - sqrt 3) / 20, about 0.413) and y = +sqrt(3) (about 0.587).
    let h = Conic::from_coefficients([1.0, 0.0, -1.0, 0.0, 0.0, -1.0]).unwrap();
    let line = Line::new(p(2.0, -10.0), p(2.0, 10.0)).unwrap();
    let ConicLineHits::Secant(first, second) = line_conic_hits(line, &h).unwrap() else {
        panic!("expected two hits");
    };
    assert_eq!(
        first.cmp_param(0.5).unwrap(),
        Sign::Negative,
        "first before t = 0.5"
    );
    assert_eq!(
        second.cmp_param(0.5).unwrap(),
        Sign::Positive,
        "second after t = 0.5"
    );
    assert_eq!(first.compare_along(&second).unwrap(), Sign::Negative);
    assert!(first.approx_point().y < 0.0 && second.approx_point().y > 0.0);
}

#[test]
fn overlapping_and_degenerate_input() {
    let a = circle(1.0, 2.0, 3.0);
    let doubled = Conic::from_coefficients({
        let c = a.coefficients();
        let mut out = [0.0; 6];
        for (o, v) in out.iter_mut().zip(c) {
            *o = v.to_f64() * 2.0;
        }
        out
    })
    .unwrap();
    assert_eq!(
        conic_intersections(&a, &doubled).unwrap(),
        ConicIntersection::Overlapping
    );
    assert!(Conic::from_coefficients([0.0, 0.0, 0.0, 1.0, 1.0, 1.0]).is_err());
    assert!(Conic::from_coefficients([f64::NAN, 0.0, 1.0, 0.0, 0.0, 0.0]).is_err());
    // Concentric circles never meet.
    assert!(points(&circle(0.0, 0.0, 1.0), &circle(0.0, 0.0, 2.0)).is_empty());
}

#[test]
fn line_conic_matches_line_circle() {
    let mut s: u64 = 0x5EED_0000_1111_2222;
    let mut next = || {
        s ^= s << 13;
        s ^= s >> 7;
        s ^= s << 17;
        ((s >> 11) as f64 / (1u64 << 53) as f64) * 20.0 - 10.0
    };
    let mut secants = 0;
    for _ in 0..2000 {
        let line = Line::new(p(next(), next()), p(next(), next())).unwrap();
        let (cx, cy, r) = (next(), next(), next().abs() + 0.1);
        let old = line_circle_hits(line, Circle::new(p(cx, cy), r).unwrap()).unwrap();
        let new = line_conic_hits(line, &circle(cx, cy, r)).unwrap();
        match (old, new) {
            (HitCount::Missed, ConicLineHits::None) => {}
            (HitCount::Tangent(_), ConicLineHits::Tangent(_)) => {}
            (HitCount::Secant(o1, o2), ConicLineHits::Secant(n1, n2)) => {
                secants += 1;
                // Same parameters: each new hit compares equal to the old
                // hit's approximate value in the same direction as the old.
                for (o, n) in [(o1, &n1), (o2, &n2)] {
                    let t = o.approx_param();
                    assert_eq!(o.cmp_param(t).unwrap(), n.cmp_param(t).unwrap());
                }
                assert_eq!(n1.compare_along(&n2).unwrap(), Sign::Negative);
            }
            (a, b) => panic!("disagree: {a:?} vs {b:?}"),
        }
    }
    assert!(secants > 500, "only {secants} secants");
}

#[test]
fn line_conic_special_cases() {
    // Parallel to an asymptote of xy = 1: exactly one hit.
    let h = Conic::from_coefficients([0.0, 1.0, 0.0, 0.0, 0.0, -1.0]).unwrap();
    let horizontal = Line::new(p(0.0, 2.0), p(1.0, 2.0)).unwrap();
    match line_conic_hits(horizontal, &h).unwrap() {
        ConicLineHits::Single(hit) => {
            // x = 1/2 at y = 2, and the line runs from x = 0: t = 1/2.
            assert_eq!(hit.cmp_param(0.5).unwrap(), Sign::Zero);
        }
        other => panic!("{other:?}"),
    }
    // A line lying on the degenerate conic xy = 0.
    let axes = Conic::from_coefficients([0.0, 1.0, 0.0, 0.0, 0.0, 0.0]).unwrap();
    let x_axis = Line::new(p(-1.0, 0.0), p(1.0, 0.0)).unwrap();
    assert_eq!(
        line_conic_hits(x_axis, &axes).unwrap(),
        ConicLineHits::OnConic
    );
    // Tangent to a rotated ellipse through its vertex.
    let e = Conic::ellipse(p(1.0, 1.0), p(1.0, 1.0), 2.0, 1.0).unwrap();
    // Minor-axis vertex: centre + b * (-1, 1)/sqrt(2) is irrational, so use
    // the major-axis end on the axis direction (1,1) scaled: the tangent at
    // the end of the major axis is perpendicular to (1,1). The end is at
    // centre + 2*(1,1)/sqrt(2) = (1 + sqrt2, 1 + sqrt2): irrational too.
    // Instead check the unrotated case with an exact vertex.
    let _ = e;
    let flat = Conic::ellipse(p(0.0, 0.0), p(1.0, 0.0), 2.0, 1.0).unwrap();
    let top = Line::new(p(-5.0, 1.0), p(5.0, 1.0)).unwrap();
    match line_conic_hits(top, &flat).unwrap() {
        ConicLineHits::Tangent(hit) => assert_eq!(hit.cmp_param(0.5).unwrap(), Sign::Zero),
        other => panic!("{other:?}"),
    }
}

/// Transversal crossings of conic `b` along ellipse `a`, by sampling:
/// the independent numeric oracle for random ellipse pairs. Plain `f64`
/// on purpose: an oracle sharing the exact code would not be independent.
fn sampled_crossings(centre: Point2, axis: Point2, ra: f64, rb: f64, b: &Conic) -> usize {
    let len = (axis.x * axis.x + axis.y * axis.y).sqrt();
    let (ux, uy) = (axis.x / len, axis.y / len);
    let k: Vec<f64> = b.coefficients().iter().map(Dyadic::to_f64).collect();
    let f = |t: f64| {
        let (c, s) = (t.cos(), t.sin());
        let x = centre.x + ra * c * ux - rb * s * uy;
        let y = centre.y + ra * c * uy + rb * s * ux;
        k[0] * x * x + k[1] * x * y + k[2] * y * y + k[3] * x + k[4] * y + k[5]
    };
    let n = 200_000;
    let mut count = 0;
    let mut prev = f(0.0);
    for i in 1..=n {
        let v = f(std::f64::consts::TAU * i as f64 / n as f64);
        if (prev < 0.0) != (v < 0.0) {
            count += 1;
        }
        prev = v;
    }
    count
}

#[test]
fn random_ellipse_pairs_match_the_sampled_count() {
    let mut s: u64 = 0xE111_95E5_0000_0001;
    let mut next = || {
        s ^= s << 13;
        s ^= s >> 7;
        s ^= s << 17;
        ((s >> 11) as f64 / (1u64 << 53) as f64) * 4.0 - 2.0
    };
    let mut histogram = [0usize; 5];
    for case in 0..60 {
        let (c1, a1) = (p(next(), next()), p(next(), next()));
        let (r1, s1) = (next().abs() + 0.5, next().abs() + 0.5);
        let (c2, a2) = (p(next(), next()), p(next(), next()));
        let (r2, s2) = (next().abs() + 0.5, next().abs() + 0.5);
        let e1 = Conic::ellipse(c1, a1, r1, s1).unwrap();
        let e2 = Conic::ellipse(c2, a2, r2, s2).unwrap();
        let pts = points(&e1, &e2);
        let expected = sampled_crossings(c1, a1, r1, s1, &e2);
        assert_eq!(pts.len(), expected, "case {case}");
        on_both(&pts, &e1, &e2);
        for q in &pts {
            let a = q.approx();
            assert!(
                e1.eval(a).to_f64().abs() < 1e-6,
                "case {case}: approx on e1"
            );
            assert!(
                e2.eval(a).to_f64().abs() < 1e-6,
                "case {case}: approx on e2"
            );
        }
        histogram[pts.len()] += 1;
    }
    // Guard against a vacuous pass: 0, 2 and 4 crossings all occur.
    assert!(
        histogram[0] > 0 && histogram[2] > 0 && histogram[4] > 0,
        "{histogram:?}"
    );
}

#[test]
fn side_of_line_at_conic_points() {
    let (a, b) = (circle(0.0, 0.0, 5.0), circle(6.0, 0.0, 5.0));
    for q in points(&a, &b) {
        // Both points are on x = 3: the vertical line through it is Zero.
        assert_eq!(q.side_of_line(p(3.0, -10.0), p(3.0, 10.0)), Sign::Zero);
        // The x axis separates them.
        let above = q.y().cmp_dyadic(&dy(0.0)) == Sign::Positive;
        let want = if above {
            Sign::Positive
        } else {
            Sign::Negative
        };
        assert_eq!(q.side_of_line(p(-1.0, 0.0), p(1.0, 0.0)), want);
    }
}
