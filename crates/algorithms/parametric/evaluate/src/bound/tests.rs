//! Every bound here is checked against dense sampling of the exact value
//! it bounds, and against a ceiling so a bound cannot pass by being huge.

use super::*;
use crate::curve::{
    derivative2, derivative3, evaluate2, evaluate3, flatten2, flatten3, second_derivative3,
};
use axiolid_core::{Interval, Point2, Point3};
use axiolid_curve::{BSplineCurve2, BSplineCurve3, Circle2, Ellipse2, KnotSpec};

fn frame2() -> Frame2 {
    Frame2 {
        origin: Point2::new(0.3, -0.2),
        x: Vec2::X,
        y: Vec2::Y,
    }
}

fn frame3() -> Frame3 {
    Frame3 {
        origin: Point3::new(0.1, 0.2, 0.3),
        x: Vec3::X,
        y: Vec3::Y,
        z: Vec3::Z,
    }
}

fn ellipse2() -> Curve2 {
    Curve2::Ellipse(Ellipse2 {
        frame: frame2(),
        semi_axis_x: 0.8,
        semi_axis_y: 0.15,
    })
}

/// A cubic and a rational quadratic, clamped, with one interior double knot.
fn splines2() -> Vec<Curve2> {
    let points = vec![
        Point2::new(0.0, 0.0),
        Point2::new(0.3, 0.5),
        Point2::new(0.7, -0.4),
        Point2::new(1.0, 0.2),
        Point2::new(1.4, 0.6),
        Point2::new(1.9, 0.0),
    ];
    let cubic = BSplineCurve2 {
        degree: 3,
        control_points: points.clone(),
        knots: vec![0.0, 0.4, 0.7, 1.0],
        multiplicities: vec![4, 1, 1, 4],
        weights: None,
        closed: false,
        self_intersect: None,
        knot_spec: KnotSpec::Unspecified,
    };
    let rational = BSplineCurve2 {
        degree: 2,
        control_points: points,
        knots: vec![0.0, 1.0, 2.0, 3.0, 4.0],
        multiplicities: vec![3, 1, 1, 1, 3],
        weights: Some(vec![1.0, 2.0, 0.5, 1.5, 0.7, 1.0]),
        closed: false,
        self_intersect: None,
        knot_spec: KnotSpec::Unspecified,
    };
    vec![Curve2::BSpline(cubic), Curve2::BSpline(rational)]
}

fn spline3() -> Curve3 {
    Curve3::BSpline(BSplineCurve3 {
        degree: 3,
        control_points: vec![
            Point3::new(0.0, 0.0, 0.0),
            Point3::new(0.3, 0.4, 0.1),
            Point3::new(0.8, -0.2, 0.3),
            Point3::new(1.0, 0.5, 0.2),
            Point3::new(1.5, 0.3, -0.2),
        ],
        knots: vec![0.0, 0.5, 1.0],
        multiplicities: vec![4, 1, 4],
        weights: Some(vec![1.0, 1.4, 0.8, 1.2, 1.0]),
        closed: false,
        self_intersect: None,
        knot_spec: KnotSpec::Unspecified,
    })
}

/// Largest distance from dense samples of `c` over `[a, b]` to its chord.
fn measured_chord2(curve: &Curve2, a: Scalar, b: Scalar) -> Scalar {
    let (pa, pb) = (evaluate2(curve, a).unwrap(), evaluate2(curve, b).unwrap());
    (0..=400)
        .map(|i| {
            let t = a + (b - a) * i as Scalar / 400.0;
            segment_distance2(evaluate2(curve, t).unwrap(), pa, pb)
        })
        .fold(0.0, Scalar::max)
}

fn segment_distance2(p: Point2, a: Point2, b: Point2) -> Scalar {
    let ab = b - a;
    let t = if ab.length_squared() > 0.0 {
        ((p - a).dot(ab) / ab.length_squared()).clamp(0.0, 1.0)
    } else {
        0.0
    };
    (p - (a + ab * t)).length()
}

#[test]
fn chord_bounds_cover_the_sampled_deviation_and_stay_close_to_it() {
    let mut curves = vec![
        ellipse2(),
        Curve2::Circle(Circle2 {
            frame: frame2(),
            radius: 0.4,
        }),
    ];
    curves.extend(splines2());
    for curve in &curves {
        let domain = crate::curve::domain2(curve);
        let breaks = continuity_breaks2(curve, 1);
        for n in [3_usize, 8, 32] {
            for i in 0..n {
                let a = domain.start + (domain.end - domain.start) * i as Scalar / n as Scalar;
                let b =
                    domain.start + (domain.end - domain.start) * (i + 1) as Scalar / n as Scalar;
                let Some(bound) = chord_bound2(curve, a, b) else {
                    assert!(
                        breaks_inside(&breaks, a, b),
                        "{curve:?} refused a smooth span"
                    );
                    continue;
                };
                let measured = measured_chord2(curve, a, b);
                assert!(
                    measured <= bound,
                    "{curve:?} [{a}, {b}]: measured {measured} above bound {bound}"
                );
            }
        }
    }
}

#[test]
fn an_ellipse_chord_bound_is_tight_on_short_spans() {
    let curve = ellipse2();
    // A short span at the end of the major axis, where c'' is normal to
    // the curve: the bound and the sagitta agree to second order.
    let (a, b) = (-0.05, 0.05);
    let bound = chord_bound2(&curve, a, b).unwrap();
    let measured = measured_chord2(&curve, a, b);
    assert!(bound / measured < 1.05, "ratio {}", bound / measured);
}

#[test]
fn derivative_bounds_cover_sampled_derivatives() {
    let mut curves = vec![ellipse2()];
    curves.extend(splines2());
    for curve in &curves {
        let domain = crate::curve::domain2(curve);
        let bounds = curve_derivative_bounds2(curve, domain.start, domain.end).unwrap();
        for i in 0..=500 {
            let t = domain.start + (domain.end - domain.start) * i as Scalar / 500.0;
            let d1 = derivative2(curve, t).unwrap().length();
            let d2 = crate::curve::second_derivative2(curve, t).unwrap().length();
            assert!(d1 <= bounds.first, "{curve:?} c' {d1} > {}", bounds.first);
            assert!(
                d2 <= bounds.second,
                "{curve:?} c'' {d2} > {}",
                bounds.second
            );
        }
    }
    let curve = spline3();
    let bounds = curve_derivative_bounds3(&curve, 0.0, 1.0).unwrap();
    for i in 0..=500 {
        let t = i as Scalar / 500.0;
        assert!(derivative3(&curve, t).unwrap().length() <= bounds.first);
        assert!(second_derivative3(&curve, t).unwrap().length() <= bounds.second);
        // The third derivative, by central differences of the exact second.
        let h = 1e-5;
        if t > h && t < 1.0 - h && (t - 0.5).abs() > h {
            let d3 = (second_derivative3(&curve, t + h).unwrap()
                - second_derivative3(&curve, t - h).unwrap())
                / (2.0 * h);
            assert!(
                d3.length() <= bounds.third * (1.0 + 1e-4),
                "c''' {} > {}",
                d3.length(),
                bounds.third
            );
        }
    }
}

#[test]
fn flatten_meets_the_tolerance_where_the_midpoint_test_alone_did_not() {
    // A flat ellipse flattened coarsely: the sagitta test can accept a
    // span whose chord the arc crosses at the midpoint.
    let tolerance = 1e-3;
    let mut curves = vec![ellipse2()];
    curves.extend(splines2());
    for curve in &curves {
        assert!(certifies_flattening2(curve));
        let domain = crate::curve::domain2(curve);
        let points = flatten2(curve, domain, tolerance, 24).unwrap();
        // Every dense sample lies within tolerance of the polyline.
        for i in 0..=4000 {
            let t = domain.start + (domain.end - domain.start) * i as Scalar / 4000.0;
            let p = evaluate2(curve, t).unwrap();
            let d = points
                .windows(2)
                .map(|w| segment_distance2(p, w[0], w[1]))
                .fold(Scalar::INFINITY, Scalar::min);
            assert!(d <= tolerance, "{curve:?}: {d} at t = {t}");
        }
    }
    let curve = spline3();
    let points = flatten3(
        &curve,
        Interval {
            start: 0.0,
            end: 1.0,
        },
        tolerance,
        24,
    )
    .unwrap();
    assert!(points.len() >= 2);
}

#[test]
fn a_corner_knot_is_kept_as_a_vertex() {
    // Degree 2 with a double interior knot: C^0 there, a corner.
    let curve = Curve2::BSpline(BSplineCurve2 {
        degree: 2,
        control_points: vec![
            Point2::new(0.0, 0.0),
            Point2::new(0.5, 1.0),
            Point2::new(1.0, 0.0),
            Point2::new(1.5, 1.0),
            Point2::new(2.0, 0.0),
        ],
        knots: vec![0.0, 0.5, 1.0],
        multiplicities: vec![3, 2, 3],
        weights: None,
        closed: false,
        self_intersect: None,
        knot_spec: KnotSpec::PiecewiseBezier,
    });
    assert_eq!(continuity_breaks2(&curve, 1), vec![0.5]);
    assert!(chord_bound2(&curve, 0.25, 0.75).is_none());
    let points = flatten2(
        &curve,
        Interval {
            start: 0.0,
            end: 1.0,
        },
        1e-3,
        24,
    )
    .unwrap();
    let corner = evaluate2(&curve, 0.5).unwrap();
    assert!(points.iter().any(|p| (*p - corner).length() == 0.0));
}

#[test]
fn a_rational_conic_lies_within_its_bound() {
    // A quarter circle as a rational quadratic.
    let w = core::f64::consts::FRAC_1_SQRT_2;
    let curve = Curve3::BSpline(BSplineCurve3 {
        degree: 2,
        control_points: vec![
            Point3::new(1.0, 0.0, 0.0),
            Point3::new(1.0, 1.0, 0.0),
            Point3::new(0.0, 1.0, 0.0),
        ],
        knots: vec![0.0, 1.0],
        multiplicities: vec![3, 3],
        weights: Some(vec![1.0, w, 1.0]),
        closed: false,
        self_intersect: None,
        knot_spec: KnotSpec::Unspecified,
    });
    for n in [1_usize, 4, 16] {
        for i in 0..n {
            let (a, b) = (i as Scalar / n as Scalar, (i + 1) as Scalar / n as Scalar);
            let bound = chord_bound3(&curve, a, b).unwrap();
            let (pa, pb) = (evaluate3(&curve, a).unwrap(), evaluate3(&curve, b).unwrap());
            let chord_angle = 2.0 * (0.5 * (pb - pa).length()).asin();
            let sagitta = 1.0 - (0.5 * chord_angle).cos();
            assert!(sagitta <= bound, "n {n}: sagitta {sagitta} > {bound}");
            assert!(bound < 20.0 * sagitta, "n {n}: bound {bound} vs {sagitta}");
        }
    }
}

#[test]
fn surface_bounds_cover_sampled_jets() {
    use axiolid_surface::{Cone, Cylinder, Sphere, Torus};
    let f = frame3();
    let surfaces = vec![
        Surface::Cylinder(Cylinder {
            frame: f,
            radius: 0.4,
        }),
        Surface::Cone(Cone {
            frame: f,
            radius: 0.5,
            semi_angle: 0.3,
        }),
        Surface::Sphere(Sphere {
            frame: f,
            radius: 0.7,
        }),
        Surface::Torus(Torus {
            frame: f,
            major_radius: 0.5,
            minor_radius: 0.2,
        }),
        Surface::BSpline(BSplineSurface {
            u_degree: 2,
            v_degree: 3,
            control_points: (0..4)
                .map(|i| {
                    (0..5)
                        .map(|j| {
                            Point3::new(
                                i as Scalar * 0.3,
                                j as Scalar * 0.25,
                                ((i * 7 + j * 3) % 5) as Scalar * 0.1,
                            )
                        })
                        .collect()
                })
                .collect(),
            u_knots: vec![0.0, 0.5, 1.0],
            u_multiplicities: vec![3, 1, 3],
            v_knots: vec![0.0, 0.6, 1.0],
            v_multiplicities: vec![4, 1, 4],
            weights: Some(
                (0..4)
                    .map(|i| {
                        (0..5)
                            .map(|j| 0.7 + 0.1 * ((i + 2 * j) % 4) as Scalar)
                            .collect()
                    })
                    .collect(),
            ),
            u_closed: false,
            v_closed: false,
            knot_spec: KnotSpec::Unspecified,
            self_intersect: None,
        }),
    ];
    for surface in &surfaces {
        let oracle = SurfaceBoundOracle::new(surface).unwrap();
        for (u, v) in [((0.1, 0.9), (0.05, 0.6)), ((0.2, 0.35), (0.3, 0.45))] {
            let bounds = oracle.bounds(u, v).unwrap();
            for i in 0..=30 {
                for j in 0..=30 {
                    let uu = u.0 + (u.1 - u.0) * i as Scalar / 30.0;
                    let vv = v.0 + (v.1 - v.0) * j as Scalar / 30.0;
                    let jet = crate::surface::jet(surface, uu, vv).unwrap();
                    for (value, bound, name) in [
                        (jet.du.length(), bounds.du, "du"),
                        (jet.dv.length(), bounds.dv, "dv"),
                        (jet.duu.length(), bounds.duu, "duu"),
                        (jet.duv.length(), bounds.duv, "duv"),
                        (jet.dvv.length(), bounds.dvv, "dvv"),
                    ] {
                        assert!(
                            value <= bound * (1.0 + 1e-9),
                            "{surface:?} {name} {value} > {bound} at ({uu}, {vv})"
                        );
                    }
                }
            }
        }
    }
}

#[test]
fn a_heavily_weighted_rational_arc_lies_within_its_chord_bound() {
    // The middle weight pulls the arc hard towards its control point; the
    // weight's second derivative is what the projective bound must carry.
    let curve = Curve3::BSpline(BSplineCurve3 {
        degree: 2,
        control_points: vec![
            Point3::new(0.0, 0.0, 0.0),
            Point3::new(0.5, 1.0, 0.0),
            Point3::new(1.0, 0.0, 0.0),
        ],
        knots: vec![0.0, 1.0],
        multiplicities: vec![3, 3],
        weights: Some(vec![1.0, 8.0, 1.0]),
        closed: false,
        self_intersect: None,
        knot_spec: KnotSpec::Unspecified,
    });
    for n in [1_usize, 2, 4, 8, 16, 64] {
        for i in 0..n {
            let (a, b) = (i as Scalar / n as Scalar, (i + 1) as Scalar / n as Scalar);
            let bound = chord_bound3(&curve, a, b).unwrap();
            let (pa, pb) = (evaluate3(&curve, a).unwrap(), evaluate3(&curve, b).unwrap());
            let measured = (0..=400)
                .map(|k| {
                    let p = evaluate3(&curve, a + (b - a) * k as Scalar / 400.0).unwrap();
                    let ab = pb - pa;
                    let t = ((p - pa).dot(ab) / ab.length_squared()).clamp(0.0, 1.0);
                    (p - (pa + ab * t)).length()
                })
                .fold(0.0, Scalar::max);
            assert!(measured <= bound, "n {n} span {i}: {measured} > {bound}");
        }
    }
}

/// Exact distance from `p` to a triangle, by the closest feature.
fn triangle_distance(p: Point3, [a, b, c]: [Point3; 3]) -> Scalar {
    let n = (b - a).cross(c - a);
    let inside = |q: Point3| {
        let s0 = (b - a).cross(q - a).dot(n);
        let s1 = (c - b).cross(q - b).dot(n);
        let s2 = (a - c).cross(q - c).dot(n);
        s0 >= 0.0 && s1 >= 0.0 && s2 >= 0.0
    };
    let segment = |x: Point3, y: Point3| {
        let d = y - x;
        let t = ((p - x).dot(d) / d.length_squared()).clamp(0.0, 1.0);
        (p - (x + d * t)).length()
    };
    let unit = n.normalize();
    let foot = p - unit * (p - a).dot(unit);
    if inside(foot) {
        (p - foot).length()
    } else {
        segment(a, b).min(segment(b, c)).min(segment(c, a))
    }
}

#[test]
fn a_heavily_weighted_rational_surface_lies_within_its_interpolation_bound() {
    let weights: Vec<Vec<Scalar>> = (0..3)
        .map(|i| {
            (0..3)
                .map(|j| if i == 1 && j == 1 { 6.0 } else { 1.0 })
                .collect()
        })
        .collect();
    let surface = Surface::BSpline(BSplineSurface {
        u_degree: 2,
        v_degree: 2,
        control_points: (0..3)
            .map(|i| {
                (0..3)
                    .map(|j| {
                        let z = if i == 1 && j == 1 { 0.8 } else { 0.0 };
                        Point3::new(i as Scalar * 0.5, j as Scalar * 0.5, z)
                    })
                    .collect()
            })
            .collect(),
        u_knots: vec![0.0, 1.0],
        u_multiplicities: vec![3, 3],
        v_knots: vec![0.0, 1.0],
        v_multiplicities: vec![3, 3],
        weights: Some(weights),
        u_closed: false,
        v_closed: false,
        knot_spec: KnotSpec::Unspecified,
        self_intersect: None,
    });
    let oracle = SurfaceBoundOracle::new(&surface).unwrap();
    for n in [2_usize, 4, 8] {
        let h = 1.0 / n as Scalar;
        for i in 0..n {
            for j in 0..n {
                let (u, v) = (i as Scalar * h, j as Scalar * h);
                let b = oracle.bounds((u, u + h), (v, v + h)).unwrap();
                let [x, y, z] = b.interpolation;
                // A right triangle: its enclosing circle is half the
                // hypotenuse, in the metric scaled by sqrt(A + B), sqrt(C + B).
                let r2 = 0.25 * ((x + y) * h * h + (z + y) * h * h);
                let bound = 0.5 * r2;
                let corners = [(u, v), (u + h, v), (u, v + h)]
                    .map(|(s, t)| crate::surface::evaluate(&surface, s, t).unwrap());
                for k in 0..=20 {
                    for l in 0..=(20 - k) {
                        let (s, t) = (u + h * k as Scalar / 20.0, v + h * l as Scalar / 20.0);
                        let p = crate::surface::evaluate(&surface, s, t).unwrap();
                        let d = triangle_distance(p, corners);
                        assert!(d <= bound, "n {n} cell ({i}, {j}): {d} > {bound}");
                    }
                }
            }
        }
    }
}

/// The unit circle `u^2 + v^2 - 1` as an implicit curve: a regular cell
/// along `u` (upper arc, `v` solved in `[0.5, 1.5]`), then a bridge cell
/// (a cubic) leaving it, as a trace into a crossing ends (ADR 0077).
fn implicit_circle() -> Curve2 {
    use axiolid_curve::{Axis, Basis, Field2, ImplicitCell, ImplicitCurve2, SeriesField2};
    let field = Field2::Series(SeriesField2 {
        u: Basis::Power,
        v: Basis::Power,
        coefficients: vec![vec![-1.0, 0.0, 1.0], vec![0.0], vec![1.0]],
    });
    let regular = ImplicitCell {
        axis: Axis::U,
        from: -0.6,
        to: 0.6,
        low: 0.5,
        high: 1.5,
        bridge: None,
    };
    let bridge = ImplicitCell::bridge(
        Point2::new(0.6, 0.8),
        Point2::new(0.9, 0.3),
        Vec2::new(0.8, -0.6),
        Vec2::new(1.0, -1.5),
    );
    Curve2::Implicit(ImplicitCurve2 {
        field,
        cells: vec![regular, bridge],
    })
}

#[test]
fn implicit_chord_bounds_cover_the_sampled_deviation_cell_by_cell() {
    // #249: the boolean's pcurves on a filleted beam are implicit; the
    // deviation path flattens them with these bounds.
    let curve = implicit_circle();
    assert_eq!(continuity_breaks2(&curve, 1), vec![1.0]);
    // Across the cell join: refused, the caller splits there.
    assert_eq!(chord_bound2(&curve, 0.5, 1.5), None);
    for (lo, hi) in [(0.0, 1.0), (1.0, 2.0)] {
        let mut previous = Scalar::INFINITY;
        for n in [1_usize, 4, 16, 64] {
            let mut worst: Scalar = 0.0;
            for i in 0..n {
                let a = lo + (hi - lo) * i as Scalar / n as Scalar;
                let b = lo + (hi - lo) * (i + 1) as Scalar / n as Scalar;
                let bound = chord_bound2(&curve, a, b).expect("bounded");
                let measured = measured_chord2(&curve, a, b);
                assert!(
                    measured <= bound,
                    "[{a}, {b}]: measured {measured} above bound {bound}"
                );
                worst = worst.max(bound);
            }
            // Halving shrinks the bound (second order: by about four).
            assert!(worst < previous, "{lo}: {worst} after {previous}");
            previous = worst;
        }
        assert!(previous < 1e-3, "{lo}: {previous}");
    }
}
