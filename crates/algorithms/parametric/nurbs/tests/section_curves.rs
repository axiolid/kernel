//! Curve/curve intersection for the section families (#119, B8) and
//! extrema between any pieces (#119, B7).
//!
//! Oracles independent of the routines: every hit is each curve's own
//! point at its own parameter; hit counts equal the sign changes of one
//! curve's implicit residual along the other, scanned densely; distances
//! match closed forms or bracket a dense sample's minimum.

use axiolid_core::{Frame2, Frame3, Interval, Point2, Point3, Tolerance, Vec2, Vec3};
use axiolid_curve::{Circle2, Circle3, Curve2, Curve3, Ellipse2, Line2, Line3, Sinusoid2};
use axiolid_evaluate::{evaluate2, evaluate3};
use axiolid_nurbs::extrema::{minimum_distance, Piece};
use axiolid_nurbs::{
    exact_surface_intersection, section_curve_curve_intersection2,
    section_curve_curve_intersection3,
};
use axiolid_surface::{Cylinder, Sphere, Surface, Torus};

const PI: f64 = std::f64::consts::PI;
const TAU: f64 = std::f64::consts::TAU;

fn tol() -> Tolerance {
    Tolerance::METRE
}

fn frame(origin: Point3, axis: Vec3) -> Frame3 {
    let z = axis.normalize();
    let helper = if z.x.abs() < 0.9 { Vec3::X } else { Vec3::Y };
    let x = helper.cross(z).normalize();
    let y = z.cross(x);
    Frame3 { origin, x, y, z }
}

fn agree2(a: &Curve2, b: &Curve2, hits: &[axiolid_nurbs::CurveCurveHit]) {
    for h in hits {
        let (p, q) = (
            evaluate2(a, h.first).unwrap(),
            evaluate2(b, h.second).unwrap(),
        );
        assert!((p - q).length() < 1e-9, "{p:?} vs {q:?}");
        assert!((Point2::new(h.point.x, h.point.y) - p).length() < 1e-9);
    }
}

/// Sign changes of `residual` along `curve` over `span`.
fn scan2(curve: &Curve2, span: Interval, residual: impl Fn(Point2) -> f64) -> usize {
    let n = 20_000;
    let mut last = residual(evaluate2(curve, span.start).unwrap()).signum();
    let mut count = 0;
    for i in 1..=n {
        let t = span.start + (span.end - span.start) * i as f64 / n as f64;
        let now = residual(evaluate2(curve, t).unwrap()).signum();
        if now != last {
            count += 1;
        }
        last = now;
    }
    count
}

#[test]
fn a_sinusoid_against_a_line_and_an_ellipse() {
    let wave = Curve2::Sinusoid(Sinusoid2 {
        mean: 0.5,
        cosine: 0.3,
        sine: 0.2,
    });
    let span = Interval::new(0.0, TAU);
    let line = Curve2::Line(Line2 {
        origin: Point2::new(0.0, 0.6),
        direction: Vec2::new(1.0, 0.01),
    });
    let hits =
        section_curve_curve_intersection2(&wave, span, &line, Interval::new(-1.0, 8.0), tol())
            .unwrap();
    agree2(&wave, &line, &hits);
    let expect = scan2(&wave, span, |p| p.y - 0.6 - 0.01 * p.x);
    assert_eq!(hits.len(), expect);
    let ellipse = Curve2::Ellipse(Ellipse2 {
        frame: Frame2 {
            origin: Point2::new(2.0, 0.5),
            x: Vec2::X,
            y: Vec2::Y,
        },
        semi_axis_x: 1.5,
        semi_axis_y: 0.25,
    });
    let hits =
        section_curve_curve_intersection2(&wave, span, &ellipse, Interval::new(0.0, TAU), tol())
            .unwrap();
    agree2(&wave, &ellipse, &hits);
    let expect = scan2(&wave, span, |p| {
        ((p.x - 2.0) / 1.5).powi(2) + ((p.y - 0.5) / 0.25).powi(2) - 1.0
    });
    assert_eq!(hits.len(), expect);
    assert!(expect >= 2);
}

#[test]
fn a_traced_pcurve_against_a_circle() {
    // The pcurve of a torus/cylinder section on the torus, against a
    // circle in the torus's parameters.
    let ring = Surface::Torus(Torus {
        frame: frame(Point3::ZERO, Vec3::Z),
        major_radius: 3.0,
        minor_radius: 1.0,
    });
    let pipe = Surface::Cylinder(Cylinder {
        frame: frame(Point3::new(0.0, 3.0, 0.1), Vec3::X),
        radius: 0.5,
    });
    let curve = exact_surface_intersection(&ring, &pipe).unwrap();
    let Curve3::ImplicitSection(s) = &curve.branches[0] else {
        panic!()
    };
    let pcurve = Curve2::Implicit(s.curve.clone());
    let span = Interval::new(0.0, s.curve.end());
    let centre = evaluate2(&pcurve, 0.3 * s.curve.end()).unwrap();
    let circle = Curve2::Circle(Circle2 {
        frame: Frame2 {
            origin: centre + Vec2::new(0.05, 0.02),
            x: Vec2::X,
            y: Vec2::Y,
        },
        radius: 0.2,
    });
    let hits =
        section_curve_curve_intersection2(&pcurve, span, &circle, Interval::new(0.0, TAU), tol())
            .unwrap();
    agree2(&pcurve, &circle, &hits);
    let c = centre + Vec2::new(0.05, 0.02);
    let expect = scan2(&pcurve, span, |p| (p - c).length() - 0.2);
    assert_eq!(hits.len(), expect);
    assert!(expect >= 2);
}

#[test]
fn a_traced_loop_meets_a_circle_on_its_own_cylinder() {
    // The circle lies on the pipe; it meets the torus/pipe loop exactly
    // where it crosses the torus.
    let ring = Surface::Torus(Torus {
        frame: frame(Point3::ZERO, Vec3::Z),
        major_radius: 3.0,
        minor_radius: 1.0,
    });
    let pipe = Surface::Cylinder(Cylinder {
        frame: frame(Point3::new(0.0, 3.0, 0.1), Vec3::X),
        radius: 0.5,
    });
    let curve = exact_surface_intersection(&ring, &pipe).unwrap();
    let off_torus = |p: Point3| (p.x.hypot(p.y) - 3.0).hypot(p.z) - 1.0;
    let mut total = 0;
    for x in [-2.5, -2.4, 2.45] {
        let circle = Curve3::Circle(Circle3 {
            frame: Frame3 {
                origin: Point3::new(x, 3.0, 0.1),
                x: Vec3::Y,
                y: Vec3::Z,
                z: Vec3::X,
            },
            radius: 0.5,
        });
        let mut expect = 0;
        let mut last = off_torus(evaluate3(&circle, 0.0).unwrap()).signum();
        for k in 1..=4000 {
            let now = off_torus(evaluate3(&circle, TAU * k as f64 / 4000.0).unwrap()).signum();
            if now != last {
                expect += 1;
            }
            last = now;
        }
        let mut found = 0;
        for branch in &curve.branches {
            let Curve3::ImplicitSection(s) = branch else {
                panic!()
            };
            let hits = section_curve_curve_intersection3(
                &circle,
                Interval::new(0.0, TAU),
                branch,
                Interval::new(0.0, s.curve.end()),
                tol(),
            )
            .unwrap();
            for h in &hits {
                let (p, q) = (
                    evaluate3(&circle, h.first).unwrap(),
                    evaluate3(branch, h.second).unwrap(),
                );
                assert!((p - q).length() < 1e-9);
            }
            found += hits.len();
        }
        assert_eq!(found, expect, "circle at x = {x}");
        total += found;
    }
    assert!(total >= 2);
}

#[test]
fn distances_with_closed_forms() {
    // Point to sphere.
    let ball = Surface::Sphere(Sphere {
        frame: frame(Point3::new(1.0, 2.0, 0.5), Vec3::Z),
        radius: 0.7,
    });
    let p = Point3::new(3.0, -1.0, 2.0);
    let e = minimum_distance(
        &Piece::Point(p),
        &Piece::Patch {
            surface: &ball,
            lo: Point2::new(-PI, -0.5 * PI),
            hi: Point2::new(PI, 0.5 * PI),
        },
        1e-9,
    )
    .unwrap();
    let want = (p - Point3::new(1.0, 2.0, 0.5)).length() - 0.7;
    assert!(
        e.lower <= want + 1e-12 && e.upper >= want - 1e-12 && e.upper - e.lower <= 1e-9,
        "{e:?} vs {want}"
    );
    // A segment beside the axis of a circle: nearest from its lower end to
    // the circle's nearest point.
    let circle = Curve3::Circle(Circle3 {
        frame: frame(Point3::ZERO, Vec3::Z),
        radius: 1.0,
    });
    let segment = Curve3::Line(Line3 {
        origin: Point3::new(0.3, 0.0, 1.0),
        direction: Vec3::Z,
    });
    let e = minimum_distance(
        &Piece::Curve {
            curve: &segment,
            span: Interval::new(0.0, 1.0),
        },
        &Piece::Curve {
            curve: &circle,
            span: Interval::new(0.0, TAU),
        },
        1e-9,
    )
    .unwrap();
    let want = (0.49f64 + 1.0).sqrt();
    assert!(
        (e.upper - want).abs() <= 1e-9 && e.lower <= want + 1e-12,
        "{e:?}"
    );
    // A ball outside a torus, off its axis: the ball's centre is 2.3 from
    // the tube circle, so the surfaces are 2.3 - 1 - 0.5 apart.
    let ring = Surface::Torus(Torus {
        frame: frame(Point3::ZERO, Vec3::Z),
        major_radius: 3.0,
        minor_radius: 1.0,
    });
    let small = Surface::Sphere(Sphere {
        frame: frame(Point3::new(5.3, 0.0, 0.0), Vec3::Z),
        radius: 0.5,
    });
    let e = minimum_distance(
        &Piece::Patch {
            surface: &ring,
            lo: Point2::new(-PI, -PI),
            hi: Point2::new(PI, PI),
        },
        &Piece::Patch {
            surface: &small,
            lo: Point2::new(-PI, -0.5 * PI),
            hi: Point2::new(PI, 0.5 * PI),
        },
        1e-7,
    )
    .unwrap();
    assert!(
        (e.upper - 0.8).abs() <= 1e-7 && e.lower <= 0.8 + 1e-12,
        "{e:?}"
    );
}

#[test]
fn distance_to_a_traced_section_brackets_sampling() {
    let ring = Surface::Torus(Torus {
        frame: frame(Point3::ZERO, Vec3::Z),
        major_radius: 3.0,
        minor_radius: 1.0,
    });
    let pipe = Surface::Cylinder(Cylinder {
        frame: frame(Point3::new(0.0, 3.0, 0.1), Vec3::X),
        radius: 0.5,
    });
    let curve = exact_surface_intersection(&ring, &pipe).unwrap();
    let branch = &curve.branches[0];
    let Curve3::ImplicitSection(s) = branch else {
        panic!()
    };
    let end = s.curve.end();
    let p = Point3::new(2.0, 4.0, 1.0);
    let e = minimum_distance(
        &Piece::Point(p),
        &Piece::Curve {
            curve: branch,
            span: Interval::new(0.0, end),
        },
        1e-8,
    )
    .unwrap();
    let sampled = (0..=20_000)
        .map(|i| (evaluate3(branch, end * i as f64 / 20_000.0).unwrap() - p).length())
        .fold(f64::INFINITY, f64::min);
    assert!(e.lower <= sampled + 1e-12, "{e:?} vs {sampled}");
    assert!(e.upper <= sampled + 1e-6, "{e:?} vs {sampled}");
    assert!(e.upper - e.lower <= 1e-8);
}

#[test]
fn distance_to_a_spline_patch_brackets_sampling() {
    use axiolid_curve::{BSplineSurface, KnotSpec};
    let heights = [[0.0, 0.4, 0.1], [0.3, 1.0, 0.2], [0.1, 0.5, 0.0]];
    let patch = Surface::BSpline(BSplineSurface {
        u_degree: 2,
        v_degree: 2,
        control_points: (0..3)
            .map(|i| {
                (0..3)
                    .map(|j| Point3::new(i as f64, j as f64, heights[i][j]))
                    .collect()
            })
            .collect(),
        u_knots: vec![0.0, 1.0],
        u_multiplicities: vec![3, 3],
        v_knots: vec![0.0, 1.0],
        v_multiplicities: vec![3, 3],
        weights: Some(vec![
            vec![1.0, 1.2, 1.0],
            vec![0.9, 1.1, 1.0],
            vec![1.0, 1.3, 1.0],
        ]),
        u_closed: false,
        v_closed: false,
        knot_spec: KnotSpec::PiecewiseBezier,
        self_intersect: None,
    });
    let p = Point3::new(0.8, 1.1, 1.6);
    let e = minimum_distance(
        &Piece::Point(p),
        &Piece::Patch {
            surface: &patch,
            lo: Point2::new(0.0, 0.0),
            hi: Point2::new(1.0, 1.0),
        },
        1e-5,
    )
    .unwrap();
    let mut sampled = f64::INFINITY;
    for i in 0..=400 {
        for j in 0..=400 {
            let q = axiolid_evaluate::surface::evaluate(&patch, i as f64 / 400.0, j as f64 / 400.0)
                .unwrap();
            sampled = sampled.min((q - p).length());
        }
    }
    assert!(e.lower <= sampled + 1e-12, "{e:?} vs {sampled}");
    assert!(e.upper <= sampled + 1e-5, "{e:?} vs {sampled}");
    assert!(e.upper - e.lower <= 1e-5);
}
