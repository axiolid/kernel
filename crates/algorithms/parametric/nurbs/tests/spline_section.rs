//! Sections of B-spline surfaces by analytic surfaces (#119, ADR 0077),
//! traced on piecewise Bernstein fields over the spline's Bezier patches.
//!
//! Oracles independent of the trace: every point lies on the analytic
//! surface (its own distance function) and is the spline's own point; on
//! every iso-line `u = const` of the spline, the traced curves cross it as
//! often as a dense scan of the analytic surface's distance changes sign.
//! The spline's own evaluator is checked against `axiolid-evaluate`.

use axiolid_core::{Frame3, Point2, Point3, Vec3};
use axiolid_curve::{BSplineSurface, Curve3, KnotSpec};
use axiolid_nurbs::implicit_surface_intersection;
use axiolid_surface::{Cylinder, Plane, Sphere, Surface, Torus};

fn frame(origin: Point3, axis: Vec3) -> Frame3 {
    let z = axis.normalize();
    let helper = if z.x.abs() < 0.9 { Vec3::X } else { Vec3::Y };
    let x = helper.cross(z).normalize();
    let y = z.cross(x);
    Frame3 { origin, x, y, z }
}

/// A rational bi-quadratic by bi-cubic sheet over `[0, 4] x [0, 4]` with two
/// interior knots each way, rising and bulging.
fn sheet() -> BSplineSurface {
    let (rows, cols) = (5, 6);
    let mut control = Vec::new();
    let mut weights = Vec::new();
    for i in 0..rows {
        let mut row = Vec::new();
        let mut wrow = Vec::new();
        for j in 0..cols {
            let (x, y) = (i as f64, 0.8 * j as f64);
            let z = 0.3 * (x - 2.0).powi(2) * 0.25
                + 0.2 * (y - 2.0)
                + 0.15 * ((i * 7 + j * 3) % 5) as f64;
            row.push(Point3::new(x, y, z));
            wrow.push(1.0 + 0.1 * ((i + 2 * j) % 3) as f64);
        }
        control.push(row);
        weights.push(wrow);
    }
    BSplineSurface {
        u_degree: 2,
        v_degree: 3,
        control_points: control,
        u_knots: vec![0.0, 1.5, 2.5, 4.0],
        u_multiplicities: vec![3, 1, 1, 3],
        v_knots: vec![0.0, 1.0, 3.0, 4.0],
        v_multiplicities: vec![4, 1, 1, 4],
        weights: Some(weights),
        u_closed: false,
        v_closed: false,
        knot_spec: KnotSpec::Unspecified,
        self_intersect: None,
    }
}

#[test]
fn the_spline_jet_matches_the_evaluator() {
    let b = sheet();
    let s = Surface::BSpline(b.clone());
    for (u, v) in [(0.3, 0.7), (1.2, 2.2), (2.7, 3.9), (3.97, 0.03), (2.2, 1.4)] {
        let jet = b.jet(u, v).unwrap();
        let p = axiolid_evaluate::surface::evaluate(&s, u, v).unwrap();
        let (su, sv) = axiolid_evaluate::surface::partials(&s, u, v).unwrap();
        assert!((jet.point - p).length() < 1e-12, "{u} {v}");
        assert!(
            (jet.u - su).length() < 1e-10 * (1.0 + su.length()),
            "{u} {v}"
        );
        assert!(
            (jet.v - sv).length() < 1e-10 * (1.0 + sv.length()),
            "{u} {v}"
        );
        // Second partials by central differences of the first.
        let h = 1e-5;
        let fd_uu = (b.jet(u + h, v).unwrap().u - b.jet(u - h, v).unwrap().u) / (2.0 * h);
        let fd_uv = (b.jet(u, v + h).unwrap().u - b.jet(u, v - h).unwrap().u) / (2.0 * h);
        let fd_vv = (b.jet(u, v + h).unwrap().v - b.jet(u, v - h).unwrap().v) / (2.0 * h);
        assert!(
            (jet.uu - fd_uu).length() < 1e-5 * (1.0 + fd_uu.length()),
            "uu at {u} {v}"
        );
        assert!(
            (jet.uv - fd_uv).length() < 1e-5 * (1.0 + fd_uv.length()),
            "uv at {u} {v}"
        );
        assert!(
            (jet.vv - fd_vv).length() < 1e-5 * (1.0 + fd_vv.length()),
            "vv at {u} {v}"
        );
    }
}

fn off(surface: &Surface, p: Point3) -> f64 {
    match surface {
        Surface::Plane(s) => s.frame.z.normalize().dot(p - s.frame.origin),
        Surface::Sphere(s) => (p - s.frame.origin).length() - s.radius,
        Surface::Cylinder(s) => {
            let axis = s.frame.z.normalize();
            let d = p - s.frame.origin;
            (d - axis * d.dot(axis)).length() - s.radius
        }
        Surface::Torus(t) => {
            let axis = t.frame.z.normalize();
            let d = p - t.frame.origin;
            let h = d.dot(axis);
            (((d - axis * h).length() - t.major_radius).hypot(h)) - t.minor_radius
        }
        other => panic!("{other:?}"),
    }
}

fn check(other: Surface) -> usize {
    let b = sheet();
    let spline = Surface::BSpline(b.clone());
    let curves = implicit_surface_intersection(&spline, &other, None).expect("a section");
    for c in &curves {
        let n = 400;
        for i in 0..=n {
            let t = c.curve.end() * i as f64 / n as f64;
            let p = c.point(t).unwrap();
            let uv = c.curve.point(t).unwrap();
            let on = axiolid_evaluate::surface::evaluate(&spline, uv.x, uv.y).unwrap();
            assert!((p - on).length() < 1e-12);
            assert!(off(&other, p).abs() < 1e-9, "{} off at {t}", off(&other, p));
        }
    }
    // Completeness along iso-lines u = const.
    let pcurves: Vec<Vec<Point2>> = curves
        .iter()
        .map(|c| {
            (0..=20_000)
                .map(|i| c.curve.point(c.curve.end() * i as f64 / 20_000.0).unwrap())
                .collect()
        })
        .collect();
    let mut bad = Vec::new();
    for k in 0..60 {
        let u = 4.0 * (k as f64 + 0.5) / 60.0;
        let mut scan = 0;
        let mut last = off(&other, b.jet(u, 0.0).unwrap().point).signum();
        for s in 1..=3000 {
            let now = off(&other, b.jet(u, 4.0 * s as f64 / 3000.0).unwrap().point).signum();
            if now != last {
                scan += 1;
            }
            last = now;
        }
        let traced: usize = pcurves
            .iter()
            .map(|pc| {
                pc.windows(2)
                    .filter(|w| (w[0].x - u) * (w[1].x - u) < 0.0 || w[1].x == u)
                    .count()
            })
            .sum();
        if scan != traced {
            bad.push((u, scan, traced));
        }
    }
    assert!(bad.len() <= 1, "{bad:?}");
    curves.len()
}

#[test]
fn a_spline_sheet_cut_by_a_plane() {
    let plane = Surface::Plane(Plane {
        frame: frame(Point3::new(2.0, 2.0, 0.6), Vec3::new(0.3, -0.2, 1.0)),
    });
    assert!(check(plane) >= 1);
}

#[test]
fn a_spline_sheet_pierced_by_a_pipe() {
    let pipe = Surface::Cylinder(Cylinder {
        frame: frame(Point3::new(2.1, 1.9, 0.0), Vec3::new(0.1, 0.05, 1.0)),
        radius: 0.7,
    });
    assert_eq!(check(pipe), 1);
}

#[test]
fn a_spline_sheet_and_a_ball() {
    let ball = Surface::Sphere(Sphere {
        frame: frame(Point3::new(2.3, 2.0, 0.8), Vec3::Z),
        radius: 1.1,
    });
    assert!(check(ball) >= 1);
}

#[test]
fn a_spline_sheet_through_a_torus() {
    let ring = Surface::Torus(Torus {
        frame: frame(Point3::new(2.0, 2.0, 0.5), Vec3::new(0.2, 0.0, 1.0)),
        major_radius: 1.4,
        minor_radius: 0.5,
    });
    assert!(check(ring) >= 1);
}

#[test]
fn a_traced_spline_section_inverts() {
    let b = sheet();
    let spline = Surface::BSpline(b);
    let pipe = Surface::Cylinder(Cylinder {
        frame: frame(Point3::new(2.1, 1.9, 0.0), Vec3::Z),
        radius: 0.7,
    });
    let curves = implicit_surface_intersection(&spline, &pipe, None).unwrap();
    let curve = Curve3::ImplicitSection(curves[0].clone());
    for t in [0.3, 5.5, 11.2] {
        if t > curves[0].curve.end() {
            continue;
        }
        let p = axiolid_evaluate::evaluate3(&curve, t).unwrap();
        let back =
            axiolid_evaluate::curve::invert3(&curve, p, axiolid_core::Tolerance::METRE).unwrap();
        let q = axiolid_evaluate::evaluate3(&curve, back).unwrap();
        assert!((p - q).length() < 1e-9);
    }
}

#[test]
fn lines_and_circles_against_a_spline_sheet() {
    use axiolid_curve::{Circle3, Line3};
    use axiolid_nurbs::{exact_curve_surface_intersection, ExactCurveIntersection};
    let b = sheet();
    let spline = Surface::BSpline(b.clone());
    let tol = axiolid_core::Tolerance::METRE;
    // A near-vertical line through the sheet meets it once.
    let line = Curve3::Line(Line3 {
        origin: Point3::new(1.7, 2.3, -5.0),
        direction: Vec3::new(0.05, -0.02, 1.0),
    });
    let ExactCurveIntersection::Points(hits) =
        exact_curve_surface_intersection(&line, &spline).unwrap()
    else {
        panic!("contained");
    };
    assert_eq!(hits.len(), 1, "{hits:?}");
    // A horizontal circle crossing the sheet.
    let circle = Curve3::Circle(Circle3 {
        frame: frame(Point3::new(2.0, 2.0, 0.55), Vec3::Z),
        radius: 1.2,
    });
    let ExactCurveIntersection::Points(round) =
        exact_curve_surface_intersection(&circle, &spline).unwrap()
    else {
        panic!("contained");
    };
    assert!(!round.is_empty());
    // Every hit is the curve's own point at its parameter, and a point of
    // the sheet.
    for (curve, hit) in hits
        .iter()
        .map(|h| (&line, h))
        .chain(round.iter().map(|h| (&circle, h)))
    {
        let at = axiolid_evaluate::evaluate3(curve, hit.parameter.approx()).unwrap();
        assert!((at - hit.point).length() < 1e-9);
        let (u, v) = axiolid_evaluate::surface::locate(&spline, hit.point, tol).unwrap();
        assert!((b.jet(u, v).unwrap().point - hit.point).length() < 1e-9);
    }
}
