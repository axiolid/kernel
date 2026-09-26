//! Where two B-spline surfaces meet (#119, ADR 0077): pair sections.
//!
//! Both sheets are graphs over the plane (their control grids are uniform
//! in x and y, so x and y are linear in the parameters): along any line of
//! the first sheet, the traced curves cross it exactly where the height
//! difference of the two sheets changes sign, which a dense scan counts
//! independently. Every node and every point between lies on both.

use axiolid_core::Point3;
use axiolid_curve::{BSplineSurface, Curve3, KnotSpec};
use axiolid_evaluate::evaluate3;
use axiolid_nurbs::exact_surface_intersection;
use axiolid_surface::Surface;

/// A bi-quadratic graph sheet over `[x0, x0 + 2] x [y0, y0 + 2]`.
fn sheet(x0: f64, y0: f64, h: [[f64; 3]; 3]) -> BSplineSurface {
    BSplineSurface {
        u_degree: 2,
        v_degree: 2,
        control_points: (0..3)
            .map(|i| {
                (0..3)
                    .map(|j| Point3::new(x0 + i as f64, y0 + j as f64, h[i][j]))
                    .collect()
            })
            .collect(),
        u_knots: vec![0.0, 1.0],
        u_multiplicities: vec![3, 3],
        v_knots: vec![0.0, 1.0],
        v_multiplicities: vec![3, 3],
        weights: None,
        u_closed: false,
        v_closed: false,
        knot_spec: KnotSpec::PiecewiseBezier,
        self_intersect: None,
    }
}

/// A sheet's height over `(x, y)`, using its linear x and y.
fn height(b: &BSplineSurface, x: f64, y: f64) -> f64 {
    let x0 = b.control_points[0][0].x;
    let y0 = b.control_points[0][0].y;
    b.jet((x - x0) / 2.0, (y - y0) / 2.0).unwrap().point.z
}

#[test]
fn two_spline_sheets_cross_in_traced_pair_sections() {
    let a = sheet(
        0.0,
        0.0,
        [[0.0, 0.4, 0.1], [0.3, 1.2, 0.2], [0.1, 0.5, 0.0]],
    );
    let b = sheet(
        0.3,
        0.2,
        [[0.6, 0.3, 0.5], [0.2, 0.1, 0.4], [0.7, 0.6, 0.9]],
    );
    let curve =
        exact_surface_intersection(&Surface::BSpline(a.clone()), &Surface::BSpline(b.clone()))
            .expect("a pair section");
    assert!(!curve.branches.is_empty());
    let mut polylines = Vec::new();
    for branch in &curve.branches {
        let Curve3::PairSection(ps) = branch else {
            panic!("expected a pair section");
        };
        let n = 200 * ps.nodes.len();
        let mut line = Vec::new();
        for i in 0..=n {
            let t = ps.end() * i as f64 / n as f64;
            let p = evaluate3(branch, t).unwrap();
            let (sa, sb, _) = ps.solve(t).unwrap();
            assert!((a.jet(sa.x, sa.y).unwrap().point - p).length() < 1e-9);
            assert!((b.jet(sb.x, sb.y).unwrap().point - p).length() < 1e-9);
            line.push(p);
        }
        polylines.push(line);
    }
    // Completeness along lines x = const over the overlap of the sheets.
    let mut bad = Vec::new();
    for k in 0..40 {
        let x = 0.3 + 1.7 * (k as f64 + 0.5) / 40.0;
        let mut scan = 0;
        let steps = 4000;
        let at = |i: usize| 0.2 + 1.8 * i as f64 / steps as f64;
        let mut last = (height(&a, x, at(0)) - height(&b, x, at(0))).signum();
        for i in 1..=steps {
            let now = (height(&a, x, at(i)) - height(&b, x, at(i))).signum();
            if now != last {
                scan += 1;
            }
            last = now;
        }
        let traced: usize = polylines
            .iter()
            .map(|pl| {
                pl.windows(2)
                    .filter(|w| (w[0].x - x) * (w[1].x - x) < 0.0)
                    .count()
            })
            .sum();
        if scan != traced {
            bad.push((x, scan, traced));
        }
    }
    assert!(bad.len() <= 1, "{bad:?}");
}

#[test]
fn a_bowl_and_a_flat_sheet_meet_in_one_closed_loop() {
    // z = x^2 + y^2 over [-1, 1]^2 (the square's Bernstein coefficients are
    // 1, -1, 1), against the flat sheet z = 1/2: a circle of radius 1/sqrt 2
    // lying inside one patch of each, which the search must still find.
    let c = [1.0, -1.0, 1.0];
    let mut bowl = sheet(-1.0, -1.0, [[0.0; 3]; 3]);
    for i in 0..3 {
        for j in 0..3 {
            let p = &mut bowl.control_points[i][j];
            *p = Point3::new(-1.0 + i as f64, -1.0 + j as f64, c[i] + c[j]);
        }
    }
    let flat = sheet(-1.2, -1.2, [[0.5; 3]; 3]);
    let flat = BSplineSurface {
        control_points: flat
            .control_points
            .iter()
            .map(|row| {
                row.iter()
                    .map(|p| Point3::new(-1.2 + 1.2 * (p.x + 1.2), -1.2 + 1.2 * (p.y + 1.2), 0.5))
                    .collect()
            })
            .collect(),
        ..flat
    };
    let curve = exact_surface_intersection(
        &Surface::BSpline(bowl.clone()),
        &Surface::BSpline(flat.clone()),
    )
    .expect("a pair section");
    assert_eq!(curve.branches.len(), 1, "{:?}", curve.branches.len());
    let Curve3::PairSection(ps) = &curve.branches[0] else {
        panic!("expected a pair section");
    };
    assert_eq!(
        ps.nodes.first().unwrap().point,
        ps.nodes.last().unwrap().point,
        "closed"
    );
    for i in 0..=400 {
        let t = ps.end() * i as f64 / 400.0;
        let p = evaluate3(&curve.branches[0], t).unwrap();
        assert!((p.z - 0.5).abs() < 1e-9);
        assert!(((p.x * p.x + p.y * p.y).sqrt() - 0.5f64.sqrt()).abs() < 1e-9);
    }
}

#[test]
fn a_spline_curve_crosses_a_spline_surface_where_a_scan_says() {
    use axiolid_curve::BSplineCurve3;
    use axiolid_nurbs::{exact_curve_surface_intersection, ExactCurveIntersection};
    let c = [1.0, -1.0, 1.0];
    let mut bowl = sheet(-1.0, -1.0, [[0.0; 3]; 3]);
    for i in 0..3 {
        for j in 0..3 {
            bowl.control_points[i][j] = Point3::new(-1.0 + i as f64, -1.0 + j as f64, c[i] + c[j]);
        }
    }
    let curve = BSplineCurve3 {
        degree: 2,
        control_points: vec![
            Point3::new(-0.9, -0.2, 1.0),
            Point3::new(0.0, 0.0, -0.5),
            Point3::new(0.9, 0.3, 1.0),
        ],
        knots: vec![0.0, 1.0],
        multiplicities: vec![3, 3],
        weights: None,
        closed: false,
        self_intersect: None,
        knot_spec: KnotSpec::PiecewiseBezier,
    };
    let as_curve = Curve3::BSpline(curve);
    let Ok(ExactCurveIntersection::Points(hits)) =
        exact_curve_surface_intersection(&as_curve, &Surface::BSpline(bowl))
    else {
        panic!("points");
    };
    let gap = |t: f64| {
        let p = evaluate3(&as_curve, t).unwrap();
        p.z - (p.x * p.x + p.y * p.y)
    };
    let changes = (0..4000)
        .filter(|&i| gap(i as f64 / 4000.0).signum() != gap((i + 1) as f64 / 4000.0).signum())
        .count();
    assert_eq!(hits.len(), changes);
    for hit in &hits {
        let t = hit.parameter.approx();
        assert!(gap(t).abs() < 1e-12, "{t}");
        assert!((evaluate3(&as_curve, t).unwrap() - hit.point).length() < 1e-12);
    }
}
