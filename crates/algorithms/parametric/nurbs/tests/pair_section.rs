//! Where two B-spline surfaces meet (#119, ADR 0077): pair sections.
//!
//! Both sheets are graphs over the plane (their control grids are uniform
//! in x and y, so x and y are linear in the parameters): along any line of
//! the first sheet, the traced curves cross it exactly where the height
//! difference of the two sheets changes sign, which a dense scan counts
//! independently. Every node and every point between lies on both.

use axiolid_core::{Point2, Point3};
use axiolid_curve::{BSplineSurface, Curve3, KnotSpec};
use axiolid_evaluate::evaluate3;
use axiolid_nurbs::{
    exact_surface_intersection, spline_pair_intersection, ExactIntersectionRefusal,
};
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
            Point3::new(0.0, 0.0, -1.5),
            Point3::new(0.9, 0.3, 1.0),
        ],
        // Knots over [0, 2], so the curve's parameter is not a Bezier
        // piece's local one.
        knots: vec![0.0, 2.0],
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
        .filter(|&i| gap(i as f64 / 2000.0).signum() != gap((i + 1) as f64 / 2000.0).signum())
        .count();
    // Down through the bowl and back up.
    assert_eq!(changes, 2);
    assert_eq!(hits.len(), changes);
    for hit in &hits {
        let t = hit.parameter.approx();
        assert!(gap(t).abs() < 1e-12, "{t}");
        assert!((evaluate3(&as_curve, t).unwrap() - hit.point).length() < 1e-12);
    }
}

/// The bowl `z = x^2 + y^2` over `[-1, 1]^2`.
fn bowl() -> BSplineSurface {
    let c = [1.0, -1.0, 1.0];
    let mut b = sheet(-1.0, -1.0, [[0.0; 3]; 3]);
    for i in 0..3 {
        for j in 0..3 {
            b.control_points[i][j] = Point3::new(-1.0 + i as f64, -1.0 + j as f64, c[i] + c[j]);
        }
    }
    b
}

/// A flat sheet at height `z` over `[-1.2, 1.2]^2`.
fn flat(z: f64) -> BSplineSurface {
    let mut f = sheet(-1.2, -1.2, [[z; 3]; 3]);
    for i in 0..3 {
        for j in 0..3 {
            f.control_points[i][j] = Point3::new(-1.2 + 1.2 * i as f64, -1.2 + 1.2 * j as f64, z);
        }
    }
    f
}

#[test]
fn a_small_window_inside_one_patch_pair_still_finds_its_stretch() {
    // Two gently tilted sheets over [-1, 1]^2 whose normal cones are apart
    // from the start, so the pair is never split: z = x / 10 against
    // z = x^2 / 20 - y / 10, meeting along y = x^2 / 2 - x. The window
    // x in [-0.1, 0.1] holds a stretch of it and no edge of either sheet:
    // only the window's own edges can seed it.
    let a = sheet(-1.0, -1.0, [[-0.1; 3], [0.0; 3], [0.1; 3]]);
    let c = [1.0, -1.0, 1.0];
    let mut h = [[0.0; 3]; 3];
    for i in 0..3 {
        for j in 0..3 {
            h[i][j] = 0.05 * c[i] - 0.1 * (j as f64 - 1.0);
        }
    }
    let b = sheet(-1.0, -1.0, h);
    let w1 = (Point2::new(0.45, 0.3), Point2::new(0.55, 0.7));
    let w2 = (Point2::new(0.0, 0.0), Point2::new(1.0, 1.0));
    let curves = spline_pair_intersection(&a, &b, Some((w1, w2))).expect("a stretch");
    assert_eq!(curves.len(), 1);
    let s = &curves[0];
    for i in 0..=200 {
        let t = s.end() * i as f64 / 200.0;
        let (uv, _, p) = s.solve(t).unwrap();
        assert!(uv.x >= w1.0.x - 1e-12 && uv.x <= w1.1.x + 1e-12);
        assert!((p.y - (0.5 * p.x * p.x - p.x)).abs() < 1e-9, "{p:?}");
    }
    // From one side of the window to the other.
    let ends = [s.nodes[0].first.x, s.nodes[s.nodes.len() - 1].first.x];
    assert!(
        (ends[0].min(ends[1]) - 0.45).abs() < 1e-12 && (ends[0].max(ends[1]) - 0.55).abs() < 1e-12,
        "{ends:?}"
    );
}

#[test]
fn a_window_edge_seed_heading_out_ends_where_it_is() {
    // The loop x^2 + y^2 = 1/2 on the bowl leaves the window u >= 0.8
    // (x >= 0.6) through its edge: the seeds there, marched outward, end
    // at once.
    let (b, f) = (bowl(), flat(0.5));
    let w1 = (Point2::new(0.8, 0.3), Point2::new(0.9, 0.7));
    let w2 = (Point2::new(0.0, 0.0), Point2::new(1.0, 1.0));
    let curves = spline_pair_intersection(&b, &f, Some((w1, w2))).expect("a stretch");
    assert_eq!(curves.len(), 1);
    for i in 0..=200 {
        let t = curves[0].end() * i as f64 / 200.0;
        let (a, _, p) = curves[0].solve(t).unwrap();
        assert!(a.x >= w1.0.x - 1e-12 && a.x <= w1.1.x + 1e-12);
        assert!((p.x.hypot(p.y) - 0.5f64.sqrt()).abs() < 1e-9 && (p.z - 0.5).abs() < 1e-9);
    }
}

#[test]
fn a_bowl_touching_a_flat_sheet_is_refused_not_guessed() {
    // The bowl's lowest point touches the plane z = 0: the section is one
    // point, where the surfaces' normals agree and no loop-free split
    // exists.
    assert_eq!(
        spline_pair_intersection(&bowl(), &flat(0.0), None),
        Err(ExactIntersectionRefusal::NotRegularCurve)
    );
}
