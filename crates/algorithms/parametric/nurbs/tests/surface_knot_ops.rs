//! Surface knot removal, degree change and iso-curves (#141).
//!
//! Oracles: removal undoes insertion exactly; elevation and iso-curves
//! agree with evaluation of the original surface on dense grids that hit
//! the knots and the domain ends; a lossy removal or reduction is refused
//! within a small tolerance, and where it is accepted the stated bound
//! covers the deviation measured on a grid.

use axiolid_contracts::{GeomError, Operation};
use axiolid_core::Point3;
use axiolid_curve::KnotSpec;
use axiolid_evaluate::curve::bspline_jet3;
use axiolid_evaluate::surface::bspline_jet;
use axiolid_nurbs::{
    elevate_surface_degree_u, elevate_surface_degree_v, insert_surface_knot_u,
    insert_surface_knot_v, iso_curve_at_u, iso_curve_at_v, reduce_surface_degree_u,
    reduce_surface_degree_v, remove_surface_knot_u, remove_surface_knot_v,
};
use axiolid_surface::BSplineSurface;

/// A bicubic-by-quadratic patch with an interior knot on each axis and a
/// net with no hidden structure; `rational` gives it uneven weights.
fn patch(rational: bool) -> BSplineSurface {
    let control_points: Vec<Vec<Point3>> = (0..5)
        .map(|i| {
            (0..4)
                .map(|j| {
                    let (x, y) = (f64::from(i), f64::from(j));
                    Point3::new(x, y, ((3 * i + 5 * j) % 7) as f64 * 0.25 - 0.5)
                })
                .collect()
        })
        .collect();
    let weights = rational.then(|| {
        (0..5)
            .map(|i| {
                (0..4)
                    .map(|j| 0.8 + f64::from((2 * i + 3 * j) % 5) * 0.1)
                    .collect()
            })
            .collect()
    });
    BSplineSurface {
        u_degree: 3,
        v_degree: 2,
        control_points,
        u_knots: vec![0.0, 0.4, 1.0],
        u_multiplicities: vec![4, 1, 4],
        v_knots: vec![0.0, 0.5, 1.0],
        v_multiplicities: vec![3, 1, 3],
        weights,
        knot_spec: KnotSpec::Unspecified,
        u_closed: false,
        v_closed: false,
        self_intersect: None,
    }
}

/// Parameters from 0 to 1 in steps of 1/40: they include 0.4, 0.5 and both
/// ends.
fn grid() -> Vec<f64> {
    (0..=40).map(|k| f64::from(k) / 40.0).collect()
}

fn largest_gap(a: &BSplineSurface, b: &BSplineSurface) -> f64 {
    let mut worst: f64 = 0.0;
    for &u in &grid() {
        for &v in &grid() {
            let p = bspline_jet(a, u, v).unwrap().point;
            let q = bspline_jet(b, u, v).unwrap().point;
            worst = worst.max(p.distance(q));
        }
    }
    worst
}

fn same_net(a: &BSplineSurface, b: &BSplineSurface) {
    assert_eq!(a.u_knots, b.u_knots);
    assert_eq!(a.u_multiplicities, b.u_multiplicities);
    assert_eq!(a.v_knots, b.v_knots);
    assert_eq!(a.v_multiplicities, b.v_multiplicities);
    assert_eq!(a.control_points.len(), b.control_points.len());
    for (ra, rb) in a.control_points.iter().zip(&b.control_points) {
        assert_eq!(ra.len(), rb.len());
        for (p, q) in ra.iter().zip(rb) {
            assert!(p.distance(*q) < 1e-12, "{p:?} vs {q:?}");
        }
    }
    match (&a.weights, &b.weights) {
        (None, None) => {}
        (Some(wa), Some(wb)) => {
            for (ra, rb) in wa.iter().zip(wb) {
                for (x, y) in ra.iter().zip(rb) {
                    assert!((x - y).abs() < 1e-12, "{x} vs {y}");
                }
            }
        }
        _ => panic!("rationality changed"),
    }
}

#[test]
fn removal_undoes_insertion_along_both_axes() {
    for rational in [false, true] {
        let surface = patch(rational);
        let once = insert_surface_knot_u(&surface, 0.7).unwrap();
        let twice = insert_surface_knot_u(&once, 0.7).unwrap();
        // Asked for five, the knot runs out after two.
        let removal = remove_surface_knot_u(&twice, 0.7, 5, 1e-9).unwrap();
        assert_eq!(removal.removed, 2);
        assert!(removal.deviation_upper_bound < 1e-12);
        same_net(&removal.surface, &surface);

        let inserted = insert_surface_knot_v(&surface, 0.25).unwrap();
        let removal = remove_surface_knot_v(&inserted, 0.25, 1, 1e-9).unwrap();
        assert_eq!(removal.removed, 1);
        same_net(&removal.surface, &surface);
    }
}

#[test]
fn a_knot_that_carries_shape_is_kept() {
    for rational in [false, true] {
        let surface = patch(rational);
        for removal in [
            remove_surface_knot_u(&surface, 0.4, 1, 1e-6).unwrap(),
            remove_surface_knot_v(&surface, 0.5, 1, 1e-6).unwrap(),
        ] {
            assert_eq!(removal.removed, 0);
            assert_eq!(removal.deviation_upper_bound, 0.0);
            assert_eq!(removal.surface, surface);
        }
    }
}

#[test]
fn an_accepted_lossy_removal_stays_within_its_bound() {
    for rational in [false, true] {
        let surface = patch(rational);
        for removal in [
            remove_surface_knot_u(&surface, 0.4, 1, 100.0).unwrap(),
            remove_surface_knot_v(&surface, 0.5, 1, 100.0).unwrap(),
        ] {
            assert_eq!(removal.removed, 1);
            let measured = largest_gap(&surface, &removal.surface);
            assert!(measured > 1e-3, "the knot carried shape");
            assert!(
                measured <= removal.deviation_upper_bound,
                "rational {rational}: measured {measured} above the bound {}",
                removal.deviation_upper_bound
            );
            // A tolerance just under the bound refuses the removal.
            let tight = removal.deviation_upper_bound * 0.999;
            let refused = if removal.surface.u_knots.len() == 2 {
                remove_surface_knot_u(&surface, 0.4, 1, tight).unwrap()
            } else {
                remove_surface_knot_v(&surface, 0.5, 1, tight).unwrap()
            };
            assert_eq!(refused.removed, 0);
        }
    }
}

#[test]
fn bad_removal_requests_are_errors() {
    let surface = patch(false);
    for (parameter, tolerance) in [(0.3, 1.0), (0.0, 1.0), (1.0, 1.0), (0.4, -1.0)] {
        assert!(matches!(
            remove_surface_knot_u(&surface, parameter, 1, tolerance),
            Err(GeomError::InvalidInput(_))
        ));
    }
}

#[test]
fn elevation_keeps_the_surface_in_bezier_form() {
    for rational in [false, true] {
        let surface = patch(rational);
        let u = elevate_surface_degree_u(&surface).unwrap();
        assert_eq!((u.u_degree, u.v_degree), (4, 2));
        assert_eq!(u.u_multiplicities, vec![5, 4, 5]);
        assert_eq!(u.control_points.len(), 9);
        assert_eq!(u.v_multiplicities, surface.v_multiplicities);
        let v = elevate_surface_degree_v(&surface).unwrap();
        assert_eq!((v.u_degree, v.v_degree), (3, 3));
        assert_eq!(v.v_multiplicities, vec![4, 3, 4]);
        assert_eq!(v.control_points[0].len(), 7);
        assert_eq!(u.weights.is_some(), rational);
        // The net spans about 4 units: 1e-12 relative.
        assert!(largest_gap(&surface, &u) < 4e-12);
        assert!(largest_gap(&surface, &v) < 4e-12);
    }
}

#[test]
fn reduction_recovers_an_elevated_surface_and_refuses_otherwise() {
    let surface = patch(false);
    let back = reduce_surface_degree_u(&elevate_surface_degree_u(&surface).unwrap(), 1e-9).unwrap();
    assert_eq!(back.surface.u_degree, 3);
    assert!(back.deviation_upper_bound < 1e-12);
    assert!(largest_gap(&surface, &back.surface) < 4e-12);
    let back = reduce_surface_degree_v(&elevate_surface_degree_v(&surface).unwrap(), 1e-9).unwrap();
    assert_eq!(back.surface.v_degree, 2);
    assert!(largest_gap(&surface, &back.surface) < 4e-12);

    // This patch needs its degree on both axes.
    for result in [
        reduce_surface_degree_u(&surface, 1e-6),
        reduce_surface_degree_v(&surface, 1e-6),
    ] {
        assert!(
            matches!(result, Err(GeomError::Degenerate(_))),
            "{result:?}"
        );
    }
    // Where a lossy reduction is accepted, its bound covers the deviation.
    let loose = reduce_surface_degree_u(&surface, 100.0).unwrap();
    let measured = largest_gap(&surface, &loose.surface);
    assert!(measured > 1e-3 && measured <= loose.deviation_upper_bound);

    // Refused for the surface, not only by the curve step beneath it.
    assert!(matches!(
        reduce_surface_degree_u(&patch(true), 1.0),
        Err(GeomError::Unsupported {
            operation: Operation::SurfaceEvaluation,
            ..
        })
    ));
}

/// A degenerate strip, both columns the segment from `x = 0` to `x = 2`,
/// linear in `u` with a knot at 0.5, whose middle point lies on the chord
/// but weighs `heavy` times the ends. Removing the knot keeps the image and
/// changes only the parametrisation: at `u = 1/4` the original is at
/// `x = heavy / (1 + heavy)`, the result at `x = 1/2`.
fn reparametrised_strip(scale: f64) -> BSplineSurface {
    let heavy = 10.0;
    BSplineSurface {
        u_degree: 1,
        v_degree: 1,
        control_points: (0..3)
            .map(|i| vec![Point3::new(f64::from(i), 0.0, 0.0); 2])
            .collect(),
        u_knots: vec![0.0, 0.5, 1.0],
        u_multiplicities: vec![2, 1, 2],
        v_knots: vec![0.0, 1.0],
        v_multiplicities: vec![2, 2],
        weights: Some(vec![
            vec![scale, scale],
            vec![heavy * scale, heavy * scale],
            vec![scale, scale],
        ]),
        knot_spec: KnotSpec::Unspecified,
        u_closed: false,
        v_closed: false,
        self_intersect: None,
    }
}

#[test]
fn the_bound_covers_a_change_made_only_by_the_weights() {
    // The homogeneous points move with the weights, so every difference is
    // a weight difference: a bound read off positions alone would be zero.
    // The weights' scale must not matter either, since the surface does not
    // change when all weights are scaled together.
    let mut bounds = Vec::new();
    for scale in [1.0, 0.01] {
        let strip = reparametrised_strip(scale);
        let removal = remove_surface_knot_u(&strip, 0.5, 1, 100.0).unwrap();
        assert_eq!(removal.removed, 1);
        let at = |s: &BSplineSurface| bspline_jet(s, 0.25, 0.5).unwrap().point;
        assert!((at(&strip).x - 10.0 / 11.0).abs() < 1e-12);
        assert!((at(&removal.surface).x - 0.5).abs() < 1e-12);
        let measured = largest_gap(&strip, &removal.surface);
        assert!(measured >= 10.0 / 11.0 - 0.5 - 1e-12);
        assert!(
            measured <= removal.deviation_upper_bound,
            "scale {scale}: measured {measured} above the bound {}",
            removal.deviation_upper_bound
        );
        bounds.push(removal.deviation_upper_bound);
    }
    assert!((bounds[0] - bounds[1]).abs() <= 1e-9 * bounds[0]);
}

#[test]
fn removal_stops_at_the_requested_count() {
    let surface = patch(true);
    let once = insert_surface_knot_u(&surface, 0.7).unwrap();
    let twice = insert_surface_knot_u(&once, 0.7).unwrap();
    let removal = remove_surface_knot_u(&twice, 0.7, 1, 1e-9).unwrap();
    assert_eq!(removal.removed, 1);
    assert_eq!(removal.surface.u_knots, vec![0.0, 0.4, 0.7, 1.0]);
    assert_eq!(removal.surface.u_multiplicities, vec![4, 1, 1, 4]);
    assert!(largest_gap(&surface, &removal.surface) < 4e-12);
}

#[test]
fn iso_curves_match_the_surface() {
    for rational in [false, true] {
        let surface = patch(rational);
        for &fixed in &[0.0, 0.123, 0.4, 0.5, 0.77, 1.0] {
            let along_v = iso_curve_at_u(&surface, fixed).unwrap();
            let along_u = iso_curve_at_v(&surface, fixed).unwrap();
            assert_eq!(along_v.degree, 2);
            assert_eq!(along_u.degree, 3);
            assert_eq!(along_v.weights.is_some(), rational);
            for &t in &grid() {
                let on_v = bspline_jet3(&along_v, t).unwrap().point;
                let on_u = bspline_jet3(&along_u, t).unwrap().point;
                let s_v = bspline_jet(&surface, fixed, t).unwrap().point;
                let s_u = bspline_jet(&surface, t, fixed).unwrap().point;
                assert!(on_v.distance(s_v) < 4e-12, "u = {fixed}, v = {t}");
                assert!(on_u.distance(s_u) < 4e-12, "v = {fixed}, u = {t}");
            }
        }
        for outside in [-0.01, 1.01, f64::NAN] {
            assert!(iso_curve_at_u(&surface, outside).is_err());
            assert!(iso_curve_at_v(&surface, outside).is_err());
        }
    }
}
