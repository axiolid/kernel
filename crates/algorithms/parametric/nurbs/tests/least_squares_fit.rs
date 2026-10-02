//! Least-squares curve and surface approximation (#150, ledger B11).
//!
//! Approximation is a different contract from `interpolate_curve3`: the
//! fitted curve or surface minimises deviation rather than passing exactly
//! through every point, so these tests check the deviation the fit itself
//! reports against an independent evaluator (`axiolid_evaluate`), never
//! against the fitter's own bookkeeping.

use axiolid_contracts::GeomError;
use axiolid_core::{Point3, Scalar};
use axiolid_curve::{BSplineCurve3, KnotSpec};
use axiolid_evaluate::curve::bspline_jet3;
use axiolid_evaluate::surface::bspline_jet;
use axiolid_nurbs::{
    fit_curve3, fit_curve3_to_tolerance, fit_surface_grid, CurveFitOptions, CurveFitTolerance,
    Parameterisation, SurfaceFitOptions,
};
use axiolid_surface::BSplineSurface;

fn sample_curve(curve: &BSplineCurve3, t: Scalar) -> Point3 {
    bspline_jet3(curve, t).expect("valid sample").point
}

fn sample_surface(surface: &BSplineSurface, u: Scalar, v: Scalar) -> Point3 {
    bspline_jet(surface, u, v).expect("valid sample").point
}

/// A single-span (Bezier) cubic: no interior knots, so the approximation
/// knot-placement formula and the curve's own knot vector are both just the
/// clamped endpoints, regardless of parameterisation. That equality is what
/// makes exact reproduction possible below.
fn bezier_cubic() -> BSplineCurve3 {
    BSplineCurve3 {
        degree: 3,
        control_points: vec![
            Point3::new(0.0, 0.0, 0.0),
            Point3::new(1.0, 3.0, 0.5),
            Point3::new(3.0, 2.0, -1.0),
            Point3::new(4.0, 0.0, 1.0),
        ],
        knots: vec![0.0, 1.0],
        multiplicities: vec![4, 4],
        weights: None,
        knot_spec: KnotSpec::Unspecified,
        closed: false,
        self_intersect: None,
    }
}

fn bicubic_patch() -> BSplineSurface {
    let row = |z: Scalar| {
        vec![
            Point3::new(0.0, 0.0, z),
            Point3::new(1.0, 0.0, z + 0.5),
            Point3::new(2.0, 0.0, z - 0.3),
            Point3::new(3.0, 0.0, z + 0.2),
        ]
    };
    let control_points = vec![
        row(0.0)
            .into_iter()
            .map(|p| Point3::new(p.x, 0.0, p.z))
            .collect::<Vec<_>>(),
        (0..4)
            .map(|i| Point3::new(i as Scalar, 1.0, 0.4 + 0.2 * i as Scalar))
            .collect(),
        (0..4)
            .map(|i| Point3::new(i as Scalar, 2.0, -0.3 + 0.1 * i as Scalar))
            .collect(),
        (0..4)
            .map(|i| Point3::new(i as Scalar, 3.0, 0.6 - 0.15 * i as Scalar))
            .collect(),
    ];
    BSplineSurface {
        u_degree: 3,
        v_degree: 3,
        control_points,
        u_knots: vec![0.0, 1.0],
        u_multiplicities: vec![4, 4],
        v_knots: vec![0.0, 1.0],
        v_multiplicities: vec![4, 4],
        weights: None,
        u_closed: false,
        v_closed: false,
        knot_spec: KnotSpec::Unspecified,
        self_intersect: None,
    }
}

#[test]
fn fitting_a_bezier_sample_reproduces_the_curve_to_near_machine_precision() {
    let original = bezier_cubic();
    let samples = 11;
    let points: Vec<Point3> = (0..=samples)
        .map(|k| sample_curve(&original, k as Scalar / samples as Scalar))
        .collect();

    let options = CurveFitOptions {
        degree: 3,
        control_point_count: 4,
        parameterisation: Parameterisation::Uniform,
        interpolate_endpoints: false,
        smoothing: 0.0,
        parameter_correction_iterations: 0,
    };
    let fit = fit_curve3(&points, &options).expect("well-posed fit");

    assert!(
        fit.max_deviation < 1e-9,
        "max deviation too large: {}",
        fit.max_deviation
    );
    assert!(fit.rms_deviation <= fit.max_deviation);

    // Cross-check against an independent evaluator at parameters other than
    // the fitting nodes, not just the fitter's own report.
    for k in 0..=20 {
        let t = k as Scalar / 20.0;
        let expected = sample_curve(&original, t);
        let actual = sample_curve(&fit.curve, t);
        assert!(
            (actual - expected).length() < 1e-8,
            "curve diverges at t={t}: expected {expected:?}, got {actual:?}"
        );
    }
}

#[test]
fn surface_fit_of_a_bicubic_sample_reproduces_the_patch() {
    let original = bicubic_patch();
    let n = 8;
    let points: Vec<Vec<Point3>> = (0..=n)
        .map(|i| {
            (0..=n)
                .map(|j| {
                    sample_surface(
                        &original,
                        i as Scalar / n as Scalar,
                        j as Scalar / n as Scalar,
                    )
                })
                .collect()
        })
        .collect();

    let options = SurfaceFitOptions {
        u_degree: 3,
        v_degree: 3,
        u_control_count: 4,
        v_control_count: 4,
        parameterisation: Parameterisation::Uniform,
    };
    let fit = fit_surface_grid(&points, &options).expect("well-posed surface fit");
    eprintln!(
        "bicubic patch reproduction: max_deviation = {:e}, rms_deviation = {:e}",
        fit.max_deviation, fit.rms_deviation
    );

    // Observed ~3.6e-15 (see the eprintln above); a single-span bicubic
    // patch resampled at a matching degree/control count is essentially an
    // exact reproduction, limited only by rounding through two sequential
    // least-squares solves (row fit, then column fit).
    assert!(
        fit.max_deviation < 5e-14,
        "max deviation too large: {}",
        fit.max_deviation
    );

    for i in 0..=4 {
        for j in 0..=4 {
            let u = i as Scalar / 4.0;
            let v = j as Scalar / 4.0;
            let expected = sample_surface(&original, u, v);
            let actual = sample_surface(&fit.surface, u, v);
            assert!(
                (actual - expected).length() < 1e-12,
                "surface diverges at ({u}, {v}): expected {expected:?}, got {actual:?}"
            );
        }
    }
}

/// A quarter-circle approximated by cubics: the Piegl and Tiller
/// approximation error for a degree-`p` B-spline should shrink quickly as
/// control points are added, since a circular arc is smooth. This checks
/// the trend, not an exact rate constant.
fn quarter_circle_points(count: usize) -> Vec<Point3> {
    (0..count)
        .map(|k| {
            let theta = std::f64::consts::FRAC_PI_2 * (k as Scalar) / (count - 1) as Scalar;
            Point3::new(theta.cos(), theta.sin(), 0.0)
        })
        .collect()
}

#[test]
fn more_control_points_reduce_arc_fit_deviation() {
    let points = quarter_circle_points(60);
    let mut previous = Scalar::INFINITY;
    for control_count in [5usize, 7, 10, 14] {
        let options = CurveFitOptions {
            degree: 3,
            control_point_count: control_count,
            parameterisation: Parameterisation::ChordLength,
            interpolate_endpoints: false,
            smoothing: 0.0,
            parameter_correction_iterations: 0,
        };
        let fit = fit_curve3(&points, &options).expect("arc fit is well-posed");
        assert!(
            fit.max_deviation < previous,
            "deviation did not shrink at {control_count} control points: {} >= {previous}",
            fit.max_deviation
        );
        previous = fit.max_deviation;
    }
}

#[test]
fn constrained_endpoints_hold_exactly() {
    let points = quarter_circle_points(20);
    let options = CurveFitOptions {
        degree: 3,
        control_point_count: 6,
        parameterisation: Parameterisation::ChordLength,
        interpolate_endpoints: true,
        smoothing: 0.0,
        parameter_correction_iterations: 0,
    };
    let fit = fit_curve3(&points, &options).expect("constrained fit is well-posed");

    let start = sample_curve(&fit.curve, 0.0);
    let end = sample_curve(&fit.curve, 1.0);
    assert!(
        (start - points[0]).length() < 1e-9,
        "start moved: {start:?}"
    );
    assert!(
        (end - points[points.len() - 1]).length() < 1e-9,
        "end moved: {end:?}"
    );
}

#[test]
fn parameter_correction_stays_close_to_the_uncorrected_fit() {
    // Newton parameter correction is not guaranteed to monotonically
    // improve the max deviation (it minimises a local, not global,
    // objective per point), but it must not blow up: it stays within the
    // same order of magnitude as the uncorrected fit on a smooth arc.
    let points = quarter_circle_points(30);
    let base = CurveFitOptions {
        degree: 3,
        control_point_count: 7,
        parameterisation: Parameterisation::ChordLength,
        interpolate_endpoints: false,
        smoothing: 0.0,
        parameter_correction_iterations: 0,
    };
    let corrected = CurveFitOptions {
        parameter_correction_iterations: 3,
        ..base.clone()
    };
    let before = fit_curve3(&points, &base).expect("base fit");
    let after = fit_curve3(&points, &corrected).expect("corrected fit");
    assert!(
        after.max_deviation <= before.max_deviation * 3.0,
        "parameter correction made the fit much worse: {} vs {}",
        after.max_deviation,
        before.max_deviation
    );
}

#[test]
fn tolerance_search_meets_the_requested_accuracy() {
    let points = quarter_circle_points(40);
    let tolerance = CurveFitTolerance {
        degree: 3,
        parameterisation: Parameterisation::ChordLength,
        interpolate_endpoints: false,
        max_deviation: 1e-4,
        max_control_points: 40,
    };
    let fit = fit_curve3_to_tolerance(&points, &tolerance).expect("budget is generous enough");
    assert!(fit.max_deviation <= 1e-4);
}

#[test]
fn tolerance_search_refuses_rather_than_silently_missing_the_budget() {
    let points = quarter_circle_points(40);
    let tolerance = CurveFitTolerance {
        degree: 3,
        parameterisation: Parameterisation::ChordLength,
        interpolate_endpoints: false,
        max_deviation: 1e-12,
        max_control_points: 6,
    };
    let err = fit_curve3_to_tolerance(&points, &tolerance).expect_err("budget is too small");
    assert!(matches!(err, GeomError::BudgetExceeded { .. }));
}

#[test]
fn too_few_points_for_the_degree_is_refused() {
    let points = vec![Point3::new(0.0, 0.0, 0.0), Point3::new(1.0, 0.0, 0.0)];
    let options = CurveFitOptions {
        degree: 3,
        control_point_count: 4,
        parameterisation: Parameterisation::ChordLength,
        ..CurveFitOptions::default()
    };
    let err = fit_curve3(&points, &options).expect_err("two points cannot fit four cubic controls");
    assert!(matches!(err, GeomError::InvalidInput(_)));
}

#[test]
fn too_many_control_points_for_the_data_is_refused() {
    let points = quarter_circle_points(5);
    let options = CurveFitOptions {
        degree: 3,
        control_point_count: 6,
        parameterisation: Parameterisation::ChordLength,
        ..CurveFitOptions::default()
    };
    let err = fit_curve3(&points, &options).expect_err("more controls than points is never a fit");
    assert!(matches!(err, GeomError::InvalidInput(_)));
}

#[test]
fn non_finite_points_are_refused() {
    let points = vec![
        Point3::new(0.0, 0.0, 0.0),
        Point3::new(f64::NAN, 1.0, 0.0),
        Point3::new(2.0, 0.0, 0.0),
        Point3::new(3.0, 1.0, 0.0),
        Point3::new(4.0, 0.0, 0.0),
    ];
    let options = CurveFitOptions {
        degree: 3,
        control_point_count: 4,
        ..CurveFitOptions::default()
    };
    let err = fit_curve3(&points, &options).expect_err("NaN point must be refused");
    assert!(matches!(err, GeomError::InvalidInput(_)));
}

#[test]
fn coincident_consecutive_points_under_chord_length_are_refused() {
    let points = vec![
        Point3::new(0.0, 0.0, 0.0),
        Point3::new(1.0, 0.0, 0.0),
        Point3::new(1.0, 0.0, 0.0),
        Point3::new(2.0, 1.0, 0.0),
        Point3::new(3.0, 0.0, 0.0),
    ];
    let options = CurveFitOptions {
        degree: 3,
        control_point_count: 4,
        parameterisation: Parameterisation::ChordLength,
        ..CurveFitOptions::default()
    };
    let err = fit_curve3(&points, &options).expect_err("coincident points undefine chord length");
    assert!(matches!(err, GeomError::Degenerate(_)));
}

#[test]
fn invalid_degree_zero_is_refused() {
    let points = quarter_circle_points(6);
    let options = CurveFitOptions {
        degree: 0,
        control_point_count: 4,
        ..CurveFitOptions::default()
    };
    let err = fit_curve3(&points, &options).expect_err("degree zero is not a curve");
    assert!(matches!(err, GeomError::InvalidInput(_)));
}

#[test]
fn rank_deficient_system_is_refused() {
    // Eight of nine points sit almost on top of each other at the chord's
    // start, with only the last point far away: every interior knot the
    // averaging formula places then falls in the crowded region, leaving a
    // control point near the far end with no data in its support, so the
    // design matrix loses column rank.
    let mut points = vec![Point3::new(0.0, 0.0, 0.0)];
    for k in 1..9 {
        points.push(Point3::new(k as Scalar * 1e-9, 0.0, 0.0));
    }
    points.push(Point3::new(100.0, 1.0, 0.0));
    let options = CurveFitOptions {
        degree: 3,
        control_point_count: 7,
        parameterisation: Parameterisation::ChordLength,
        ..CurveFitOptions::default()
    };
    let err = fit_curve3(&points, &options)
        .expect_err("a control point with no data in its support cannot be determined");
    assert!(matches!(err, GeomError::Degenerate(_)));
    let message = err.to_string();
    assert!(
        message.contains("numerical rank"),
        "refusal dropped the rank diagnostic: {message}"
    );
}

#[test]
fn coincident_consecutive_points_under_centripetal_are_refused() {
    let points = vec![
        Point3::new(0.0, 0.0, 0.0),
        Point3::new(1.0, 0.0, 0.0),
        Point3::new(1.0, 0.0, 0.0),
        Point3::new(2.0, 1.0, 0.0),
        Point3::new(3.0, 0.0, 0.0),
    ];
    let options = CurveFitOptions {
        degree: 3,
        control_point_count: 4,
        parameterisation: Parameterisation::Centripetal,
        ..CurveFitOptions::default()
    };
    let err =
        fit_curve3(&points, &options).expect_err("coincident points undefine centripetal spacing");
    assert!(matches!(err, GeomError::Degenerate(_)));
}

#[test]
fn negative_smoothing_weight_is_refused() {
    let points = quarter_circle_points(8);
    let options = CurveFitOptions {
        degree: 3,
        control_point_count: 5,
        smoothing: -1.0,
        ..CurveFitOptions::default()
    };
    let err = fit_curve3(&points, &options).expect_err("negative smoothing weight is invalid");
    assert!(matches!(err, GeomError::InvalidInput(_)));
}

#[test]
fn tolerance_search_refuses_at_the_full_point_budget_too() {
    // `max_control_points` is generous enough to never trip on its own, so
    // this exercises the other refusal: running out of points to add
    // without having met the tolerance.
    let points = quarter_circle_points(10);
    let tolerance = CurveFitTolerance {
        degree: 3,
        parameterisation: Parameterisation::ChordLength,
        interpolate_endpoints: false,
        max_deviation: 0.0,
        max_control_points: 1000,
    };
    let err = fit_curve3_to_tolerance(&points, &tolerance)
        .expect_err("a quarter circle is never fit to exactly zero deviation by a B-spline");
    assert!(matches!(err, GeomError::BudgetExceeded { .. }));
}

#[test]
fn surface_fit_refuses_ragged_rows() {
    // Control counts valid for a 3 x 4 grid, so the ragged-row check is the
    // only thing standing between this input and an out-of-bounds column
    // read (row 1 is one point short of the others).
    let points = vec![
        vec![
            Point3::new(0.0, 0.0, 0.0),
            Point3::new(1.0, 0.0, 0.0),
            Point3::new(2.0, 0.0, 0.0),
            Point3::new(3.0, 0.0, 0.0),
        ],
        vec![
            Point3::new(0.0, 1.0, 0.0),
            Point3::new(1.0, 1.0, 0.0),
            Point3::new(2.0, 1.0, 0.0),
        ],
        vec![
            Point3::new(0.0, 2.0, 0.0),
            Point3::new(1.0, 2.0, 0.0),
            Point3::new(2.0, 2.0, 0.0),
            Point3::new(3.0, 2.0, 0.0),
        ],
    ];
    let options = SurfaceFitOptions {
        u_degree: 1,
        v_degree: 1,
        u_control_count: 3,
        v_control_count: 3,
        parameterisation: Parameterisation::Uniform,
    };
    let err = fit_surface_grid(&points, &options).expect_err("ragged grid must be refused");
    assert!(matches!(err, GeomError::InvalidInput(_)));
}
