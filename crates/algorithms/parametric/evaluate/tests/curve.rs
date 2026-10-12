//! Curve evaluation checked against closed-form truth (ADR 0012).
//!
//! Every assertion here is a value derivable on paper: arc length by numeric
//! quadrature of the analytic derivative, curvature by the standard formula,
//! exact points at parameters where the answer is known by symmetry.
//!
//! The point is that none of these compare the evaluator to itself.

use axiolid_contracts::GeomError;
use axiolid_core::{Frame2, Frame3, Interval, Point2, Point3, Scalar, Vec2, Vec3};
use axiolid_curve::{
    BSplineCurve2, Circle2, Circle3, Curve2, Curve3, Ellipse2, KnotSpec, Line2, Polyline2,
};
use axiolid_evaluate::curve::{derivative2, derivative3, domain2, evaluate2, evaluate3, flatten2};

const TAU: Scalar = core::f64::consts::TAU;

fn frame2() -> Frame2 {
    Frame2 {
        origin: Point2::new(0.0, 0.0),
        x: Vec2::new(1.0, 0.0),
        y: Vec2::new(0.0, 1.0),
    }
}

fn unit_circle() -> Curve2 {
    Curve2::Circle(Circle2 {
        frame: frame2(),
        radius: 1.0,
    })
}

/// Arc length by Gauss-Legendre quadrature of |C'(t)|.
///
/// Independent of the evaluator's *position* code: it integrates the
/// derivative, so agreement between the two is real evidence.
fn arc_length(curve: &Curve2, a: Scalar, b: Scalar, panels: usize) -> Scalar {
    // 4-point Gauss-Legendre on [-1, 1].
    const X: [Scalar; 4] = [
        -0.861_136_311_594_052_6,
        -0.339_981_043_584_856_3,
        0.339_981_043_584_856_3,
        0.861_136_311_594_052_6,
    ];
    const W: [Scalar; 4] = [
        0.347_854_845_137_453_9,
        0.652_145_154_862_546_1,
        0.652_145_154_862_546_1,
        0.347_854_845_137_453_9,
    ];
    let h = (b - a) / panels as Scalar;
    let mut total = 0.0;
    for k in 0..panels {
        let lo = a + h * k as Scalar;
        let mid = lo + h * 0.5;
        for (x, w) in X.iter().zip(W.iter()) {
            let t = mid + x * h * 0.5;
            total += w * derivative2(curve, t).unwrap().length() * h * 0.5;
        }
    }
    total
}

// --- exact points -----------------------------------------------------------

#[test]
fn a_circle_hits_its_quadrant_points_exactly() {
    let c = unit_circle();
    for (t, want) in [
        (0.0, Point2::new(1.0, 0.0)),
        (TAU * 0.25, Point2::new(0.0, 1.0)),
        (TAU * 0.5, Point2::new(-1.0, 0.0)),
        (TAU * 0.75, Point2::new(0.0, -1.0)),
    ] {
        let got = evaluate2(&c, t).unwrap();
        assert!(
            (got - want).length() < 1e-15,
            "circle at {t}: want {want:?}, got {got:?}"
        );
    }
}

#[test]
fn a_circle_stays_on_its_radius() {
    let c = Curve2::Circle(Circle2 {
        frame: Frame2 {
            origin: Point2::new(3.0, -2.0),
            ..frame2()
        },
        radius: 2.5,
    });
    for k in 0..64 {
        let t = TAU * k as Scalar / 64.0;
        let p = evaluate2(&c, t).unwrap();
        let r = (p - Point2::new(3.0, -2.0)).length();
        assert!((r - 2.5).abs() < 1e-14, "radius drift at {t}: {r}");
    }
}

#[test]
fn a_circle_tangent_is_perpendicular_to_its_radius() {
    let c = unit_circle();
    for k in 0..32 {
        let t = TAU * k as Scalar / 32.0;
        let p = evaluate2(&c, t).unwrap();
        let d = derivative2(&c, t).unwrap();
        // Radius . tangent == 0 for a circle, at every parameter.
        assert!(p.dot(d).abs() < 1e-14, "not perpendicular at {t}");
    }
}

// --- arc length -------------------------------------------------------------

#[test]
fn a_unit_circle_has_circumference_tau() {
    let len = arc_length(&unit_circle(), 0.0, TAU, 64);
    assert!((len - TAU).abs() < 1e-12, "circumference {len}, want {TAU}");
}

#[test]
fn a_circle_arc_length_scales_with_radius() {
    let c = Curve2::Circle(Circle2 {
        frame: frame2(),
        radius: 7.0,
    });
    let quarter = arc_length(&c, 0.0, TAU * 0.25, 64);
    let want = 7.0 * TAU * 0.25;
    assert!((quarter - want).abs() < 1e-11, "arc {quarter}, want {want}");
}

#[test]
fn an_ellipse_matches_its_known_perimeter() {
    // a=2, b=1: perimeter = 4*a*E(e) with e^2 = 1 - b^2/a^2.
    //
    // The literal below was verified independently of this crate: a periodic
    // trapezoid rule over |C'(t)| converges spectrally for a smooth closed
    // curve, and agrees to machine precision at 2e6 samples. An earlier
    // hand-recalled value (...216130) was wrong in the 9th digit and the
    // evaluator was correct -- hence the provenance note.
    let e = Curve2::Ellipse(Ellipse2 {
        frame: frame2(),
        semi_axis_x: 2.0,
        semi_axis_y: 1.0,
    });
    let len = arc_length(&e, 0.0, TAU, 256);
    assert!(
        (len - 9.688_448_220_547_675).abs() < 1e-12,
        "ellipse perimeter {len}"
    );
}

// --- curvature --------------------------------------------------------------

/// Curvature by the analytic formula, using a central second difference of the
/// *derivative* function (not of positions).
fn curvature(curve: &Curve2, t: Scalar) -> Scalar {
    let h = 1e-6;
    let d = derivative2(curve, t).unwrap();
    let d2 = (derivative2(curve, t + h).unwrap() - derivative2(curve, t - h).unwrap()) / (2.0 * h);
    let cross = d.x * d2.y - d.y * d2.x;
    cross.abs() / d.length().powi(3)
}

#[test]
fn a_circle_has_constant_curvature_one_over_radius() {
    for r in [0.5, 1.0, 4.0] {
        let c = Curve2::Circle(Circle2 {
            frame: frame2(),
            radius: r,
        });
        for k in 0..16 {
            let t = TAU * k as Scalar / 16.0;
            let got = curvature(&c, t);
            assert!(
                (got - 1.0 / r).abs() < 1e-6,
                "curvature {got} at r={r}, want {}",
                1.0 / r
            );
        }
    }
}

#[test]
fn an_ellipse_is_most_curved_at_its_major_axis_end() {
    let e = Curve2::Ellipse(Ellipse2 {
        frame: frame2(),
        semi_axis_x: 3.0,
        semi_axis_y: 1.0,
    });
    // At t=0 (major-axis end) curvature is a/b^2; at t=pi/2 it is b/a^2.
    let at_major = curvature(&e, 0.0);
    let at_minor = curvature(&e, TAU * 0.25);
    assert!((at_major - 3.0 / 1.0).abs() < 1e-5, "major {at_major}");
    assert!((at_minor - 1.0 / 9.0).abs() < 1e-5, "minor {at_minor}");
    assert!(at_major > at_minor);
}

// --- lines and polylines ----------------------------------------------------

#[test]
fn a_line_is_exact_and_has_constant_derivative() {
    let l = Curve2::Line(Line2 {
        origin: Point2::new(1.0, 2.0),
        direction: Vec2::new(3.0, -4.0),
    });
    assert_eq!(evaluate2(&l, 0.0).unwrap(), Point2::new(1.0, 2.0));
    assert_eq!(evaluate2(&l, 2.0).unwrap(), Point2::new(7.0, -6.0));
    // The parameter direction is preserved, including its non-unit length.
    assert_eq!(derivative2(&l, 0.0).unwrap(), Vec2::new(3.0, -4.0));
    assert_eq!(derivative2(&l, 99.0).unwrap(), Vec2::new(3.0, -4.0));
}

#[test]
fn a_polyline_interpolates_its_vertices() {
    let p = Curve2::Polyline(Polyline2 {
        points: vec![
            Point2::new(0.0, 0.0),
            Point2::new(2.0, 0.0),
            Point2::new(2.0, 3.0),
        ],
        closed: false,
    });
    assert_eq!(evaluate2(&p, 0.0).unwrap(), Point2::new(0.0, 0.0));
    assert_eq!(evaluate2(&p, 0.5).unwrap(), Point2::new(1.0, 0.0));
    assert_eq!(evaluate2(&p, 1.0).unwrap(), Point2::new(2.0, 0.0));
    assert_eq!(evaluate2(&p, 1.5).unwrap(), Point2::new(2.0, 1.5));
    // The final parameter is the last vertex, not one past the end.
    assert_eq!(evaluate2(&p, 2.0).unwrap(), Point2::new(2.0, 3.0));
}

#[test]
fn a_closed_polyline_returns_to_its_start() {
    let p = Curve2::Polyline(Polyline2 {
        points: vec![
            Point2::new(0.0, 0.0),
            Point2::new(1.0, 0.0),
            Point2::new(0.0, 1.0),
        ],
        closed: true,
    });
    // Three segments when closed, so t=3 wraps to the origin.
    assert_eq!(evaluate2(&p, 3.0).unwrap(), Point2::new(0.0, 0.0));
}

// --- B-spline ---------------------------------------------------------------

/// Clamped cubic through 4 control points == a Bezier curve.
fn cubic_bezier(pts: [Point2; 4]) -> BSplineCurve2 {
    BSplineCurve2 {
        degree: 3,
        control_points: pts.to_vec(),
        knots: vec![0.0, 1.0],
        multiplicities: vec![4, 4],
        weights: None,
        closed: false,
        self_intersect: None,
        knot_spec: KnotSpec::PiecewiseBezier,
    }
}

#[test]
fn a_clamped_bspline_interpolates_its_end_control_points() {
    let pts = [
        Point2::new(0.0, 0.0),
        Point2::new(1.0, 2.0),
        Point2::new(3.0, 2.0),
        Point2::new(4.0, 0.0),
    ];
    let c = Curve2::BSpline(cubic_bezier(pts));
    assert!((evaluate2(&c, 0.0).unwrap() - pts[0]).length() < 1e-15);
    assert!((evaluate2(&c, 1.0).unwrap() - pts[3]).length() < 1e-15);
}

#[test]
fn a_cubic_bspline_matches_the_bernstein_form() {
    let pts = [
        Point2::new(0.0, 0.0),
        Point2::new(1.0, 2.0),
        Point2::new(3.0, 2.0),
        Point2::new(4.0, 0.0),
    ];
    let c = Curve2::BSpline(cubic_bezier(pts));
    for k in 0..=20 {
        let t = k as Scalar / 20.0;
        let u = 1.0 - t;
        // Independent Bernstein evaluation.
        let want = pts[0] * (u * u * u)
            + pts[1] * (3.0 * u * u * t)
            + pts[2] * (3.0 * u * t * t)
            + pts[3] * (t * t * t);
        let got = evaluate2(&c, t).unwrap();
        assert!(
            (got - want).length() < 1e-13,
            "de Boor vs Bernstein at {t}: {got:?} vs {want:?}"
        );
    }
}

#[test]
fn a_bspline_derivative_matches_the_bernstein_hodograph() {
    let pts = [
        Point2::new(0.0, 0.0),
        Point2::new(1.0, 2.0),
        Point2::new(3.0, 2.0),
        Point2::new(4.0, 0.0),
    ];
    let c = Curve2::BSpline(cubic_bezier(pts));
    for k in 0..=20 {
        let t = k as Scalar / 20.0;
        let u = 1.0 - t;
        // B'(t) = 3 * sum (P[i+1]-P[i]) * Bernstein_2,i(t)
        let want = (pts[1] - pts[0]) * (3.0 * u * u)
            + (pts[2] - pts[1]) * (6.0 * u * t)
            + (pts[3] - pts[2]) * (3.0 * t * t);
        let got = derivative2(&c, t).unwrap();
        assert!(
            (got - want).length() < 1e-11,
            "hodograph at {t}: {got:?} vs {want:?}"
        );
    }
}

#[test]
fn a_rational_bspline_represents_a_circular_arc_exactly() {
    // The classic quarter-circle NURBS: degree 2, weights (1, 1/sqrt2, 1).
    // Any point must lie exactly on the unit circle -- this is the test that
    // separates true rational evaluation from averaging projected points.
    let w = 1.0 / Scalar::sqrt(2.0);
    let c = Curve2::BSpline(BSplineCurve2 {
        degree: 2,
        control_points: vec![
            Point2::new(1.0, 0.0),
            Point2::new(1.0, 1.0),
            Point2::new(0.0, 1.0),
        ],
        knots: vec![0.0, 1.0],
        multiplicities: vec![3, 3],
        weights: Some(vec![1.0, w, 1.0]),
        closed: false,
        self_intersect: None,
        knot_spec: KnotSpec::PiecewiseBezier,
    });
    for k in 0..=20 {
        let t = k as Scalar / 20.0;
        let p = evaluate2(&c, t).unwrap();
        assert!(
            (p.length() - 1.0).abs() < 1e-14,
            "rational arc left the unit circle at {t}: |p| = {}",
            p.length()
        );
    }
}

#[test]
fn a_multi_span_bspline_is_continuous_across_its_interior_knot() {
    let c = Curve2::BSpline(BSplineCurve2 {
        degree: 2,
        control_points: vec![
            Point2::new(0.0, 0.0),
            Point2::new(1.0, 2.0),
            Point2::new(2.0, 0.0),
            Point2::new(3.0, 2.0),
        ],
        knots: vec![0.0, 0.5, 1.0],
        multiplicities: vec![3, 1, 3],
        weights: None,
        closed: false,
        self_intersect: None,
        knot_spec: KnotSpec::Unspecified,
    });
    let eps = 1e-9;
    let before = evaluate2(&c, 0.5 - eps).unwrap();
    let at = evaluate2(&c, 0.5).unwrap();
    let after = evaluate2(&c, 0.5 + eps).unwrap();
    assert!(
        (at - before).length() < 1e-7,
        "discontinuous below the knot"
    );
    assert!((after - at).length() < 1e-7, "discontinuous above the knot");
}

#[test]
fn mismatched_compact_knots_are_refused_instead_of_truncated() {
    let mut spline = cubic_bezier([
        Point2::ZERO,
        Point2::new(1.0, 0.0),
        Point2::new(1.0, 1.0),
        Point2::new(2.0, 1.0),
    ]);
    // `zip` used to ignore this extra distinct knot and evaluate a different
    // curve whenever the truncated multiplicities happened to sum correctly.
    spline.knots.push(2.0);
    assert!(evaluate2(&Curve2::BSpline(spline), 0.5).is_err());
}

#[test]
fn hostile_multiplicities_are_bounded_before_expansion() {
    let mut spline = cubic_bezier([
        Point2::ZERO,
        Point2::new(1.0, 0.0),
        Point2::new(1.0, 1.0),
        Point2::new(2.0, 1.0),
    ]);
    spline.multiplicities[0] = u32::MAX;
    let curve = Curve2::BSpline(spline);
    assert_eq!(
        domain2(&curve),
        Interval {
            start: 0.0,
            end: 0.0
        }
    );
    assert!(evaluate2(&curve, 0.5).is_err());
}

#[test]
fn zero_knot_multiplicity_is_refused() {
    let spline = Curve2::BSpline(BSplineCurve2 {
        degree: 1,
        control_points: vec![Point2::ZERO, Point2::new(1.0, 0.0)],
        knots: vec![0.0, 0.5, 1.0],
        multiplicities: vec![2, 0, 2],
        weights: None,
        closed: false,
        self_intersect: None,
        knot_spec: KnotSpec::Unspecified,
    });
    assert!(evaluate2(&spline, 0.5).is_err());
}

#[test]
fn non_monotone_compact_knots_are_refused() {
    let spline = Curve2::BSpline(BSplineCurve2 {
        degree: 1,
        control_points: vec![Point2::ZERO, Point2::new(1.0, 0.0), Point2::new(2.0, 0.0)],
        knots: vec![0.0, 2.0, 1.0],
        multiplicities: vec![2, 1, 1],
        weights: None,
        closed: false,
        self_intersect: None,
        knot_spec: KnotSpec::Unspecified,
    });
    assert!(evaluate2(&spline, 0.5).is_err());
}

#[test]
fn non_positive_rational_weights_are_refused() {
    let mut spline = cubic_bezier([
        Point2::ZERO,
        Point2::new(1.0, 0.0),
        Point2::new(1.0, 1.0),
        Point2::new(2.0, 1.0),
    ]);
    spline.weights = Some(vec![1.0, 0.0, 1.0, 1.0]);
    assert!(evaluate2(&Curve2::BSpline(spline), 0.5).is_err());
}

#[test]
fn non_finite_bspline_control_points_are_refused() {
    let mut spline = cubic_bezier([
        Point2::ZERO,
        Point2::new(1.0, 0.0),
        Point2::new(1.0, 1.0),
        Point2::new(2.0, 1.0),
    ]);
    spline.control_points[1].x = Scalar::NAN;
    assert!(evaluate2(&Curve2::BSpline(spline), 0.5).is_err());
}

#[test]
fn a_degenerate_bspline_is_refused_not_guessed() {
    // Degree 3 with only 3 control points cannot define a span.
    let c = Curve2::BSpline(BSplineCurve2 {
        degree: 3,
        control_points: vec![Point2::ZERO, Point2::new(1.0, 0.0), Point2::new(1.0, 1.0)],
        knots: vec![0.0, 1.0],
        multiplicities: vec![4, 4],
        weights: None,
        closed: false,
        self_intersect: None,
        knot_spec: KnotSpec::PiecewiseBezier,
    });
    assert!(evaluate2(&c, 0.5).is_err(), "must refuse, not extrapolate");
}

#[test]
fn a_spline_with_an_empty_domain_is_refused_as_degenerate() {
    // Every count check passes: degree 2 with 3 control points needs 6 knot
    // entries, and 1 + 3 + 2 is 6 with no multiplicity above degree + 1. But
    // the expansion is [0, 1, 1, 1, 2, 2], so the domain runs from knot[2] = 1
    // to knot[3] = 1 and is empty. The test above no longer reaches this
    // guard, because the control-count check refuses it first (#296).
    let c = Curve2::BSpline(BSplineCurve2 {
        degree: 2,
        control_points: vec![Point2::ZERO, Point2::new(1.0, 0.0), Point2::new(1.0, 1.0)],
        knots: vec![0.0, 1.0, 2.0],
        multiplicities: vec![1, 3, 2],
        weights: None,
        closed: false,
        self_intersect: None,
        knot_spec: KnotSpec::Unspecified,
    });
    match evaluate2(&c, 1.0) {
        Err(GeomError::Degenerate(message)) => assert!(
            message.contains("domain is empty"),
            "refused for the wrong reason: {message}"
        ),
        other => panic!("an empty domain must be refused as degenerate, got {other:?}"),
    }
    // The infallible domain reports the refusal as the empty interval at 0,
    // not as the collapsed span [1, 1] of the expanded knots.
    assert_eq!(
        domain2(&c),
        Interval {
            start: 0.0,
            end: 0.0
        }
    );
}

// --- 3D ---------------------------------------------------------------------

#[test]
fn a_3d_circle_lies_in_its_frame_plane() {
    // A tilted frame: the circle must stay in the plane spanned by x and y.
    let normal = Vec3::new(1.0, 1.0, 1.0).normalize();
    let x = Vec3::new(1.0, -1.0, 0.0).normalize();
    let y = normal.cross(x);
    let c = Curve3::Circle(Circle3 {
        frame: Frame3 {
            origin: Point3::new(5.0, 0.0, -1.0),
            x,
            y,
            z: normal,
        },
        radius: 2.0,
    });
    for k in 0..32 {
        let t = TAU * k as Scalar / 32.0;
        let p = evaluate3(&c, t).unwrap();
        let offset = p - Point3::new(5.0, 0.0, -1.0);
        assert!(
            offset.dot(normal).abs() < 1e-13,
            "left the plane at {t} by {}",
            offset.dot(normal)
        );
        assert!((offset.length() - 2.0).abs() < 1e-13, "radius drift at {t}");
        // Tangent stays in-plane too.
        let d = derivative3(&c, t).unwrap();
        assert!(d.dot(normal).abs() < 1e-13, "tangent left the plane at {t}");
    }
}

// --- adaptive flattening ----------------------------------------------------

/// Largest distance from any flattened chord to the true curve.
fn max_sagitta(curve: &Curve2, pts: &[Point2], domain: Interval) -> Scalar {
    // Sample the true curve densely and measure to the nearest chord.
    let mut worst: Scalar = 0.0;
    let n = 2000;
    for k in 0..=n {
        let t = domain.start + (domain.end - domain.start) * k as Scalar / n as Scalar;
        let p = evaluate2(curve, t).unwrap();
        let mut best = Scalar::INFINITY;
        for w in pts.windows(2) {
            let (a, b) = (w[0], w[1]);
            let ab = b - a;
            let len2 = ab.length_squared();
            let d = if len2 <= 0.0 {
                (p - a).length()
            } else {
                let s = ((p - a).dot(ab) / len2).clamp(0.0, 1.0);
                (p - (a + ab * s)).length()
            };
            best = best.min(d);
        }
        worst = worst.max(best);
    }
    worst
}

#[test]
fn flattening_a_circle_respects_the_chord_tolerance() {
    let c = unit_circle();
    let domain = Interval {
        start: 0.0,
        end: TAU,
    };
    for tol in [0.1, 0.01, 0.001] {
        let pts = flatten2(&c, domain, tol, 24).unwrap();
        let err = max_sagitta(&c, &pts, domain);
        assert!(
            err <= tol * 1.05,
            "tolerance {tol} violated: max deviation {err} with {} points",
            pts.len()
        );
    }
}

#[test]
fn a_tighter_tolerance_produces_more_points() {
    let c = unit_circle();
    let domain = Interval {
        start: 0.0,
        end: TAU,
    };
    let coarse = flatten2(&c, domain, 0.1, 24).unwrap().len();
    let fine = flatten2(&c, domain, 0.001, 24).unwrap().len();
    assert!(
        fine > coarse,
        "tolerance must drive density: {coarse} vs {fine}"
    );
}

#[test]
fn flattening_a_line_emits_only_its_endpoints() {
    // Subdividing a straight line adds vertices that carry no information.
    let l = Curve2::Line(Line2 {
        origin: Point2::ZERO,
        direction: Vec2::new(1.0, 1.0),
    });
    let pts = flatten2(&l, Interval::UNIT, 1e-6, 24).unwrap();
    assert_eq!(pts.len(), 2, "a line needs no interior points");
}

#[test]
fn flattening_a_polyline_returns_its_own_vertices() {
    let p = Curve2::Polyline(Polyline2 {
        points: vec![
            Point2::new(0.0, 0.0),
            Point2::new(1.0, 0.0),
            Point2::new(1.0, 1.0),
        ],
        closed: false,
    });
    let pts = flatten2(
        &p,
        Interval {
            start: 0.0,
            end: 2.0,
        },
        1e-9,
        24,
    )
    .unwrap();
    assert_eq!(pts.len(), 3, "polyline breakpoints, no extras");
    assert_eq!(pts[1], Point2::new(1.0, 0.0));
}

#[test]
fn flattening_an_ellipse_respects_the_tolerance_where_curvature_peaks() {
    // A flat ellipse: curvature at the major-axis ends is ~a/b^2 = 100, so a
    // uniform segment count would fail there. Adaptive subdivision must not.
    let e = Curve2::Ellipse(Ellipse2 {
        frame: frame2(),
        semi_axis_x: 10.0,
        semi_axis_y: 1.0,
    });
    let domain = Interval {
        start: 0.0,
        end: TAU,
    };
    let tol = 0.01;
    let pts = flatten2(&e, domain, tol, 28).unwrap();
    let err = max_sagitta(&e, &pts, domain);
    assert!(err <= tol * 1.05, "high-curvature region violated: {err}");
}

#[test]
fn a_non_positive_chord_tolerance_is_refused() {
    let c = unit_circle();
    let domain = Interval {
        start: 0.0,
        end: TAU,
    };
    assert!(flatten2(&c, domain, 0.0, 24).is_err());
    assert!(flatten2(&c, domain, -1.0, 24).is_err());
    assert!(flatten2(&c, domain, Scalar::NAN, 24).is_err());
}

#[test]
fn a_non_finite_parameter_is_refused() {
    let c = unit_circle();
    assert!(evaluate2(&c, Scalar::NAN).is_err());
    assert!(evaluate2(&c, Scalar::INFINITY).is_err());
    assert!(derivative2(&c, Scalar::NAN).is_err());
}

// --- gaps found by mutation probing ----------------------------------------

/// A rational curve's derivative needs the quotient rule, not just A'/w.
///
/// The Bernstein hodograph test uses a *polynomial* spline, where w is
/// constant and w' is zero -- so it cannot see a missing quotient term. This
/// uses the rational quarter-arc, where the weight genuinely varies.
#[test]
fn a_rational_derivative_is_tangent_to_the_arc() {
    let w = 1.0 / Scalar::sqrt(2.0);
    let c = Curve2::BSpline(BSplineCurve2 {
        degree: 2,
        control_points: vec![
            Point2::new(1.0, 0.0),
            Point2::new(1.0, 1.0),
            Point2::new(0.0, 1.0),
        ],
        knots: vec![0.0, 1.0],
        multiplicities: vec![3, 3],
        weights: Some(vec![1.0, w, 1.0]),
        closed: false,
        self_intersect: None,
        knot_spec: KnotSpec::PiecewiseBezier,
    });
    for k in 0..=20 {
        let t = k as Scalar / 20.0;
        let p = evaluate2(&c, t).unwrap();
        let d = derivative2(&c, t).unwrap();
        // On a circle centred at the origin the tangent is perpendicular to
        // the radius. Dropping the quotient rule breaks this immediately.
        assert!(
            p.dot(d).abs() < 1e-9,
            "rational tangent not perpendicular at {t}: {}",
            p.dot(d)
        );
    }
}

/// The degree/control-point guard must fire, not merely exist.
///
/// The original test only checked degree 3 with 3 points. That case also fails
/// the knot-count check, so removing the control-point guard changed nothing.
#[test]
fn every_degenerate_spline_shape_is_refused() {
    // Degree 2, two control points, knot vector consistent with n + d + 1 = 5.
    let too_few = Curve2::BSpline(BSplineCurve2 {
        degree: 2,
        control_points: vec![Point2::ZERO, Point2::new(1.0, 0.0)],
        knots: vec![0.0, 1.0],
        multiplicities: vec![3, 2],
        weights: None,
        closed: false,
        self_intersect: None,
        knot_spec: KnotSpec::Unspecified,
    });
    assert!(
        evaluate2(&too_few, 0.5).is_err(),
        "degree 2 with 2 control points must be refused"
    );

    // Degree 0 is not a curve.
    let degree_zero = Curve2::BSpline(BSplineCurve2 {
        degree: 0,
        control_points: vec![Point2::ZERO, Point2::new(1.0, 0.0)],
        knots: vec![0.0, 1.0],
        multiplicities: vec![1, 2],
        weights: None,
        closed: false,
        self_intersect: None,
        knot_spec: KnotSpec::Unspecified,
    });
    assert!(
        evaluate2(&degree_zero, 0.5).is_err(),
        "degree 0 must be refused"
    );

    // Weight count must match control-point count.
    let bad_weights = Curve2::BSpline(BSplineCurve2 {
        degree: 2,
        control_points: vec![Point2::ZERO, Point2::new(1.0, 0.0), Point2::new(1.0, 1.0)],
        knots: vec![0.0, 1.0],
        multiplicities: vec![3, 3],
        weights: Some(vec![1.0, 1.0]),
        closed: false,
        self_intersect: None,
        knot_spec: KnotSpec::PiecewiseBezier,
    });
    assert!(
        evaluate2(&bad_weights, 0.5).is_err(),
        "mismatched weight count must be refused"
    );
}

/// Sagitta must be the distance to the chord, not to its midpoint.
///
/// An asymmetric curve deviates most away from the chord centre. Measuring
/// only to the midpoint under-reports there and the subdivision stops early.
#[test]
fn flattening_bounds_error_on_an_asymmetric_curve() {
    // A quarter ellipse: curvature varies sharply along the span, so the
    // maximum deviation is nowhere near the chord midpoint.
    let e = Curve2::Ellipse(Ellipse2 {
        frame: frame2(),
        semi_axis_x: 8.0,
        semi_axis_y: 1.0,
    });
    let domain = Interval {
        start: 0.0,
        end: TAU * 0.25,
    };
    let tol = 0.005;
    let pts = flatten2(&e, domain, tol, 28).unwrap();
    let err = max_sagitta(&e, &pts, domain);
    assert!(
        err <= tol * 1.05,
        "asymmetric span exceeded its budget: {err} > {tol}"
    );
}

/// A polyline domain that would silently discard vertices is refused.
#[test]
fn a_polyline_domain_that_drops_vertices_is_refused() {
    let p = Curve2::Polyline(Polyline2 {
        points: vec![
            Point2::new(0.0, 0.0),
            Point2::new(1.0, 0.0),
            Point2::new(1.0, 1.0),
            Point2::new(0.0, 1.0),
        ],
        closed: true,
    });
    // A closed 4-point ring spans (0, 4). Asking for (0, 1) would return just
    // the first edge -- data loss disguised as success.
    let narrow = flatten2(
        &p,
        Interval {
            start: 0.0,
            end: 1.0,
        },
        1e-6,
        24,
    );
    assert!(
        narrow.is_err(),
        "a domain covering 1 of 4 segments must be refused, got {:?}",
        narrow.map(|v| v.len())
    );

    // The true domain still works.
    let full = flatten2(
        &p,
        Interval {
            start: 0.0,
            end: 4.0,
        },
        1e-6,
        24,
    )
    .expect("the natural domain must be accepted");
    assert_eq!(full.len(), 5, "four corners plus the closing point");
}

/// Flattening must fail fast rather than exhaust memory.
///
/// Regression test for a real incident: a broken sagitta measure made the
/// tolerance test never succeed, and the depth-24 bound still permitted 2^24
/// segments. One test binary reached 189 CPU-minutes and 8 GB RSS before it
/// was killed. A depth bound is not a resource bound.
#[test]
fn flattening_refuses_rather_than_exhausting_memory() {
    // An impossible ask: a tolerance far below the representable resolution
    // of the coordinates, on a curve with real curvature.
    let c = Curve2::Circle(Circle2 {
        frame: frame2(),
        radius: 1e6,
    });
    let domain = Interval {
        start: 0.0,
        end: TAU,
    };
    // Deliberately absurd depth: the point budget, not the depth, must stop it.
    let result = flatten2(&c, domain, 1e-300, u32::MAX);
    assert!(
        result.is_err(),
        "an unsatisfiable tolerance must error, not allocate without bound"
    );
}

/// A satisfiable request near the budget still succeeds.
#[test]
fn a_demanding_but_satisfiable_tolerance_still_flattens() {
    let c = Curve2::Circle(Circle2 {
        frame: frame2(),
        radius: 1.0,
    });
    let domain = Interval {
        start: 0.0,
        end: TAU,
    };
    let pts = flatten2(&c, domain, 1e-7, 24).expect("a reachable tolerance must succeed");
    assert!(
        pts.len() > 100,
        "expected real subdivision, got {}",
        pts.len()
    );
    assert!(pts.len() < 65_536, "must stay inside the budget");
}

#[test]
fn flattening_fails_closed_when_depth_is_exhausted() {
    let curve = Curve2::Circle(Circle2 {
        frame: frame2(),
        radius: 1.0,
    });
    let error = flatten2(&curve, Interval::UNIT, 1.0e-12, 0)
        .expect_err("depth exhaustion must not return an unverified chord");
    assert!(matches!(error, GeomError::BudgetExceeded { .. }));
}

#[test]
fn flattening_fails_closed_when_the_parameter_midpoint_collapses() {
    let curve = Curve2::Circle(Circle2 {
        frame: frame2(),
        radius: 1.0,
    });
    let start = 2_f64.powi(53);
    let domain = Interval {
        start,
        end: start + 2.0,
    };
    let error = flatten2(&curve, domain, 1.0e-12, 20)
        .expect_err("an unbisectable curved interval must fail closed");
    assert!(matches!(error, GeomError::Degenerate(_)));
}

#[test]
fn rational_curve_overflow_is_refused() {
    let mut spline = cubic_bezier([Point2::splat(Scalar::MAX); 4]);
    spline.weights = Some(vec![Scalar::MAX; 4]);
    let curve = Curve2::BSpline(spline);
    assert!(evaluate2(&curve, 0.5).is_err());
    assert!(derivative2(&curve, 0.5).is_err());
}
