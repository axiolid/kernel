//! The signed bend of a section frame law (#361, ADR 0082 amendment):
//! `FrameLaw::bend` bounds `O'.N`, `L.N'` and `U.N'` over a span.
//!
//! The references are independent of the bend's closed-form jets: the
//! forward rate from the geometry (arc length, or `1 / |T_xy|` in plan
//! distance), the two turning coefficients from central differences of the
//! frames' own tangent, which every family reads in closed form. Each
//! sampled value must lie inside the reported range, and the ranges must
//! be tight enough to tell the sign of the bend.

use axiolid_core::{Frame2, Point2, Point3, Scalar, Transform3, Vec2, Vec3};
use axiolid_curve::{
    BSplineCurve, BankConvention, Banked3, CantLaw, CantPiece, Chain2, ChainPiece2, Circle2,
    CurvatureLaw, Curve2, Curve3, Elevated3, ElevationLaw, Ellipse2, Intrinsic2, KnotSpec, Line2,
    Polyline2,
};
use axiolid_evaluate::station::{
    station_section2, station_section3, BendBounds, CompositeBasis, FrameLaw, FrameReading,
    SectionFrame, StationCurve, StationPiece,
};

const TOLERANCE: Scalar = 1e-7;

fn contains(range: [Scalar; 2], value: Scalar) -> bool {
    range[0] - 1e-7 <= value && value <= range[1] + 1e-7
}

/// The coefficients at `s` from `frame`: `L.T'` and `U.T'` by central
/// differences of the tangent in distance (the frames' axes are closed
/// form, so these are accurate to about `1e-9`).
fn turning(frame: &dyn Fn(Scalar) -> SectionFrame, s: Scalar) -> (Scalar, Scalar) {
    let h = 1e-5;
    let rate = (frame(s + h).tangent - frame(s - h).tangent) / (2.0 * h);
    let at = frame(s);
    (at.lateral.dot(rate), at.up.dot(rate))
}

/// Each coefficient sampled over `[lo, hi]` lies within the bend's range.
fn check(
    bend: &BendBounds,
    frame: &dyn Fn(Scalar) -> SectionFrame,
    forward: &dyn Fn(Scalar) -> Scalar,
    (lo, hi): (Scalar, Scalar),
    plan: bool,
) {
    for i in 0..=40 {
        let s = lo + (hi - lo) * Scalar::from(i) / 40.0;
        let (lateral, up) = turning(frame, s);
        assert!(contains(bend.forward, forward(s)), "forward at {s}");
        assert!(contains(bend.lateral, lateral), "lateral {lateral} at {s}");
        let up = if plan { 0.0 } else { up };
        assert!(contains(bend.up, up), "up {up} at {s}");
    }
    // The middle value is within its own error of the coefficient.
    let mid = 0.5 * (lo + hi);
    let (lateral, up) = turning(frame, mid);
    assert!((bend.middle[0] - forward(mid)).abs() <= bend.middle_error[0] + 1e-7);
    assert!((bend.middle[1] - lateral).abs() <= bend.middle_error[1] + 1e-7);
    if !plan {
        assert!((bend.middle[2] - up).abs() <= bend.middle_error[2] + 1e-7);
    }
}

fn circle(radius: Scalar) -> Curve2 {
    Curve2::Circle(Circle2 {
        frame: Frame2 {
            origin: Point2::new(1.0, 2.0),
            x: Vec2::X,
            y: Vec2::Y,
        },
        radius,
    })
}

#[test]
fn a_counter_clockwise_circle_turns_toward_its_lateral_axis_at_its_curvature() {
    let curve = circle(10.0);
    let law = FrameLaw::atomic(StationCurve::Two(&curve), None, TOLERANCE).unwrap();
    for reading in [FrameReading::Section, FrameReading::Plan] {
        let span = (2.0, 6.0);
        let bend = law.bend(span.0, span.1, reading).unwrap();
        // Tight enough to tell the sign of the bend.
        assert!(contains(bend.lateral, 0.1) && bend.lateral[0] > 0.0);
        assert!(contains(bend.forward, 1.0) && bend.forward[0] > 0.9);
        assert!(contains(bend.up, 0.0));
        check(
            &bend,
            &|s| station_section2(&curve, s).unwrap(),
            &|_| 1.0,
            span,
            false,
        );
    }
}

#[test]
fn a_reversed_piece_turns_the_other_way() {
    let curve = circle(10.0);
    let piece = StationPiece::between(StationCurve::Two(&curve), 0.0, 20.0)
        .unwrap()
        .reversed();
    let basis = CompositeBasis::new(vec![piece]).unwrap();
    let law = FrameLaw::composite(&basis, TOLERANCE).unwrap();
    let span = (3.0, 8.0);
    let bend = law.bend(span.0, span.1, FrameReading::Section).unwrap();
    assert!(contains(bend.lateral, -0.1) && bend.lateral[1] < 0.0);
    assert!(contains(bend.forward, 1.0));
    let reversed = |s: Scalar| {
        piece
            .section_on(s, axiolid_curve::SeamSide::Outgoing)
            .unwrap()
    };
    check(&bend, &reversed, &|_| 1.0, span, false);
}

/// `z = s^2 / 100` over a line along `x`.
fn sag() -> Curve3 {
    Curve3::Elevated(Elevated3::new(
        Curve2::Line(Line2 {
            origin: Point2::ZERO,
            direction: Vec2::X,
        }),
        ElevationLaw::Polynomial {
            coefficients: vec![0.0, 0.0, 0.01],
        },
    ))
}

#[test]
fn a_gradient_curve_bends_up_by_its_vertical_curvature_in_the_section_frame_only() {
    let curve = sag();
    let law = FrameLaw::atomic(StationCurve::Three(&curve), None, TOLERANCE).unwrap();
    let span = (5.0, 15.0);
    let section = law.bend(span.0, span.1, FrameReading::Section).unwrap();
    let grade = |s: Scalar| 0.02 * s;
    let speed = |s: Scalar| (1.0 + grade(s) * grade(s)).sqrt();
    for i in 0..=20 {
        let s = span.0 + (span.1 - span.0) * Scalar::from(i) / 20.0;
        assert!(contains(section.forward, speed(s)), "forward at {s}");
        assert!(
            contains(section.up, 0.02 / (speed(s) * speed(s))),
            "up at {s}"
        );
        assert!(contains(section.lateral, 0.0), "lateral at {s}");
    }
    // Narrow enough to tell the sign of the bend.
    let narrow = law.bend(9.0, 11.0, FrameReading::Section).unwrap();
    assert!(narrow.up[0] > 0.0);
    check(
        &section,
        &|s| station_section3(&curve, s).unwrap(),
        &|s| 1.0 / station_section3(&curve, s).unwrap().tangent.x.hypot(0.0),
        span,
        false,
    );
    let plan = law.bend(span.0, span.1, FrameReading::Plan).unwrap();
    assert_eq!(plan.up, [0.0, 0.0]);
    assert!(contains(plan.forward, 1.0) && plan.forward[0] > 0.9);
    assert!(contains(plan.lateral, 0.0));
}

/// A clothoid from a straight start whose curvature grows as `s / (R L)`.
fn clothoid() -> Curve2 {
    Curve2::Intrinsic(Intrinsic2::new(
        Frame2 {
            origin: Point2::new(3.0, -1.0),
            x: Vec2::new(0.6, 0.8),
            y: Vec2::new(-0.8, 0.6),
        },
        CurvatureLaw::clothoid(0.0, 1.0 / 30.0, 60.0),
        60.0,
    ))
}

#[test]
fn a_clothoid_turns_at_its_curvature() {
    let curve = clothoid();
    let law = FrameLaw::atomic(StationCurve::Two(&curve), None, TOLERANCE).unwrap();
    let span = (10.0, 20.0);
    let bend = law.bend(span.0, span.1, FrameReading::Section).unwrap();
    for i in 0..=20 {
        let s = span.0 + (span.1 - span.0) * Scalar::from(i) / 20.0;
        assert!(contains(bend.lateral, s / (30.0 * 60.0)), "at {s}");
    }
    assert!(contains(bend.forward, 1.0));
    check(
        &bend,
        &|s| station_section2(&curve, s).unwrap(),
        &|_| 1.0,
        span,
        false,
    );
}

#[test]
fn a_gradient_curve_over_a_clothoid_bends_in_the_plan_and_the_profile() {
    let curve = Curve3::Elevated(Elevated3::new(
        clothoid(),
        ElevationLaw::circular_arc(10.0, 0.01, 800.0),
    ));
    let law = FrameLaw::atomic(StationCurve::Three(&curve), None, TOLERANCE).unwrap();
    let span = (10.0, 20.0);
    // The profile arc's grade at plan distance `d`, in closed form.
    let grade = |d: Scalar| {
        let t0 = 0.01_f64.atan();
        (t0.sin() + d / 800.0).asin().tan()
    };
    for reading in [FrameReading::Section, FrameReading::Plan] {
        let bend = law.bend(span.0, span.1, reading).unwrap();
        let plan = reading == FrameReading::Plan;
        let frame = |s: Scalar| {
            let frame = station_section3(&curve, s).unwrap();
            if plan {
                frame.plan().unwrap()
            } else {
                frame
            }
        };
        let forward = |s: Scalar| {
            if plan {
                1.0
            } else {
                (1.0 + grade(s) * grade(s)).sqrt()
            }
        };
        check(&bend, &frame, &forward, span, plan);
    }
}

#[test]
fn frames_read_by_quadrature_have_a_bend_and_chains_none() {
    let ellipse = Curve2::Ellipse(Ellipse2 {
        frame: Frame2 {
            origin: Point2::new(30.0, 0.0),
            x: Vec2::X,
            y: Vec2::Y,
        },
        semi_axis_x: 10.0,
        semi_axis_y: 6.0,
    });
    let spline = Curve3::BSpline(BSplineCurve {
        degree: 3,
        control_points: vec![
            Point3::new(0.0, 0.0, 0.0),
            Point3::new(8.0, 2.0, 0.5),
            Point3::new(16.0, -3.0, 1.5),
            Point3::new(24.0, 2.0, 2.0),
            Point3::new(30.0, 6.0, 2.2),
        ],
        knots: vec![0.0, 0.5, 1.0],
        multiplicities: vec![4, 1, 4],
        weights: None,
        closed: false,
        self_intersect: None,
        knot_spec: KnotSpec::Unspecified,
    });
    let chain = Curve2::Chain(Chain2::new(
        Frame2 {
            origin: Point2::ZERO,
            x: Vec2::X,
            y: Vec2::Y,
        },
        vec![ChainPiece2::Intrinsic {
            curvature: CurvatureLaw::clothoid(0.0, 1.0 / 30.0, 60.0),
            length: 60.0,
        }],
    ));
    for (name, curve, bent) in [
        ("ellipse", StationCurve::Two(&ellipse), true),
        ("b-spline", StationCurve::Three(&spline), true),
        ("chain", StationCurve::Two(&chain), false),
    ] {
        let law = FrameLaw::atomic(curve, None, TOLERANCE).unwrap();
        for reading in [FrameReading::Section, FrameReading::Plan] {
            let bend = law.bend(5.0, 8.0, reading);
            assert_eq!(bend.is_some(), bent, "{name}");
            if let Some(bend) = bend {
                assert!(contains(bend.forward, 1.0), "{name}");
            }
        }
        assert_eq!(law.has_bend(5.0, 8.0), bent, "{name}");
    }
}

#[test]
fn a_polyline_segment_is_straight_and_a_vertex_has_no_bend() {
    let curve = Curve2::Polyline(Polyline2 {
        points: vec![Point2::ZERO, Point2::new(5.0, 0.0), Point2::new(5.0, 5.0)],
        closed: false,
    });
    let law = FrameLaw::atomic(StationCurve::Two(&curve), None, TOLERANCE).unwrap();
    for reading in [FrameReading::Section, FrameReading::Plan] {
        let bend = law.bend(1.0, 3.0, reading).unwrap();
        assert!(contains(bend.forward, 1.0) && contains(bend.lateral, 0.0));
        assert!(contains(bend.up, 0.0));
        let bend = law.bend(6.0, 9.0, reading).unwrap();
        assert!(contains(bend.forward, 1.0) && contains(bend.lateral, 0.0));
        assert!(law.bend(4.0, 6.0, reading).is_none());
    }
}

#[test]
fn a_banked_curve_turns_into_its_rolled_axes() {
    let plan = Curve2::Circle(Circle2 {
        frame: Frame2 {
            origin: Point2::ZERO,
            x: Vec2::X,
            y: Vec2::Y,
        },
        radius: 100.0,
    });
    let curve = Curve3::Banked(Banked3::new(
        Elevated3::new(plan, ElevationLaw::level(0.0)),
        CantLaw::new(vec![CantPiece::constant(60.0, 0.1)]),
        CantLaw::zero(60.0),
        1.5,
        BankConvention::TangentRotation,
    ));
    let law = FrameLaw::atomic(StationCurve::Three(&curve), None, TOLERANCE).unwrap();
    let span = (10.0, 20.0);
    for reading in [FrameReading::Section, FrameReading::Plan] {
        let plan = reading == FrameReading::Plan;
        let bend = law.bend(span.0, span.1, reading).unwrap();
        let frame = |s: Scalar| {
            let frame = station_section3(&curve, s).unwrap();
            if plan {
                frame.plan().unwrap()
            } else {
                frame
            }
        };
        check(&bend, &frame, &|_| 1.0, span, plan);
    }
}

#[test]
fn a_placed_piece_scales_the_forward_rate_and_keeps_the_turning() {
    let curve = circle(10.0);
    let squeeze = Transform3::from_cols(
        Vec3::new(0.0, 2.0, 0.0),
        Vec3::new(-2.0, 0.0, 0.0),
        Vec3::new(0.0, 0.0, 2.0),
        Point3::new(5.0, 0.0, 0.0),
    );
    let piece = StationPiece::between(StationCurve::Two(&curve), 0.0, 20.0)
        .unwrap()
        .placed(squeeze, true);
    let law = FrameLaw::atomic(StationCurve::Two(&curve), Some(squeeze), TOLERANCE).unwrap();
    let span = (3.0, 8.0);
    let bend = law.bend(span.0, span.1, FrameReading::Section).unwrap();
    assert!(contains(bend.forward, 2.0) && bend.forward[0] > 1.5);
    assert!(contains(bend.lateral, 0.1));
    let placed = |s: Scalar| {
        piece
            .section_on(s, axiolid_curve::SeamSide::Outgoing)
            .unwrap()
    };
    check(&bend, &placed, &|_| 2.0, span, false);
}

#[test]
fn the_plan_reading_of_a_tilted_space_circle_advances_at_its_horizontal_part() {
    let (x, y) = (Vec3::X, Vec3::new(0.0, 0.8, 0.6));
    let curve = Curve3::Circle(axiolid_curve::Circle3 {
        frame: axiolid_core::Frame3 {
            origin: Point3::ZERO,
            x,
            y,
            z: x.cross(y),
        },
        radius: 5.0,
    });
    let law = FrameLaw::atomic(StationCurve::Three(&curve), None, TOLERANCE).unwrap();
    let span = (0.5, 2.5);
    let bend = law.bend(span.0, span.1, FrameReading::Plan).unwrap();
    let frame = |s: Scalar| station_section3(&curve, s).unwrap().plan().unwrap();
    let horizontal = |s: Scalar| {
        let t = station_section3(&curve, s).unwrap().tangent;
        t.x.hypot(t.y)
    };
    check(&bend, &frame, &horizontal, span, true);
    assert!(bend.middle[0] < 1.0);
}

#[test]
fn a_span_too_long_for_a_bound_still_has_a_closed_form_bend() {
    let curve = Curve3::Circle(axiolid_curve::Circle3 {
        frame: axiolid_core::Frame3 {
            origin: Point3::ZERO,
            x: Vec3::X,
            y: Vec3::Y,
            z: Vec3::Z,
        },
        radius: 10.0,
    });
    let law = FrameLaw::atomic(StationCurve::Three(&curve), None, TOLERANCE).unwrap();
    // Longer than twice the radius: no bound over it, but halving gives one.
    assert!(law.bend(0.0, 30.0, FrameReading::Plan).is_none());
    assert!(law.has_bend(0.0, 30.0));
    assert!(law.bend(0.0, 7.5, FrameReading::Plan).is_some());
    let ellipse = Curve2::Ellipse(Ellipse2 {
        frame: Frame2 {
            origin: Point2::ZERO,
            x: Vec2::X,
            y: Vec2::Y,
        },
        semi_axis_x: 10.0,
        semi_axis_y: 6.0,
    });
    let law = FrameLaw::atomic(StationCurve::Two(&ellipse), None, TOLERANCE).unwrap();
    assert!(law.has_bend(1.0, 2.0));
}
