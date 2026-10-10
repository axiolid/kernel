//! Offset pieces of a station basis (#289): stations along an offset are
//! measured in the offset's own length. Expected values are closed forms
//! (a line beside a line, a circle of radius `r - d` beside a circle) or an
//! independent evaluation of the offset curve written out here (its
//! analytic point and speed, integrated by Simpson's rule), never the
//! offset code under test.

use axiolid_contracts::GeomResult;
use axiolid_core::{Frame2, Point2, Scalar, Vec2, Vec3};
use axiolid_curve::{
    BankConvention, Banked3, CantLaw, CantPiece, Circle2, CurvatureLaw, Curve2, Curve3, CurvePath,
    Elevated3, ElevationLaw, Intrinsic2, Line2, OffsetFrame, OffsetLaw, PathCurve, PathOffsets,
    Polyline2, SeamSide,
};
use axiolid_curve_evaluate_contract::{CurveEvaluator, CurveMeasure, DistanceConvention};
use axiolid_evaluate::station::{
    offset_pieces, CompositeBasis, StationCurve, StationOffset, StationPiece,
};
use axiolid_evaluate::ReferenceCurveEvaluator;

const PI: Scalar = core::f64::consts::PI;
const SIDES: [SeamSide; 2] = [SeamSide::Incoming, SeamSide::Outgoing];

fn close(actual: Vec3, expected: Vec3, eps: Scalar, what: &str) {
    let off = (actual - expected).length();
    assert!(
        off <= eps,
        "{what}: {actual:?} != {expected:?} (off by {off:e})"
    );
}

fn refused<T: core::fmt::Debug>(result: GeomResult<T>, needle: &str) {
    match result {
        Err(error) => assert!(
            error.to_string().contains(needle),
            "expected a refusal naming {needle:?}, got {error}"
        ),
        Ok(value) => panic!("expected a refusal naming {needle:?}, got {value:?}"),
    }
}

fn line() -> Curve2 {
    Curve2::Line(Line2 {
        origin: Point2::ZERO,
        direction: Vec2::X,
    })
}

fn circle(centre: Point2, radius: Scalar) -> Curve2 {
    Curve2::Circle(Circle2 {
        frame: Frame2 {
            origin: centre,
            x: Vec2::X,
            y: Vec2::Y,
        },
        radius,
    })
}

fn planar(distance: Scalar) -> OffsetLaw {
    OffsetLaw::Planar { distance }
}

fn linear(start: PathOffsets, end: PathOffsets) -> OffsetLaw {
    OffsetLaw::Linear {
        start,
        end,
        frame: OffsetFrame::Section,
    }
}

fn offset<'c>(base: &[StationPiece<'c>], law: OffsetLaw) -> CompositeBasis<'c> {
    CompositeBasis::new(offset_pieces(base, law).unwrap()).unwrap()
}

#[test]
fn a_constant_offset_of_a_line_is_a_line_beside_it_read_exactly() {
    let line = line();
    let base = StationPiece::between(StationCurve::Two(&line), 0.0, 10.0).unwrap();
    let basis = offset(&[base], planar(1.5));
    assert_eq!(basis.length(), 10.0);
    assert_eq!(basis.convention(), DistanceConvention::ArcLength3d);
    for s in [0.0, 3.3, 10.0] {
        let section = basis.section_on(s, SeamSide::Outgoing).unwrap();
        close(section.point, Vec3::new(s, 1.5, 0.0), 1e-15, "point");
        assert_eq!(
            (section.tangent, section.lateral, section.up),
            (Vec3::X, Vec3::Y, Vec3::Z)
        );
        assert!(basis.frame_is_exact_at(s, SeamSide::Outgoing));
    }
    // To the right: negative.
    let right = offset(&[base], planar(-2.0));
    close(
        right.section_on(4.0, SeamSide::Outgoing).unwrap().point,
        Vec3::new(4.0, -2.0, 0.0),
        1e-15,
        "right",
    );
    // A 3D offset along +Z x T is the same left offset, also exact.
    let directed = offset(
        &[base],
        OffsetLaw::Directed {
            distance: 1.5,
            reference_direction: Vec3::Z,
        },
    );
    close(
        directed.section_on(3.3, SeamSide::Outgoing).unwrap().point,
        Vec3::new(3.3, 1.5, 0.0),
        1e-15,
        "directed",
    );
    assert!(directed.frame_is_exact_at(3.3, SeamSide::Outgoing));
}

#[test]
fn a_constant_offset_of_an_arc_is_an_arc_of_radius_r_minus_d() {
    // Radius 5 about the origin, a quarter from (5, 0) heading +y.
    let circle = circle(Point2::ZERO, 5.0);
    let arc = StationPiece::between_parameters(StationCurve::Two(&circle), 0.0, 0.5 * PI).unwrap();
    for (piece, distance, radius, sense) in [
        // Left of a counter-clockwise arc is towards the centre.
        (arc, 1.0, 4.0, 1.0),
        (arc, -1.0, 6.0, 1.0),
        // Reversed, it runs clockwise and the left is outwards.
        (arc.reversed(), 1.0, 6.0, -1.0),
    ] {
        let basis = offset(&[piece], planar(distance));
        let length = radius * 0.5 * PI;
        assert!(
            (basis.length() - length).abs() <= 1e-12,
            "length {} != {length}",
            basis.length()
        );
        for s in [0.0, 0.3 * length, length] {
            let theta = if sense > 0.0 {
                s / radius
            } else {
                0.5 * PI - s / radius
            };
            let section = basis.section_on(s, SeamSide::Outgoing).unwrap();
            let (sin, cos) = theta.sin_cos();
            close(
                section.point,
                Vec3::new(radius * cos, radius * sin, 0.0),
                1e-11,
                "point",
            );
            close(
                section.tangent,
                sense * Vec3::new(-sin, cos, 0.0),
                1e-11,
                "tangent",
            );
            close(section.up, Vec3::Z, 0.0, "up");
            // A circle's frame is read through the arc-length inverse.
            assert!(!basis.frame_is_exact_at(s, SeamSide::Outgoing));
        }
        // The seams after it are exact: its length is a closed form.
        assert!(basis.exact_seams().is_ok());
    }
}

#[test]
fn a_constant_offset_of_a_line_and_arc_agrees_on_both_sides_of_the_joint() {
    // (0, 0) to (10, 0), then a quarter of radius 5 about (10, 5) turning
    // left, tangent at the joint.
    let line = line();
    let circle = circle(Point2::new(10.0, 5.0), 5.0);
    let base = [
        StationPiece::between(StationCurve::Two(&line), 0.0, 10.0).unwrap(),
        StationPiece::between_parameters(StationCurve::Two(&circle), 1.5 * PI, 2.0 * PI).unwrap(),
    ];
    let basis = offset(&base, planar(1.0));
    assert!((basis.length() - (10.0 + 2.0 * PI)).abs() <= 1e-12);
    let seams = basis.exact_seams().unwrap();
    assert_eq!(seams.len(), 1);
    assert!((seams[0].distance - 10.0).abs() <= 1e-12 && !seams[0].smooth);
    for side in SIDES {
        let joint = basis.section_on(10.0, side).unwrap();
        close(joint.point, Vec3::new(10.0, 1.0, 0.0), 1e-11, "joint point");
        close(joint.tangent, Vec3::X, 1e-11, "joint tangent");
        // Within the seam tolerance either side, the side names the piece.
        let hair = basis.section_on(10.0 + 0.4e-11, side).unwrap();
        close(hair.point, joint.point, 1e-11, "hair after");
    }
    assert!(basis.frame_is_exact_at(10.0, SeamSide::Incoming));
    assert!(!basis.frame_is_exact_at(10.0, SeamSide::Outgoing));
    for t in [0.5, 2.0 * PI] {
        let theta = 1.5 * PI + t / 4.0;
        let section = basis.section_on(10.0 + t, SeamSide::Outgoing).unwrap();
        close(
            section.point,
            Vec3::new(10.0 + 4.0 * theta.cos(), 5.0 + 4.0 * theta.sin(), 0.0),
            1e-11,
            "on the arc",
        );
    }
}

#[test]
fn an_offset_by_distances_along_a_line_is_the_chord_of_its_offsets() {
    let line = line();
    let base = StationPiece::between(StationCurve::Two(&line), 0.0, 10.0).unwrap();
    let basis = offset(
        &[base],
        linear(
            PathOffsets::new(1.0, 0.0, 0.0),
            PathOffsets::new(3.0, 0.5, 0.0),
        ),
    );
    // From (0, 1, 0) to (10, 3, 0.5): its own length.
    let chord = Vec3::new(10.0, 2.0, 0.5);
    let length = chord.length();
    assert!(
        (basis.length() - length).abs() <= 1e-9 * length,
        "length {} != {length}",
        basis.length()
    );
    for s in [0.0, 0.25 * length, 0.5 * length, length] {
        let section = basis.section_on(s, SeamSide::Outgoing).unwrap();
        close(
            section.point,
            Vec3::new(0.0, 1.0, 0.0) + chord * (s / length),
            1e-8,
            "point",
        );
        close(section.tangent, chord / length, 1e-8, "tangent");
        assert!(!basis.frame_is_exact_at(s, SeamSide::Outgoing));
    }
}

/// The offset of a counter-clockwise arc of radius `r` about the origin,
/// from angle 0, whose lateral offset grows linearly from 0 to `w` over
/// `span` radians: radius `rho(theta) = r - w theta / span`. Its own arc
/// length to `theta` is `int sqrt(rho^2 + rho'^2)`, by Simpson's rule.
fn spiral_length(r: Scalar, w: Scalar, span: Scalar, theta: Scalar) -> Scalar {
    let speed = |t: Scalar| {
        let rho = r - w * t / span;
        (rho * rho + (w / span) * (w / span)).sqrt()
    };
    let n = 20_000;
    let h = theta / Scalar::from(n);
    let mut sum = speed(0.0) + speed(theta);
    for k in 1..n {
        sum += speed(h * Scalar::from(k)) * if k % 2 == 1 { 4.0 } else { 2.0 };
    }
    sum * h / 3.0
}

#[test]
fn an_offset_by_distances_along_an_arc_agrees_with_an_independent_evaluation() {
    let (r, w, span) = (10.0, 2.0, 0.5 * PI);
    let circle = circle(Point2::ZERO, r);
    let arc = StationPiece::between_parameters(StationCurve::Two(&circle), 0.0, span).unwrap();
    let basis = offset(
        &[arc],
        linear(PathOffsets::default(), PathOffsets::new(w, 0.0, 0.0)),
    );
    let total = spiral_length(r, w, span, span);
    assert!(
        (basis.length() - total).abs() <= 1e-9 * total,
        "length {} != {total}",
        basis.length()
    );
    for fraction in [0.2, 0.5, 0.9] {
        // The angle whose own length is `s`, by bisection.
        let s = fraction * total;
        let (mut lo, mut hi) = (0.0, span);
        for _ in 0..60 {
            let mid = 0.5 * (lo + hi);
            if spiral_length(r, w, span, mid) < s {
                lo = mid;
            } else {
                hi = mid;
            }
        }
        let theta = 0.5 * (lo + hi);
        let rho = r - w * theta / span;
        let (sin, cos) = theta.sin_cos();
        let section = basis.section_on(s, SeamSide::Outgoing).unwrap();
        close(
            section.point,
            Vec3::new(rho * cos, rho * sin, 0.0),
            1e-7,
            "point",
        );
        let velocity = Vec3::new(
            -rho * sin - (w / span) * cos,
            rho * cos - (w / span) * sin,
            0.0,
        );
        close(section.tangent, velocity.normalize(), 1e-7, "tangent");
        close(section.up, Vec3::Z, 0.0, "up");
    }
}

#[test]
fn an_offset_of_an_elevated_curve_is_measured_in_its_own_plan_length() {
    // A plan circle of radius 10 on a 5% grade from height 2; 1 m to the
    // left (towards the centre) is a plan circle of radius 9, at the
    // base's height for the same angle.
    let elevated = Curve3::Elevated(Elevated3::new(
        circle(Point2::ZERO, 10.0),
        ElevationLaw::constant_grade(2.0, 0.05),
    ));
    let quarter = 5.0 * PI;
    let base = StationPiece::between(StationCurve::Three(&elevated), 0.0, quarter).unwrap();
    let basis = offset(
        &[base],
        linear(
            PathOffsets::new(1.0, 0.0, 0.0),
            PathOffsets::new(1.0, 0.0, 0.0),
        ),
    );
    assert_eq!(basis.convention(), DistanceConvention::PlanDistance);
    let plan = 9.0 * 0.5 * PI;
    assert!(
        (basis.length() - plan).abs() <= 1e-9 * plan,
        "plan length {} != {plan}",
        basis.length()
    );
    for fraction in [0.0, 0.4, 1.0] {
        let s = fraction * plan;
        let theta = s / 9.0;
        let section = basis.section_on(s, SeamSide::Outgoing).unwrap();
        let (sin, cos) = theta.sin_cos();
        close(
            section.point,
            Vec3::new(9.0 * cos, 9.0 * sin, 2.0 + 0.05 * 10.0 * theta),
            1e-8,
            "point",
        );
        // Its own grade is steeper: the base's rise over a shorter plan.
        let heading = Vec3::new(-sin, cos, 0.05 * 10.0 / 9.0).normalize();
        close(section.tangent, heading, 1e-8, "tangent");
        assert!(section.lateral.z.abs() <= 1e-12, "lateral stays level");
    }
    // A plan-measured offset does not join an arc-length-measured piece.
    let line = line();
    let after = StationPiece::between(StationCurve::Two(&line), 0.0, 1.0).unwrap();
    let mut pieces = offset_pieces(&[base], planar_section(1.0)).unwrap();
    pieces.push(after);
    refused(
        CompositeBasis::new(pieces),
        "no one distance runs through it",
    );
}

#[test]
fn an_offset_reads_its_frame_upright_in_plan_and_rolled_beside_a_banked_curve() {
    // A 10% grade on a straight plan: one up in the plan frame is +Z, in the
    // section frame the leaning up.
    let grade = 0.1;
    let elevated = Curve3::Elevated(Elevated3::new(
        line(),
        ElevationLaw::constant_grade(2.0, grade),
    ));
    let base = StationPiece::between(StationCurve::Three(&elevated), 0.0, 20.0).unwrap();
    let up = PathOffsets::new(0.0, 1.0, 0.0);
    let lean = (1.0 + grade * grade).sqrt();
    for (frame, lift) in [
        (OffsetFrame::Plan, Vec3::Z),
        (OffsetFrame::Section, Vec3::new(-grade, 0.0, 1.0) / lean),
    ] {
        let basis = offset(
            &[base],
            OffsetLaw::Linear {
                start: up,
                end: up,
                frame,
            },
        );
        assert!((basis.length() - 20.0).abs() <= 1e-9 * 20.0, "{frame:?}");
        let section = basis.section_on(5.0, SeamSide::Outgoing).unwrap();
        // A constant lift moves the offset rigidly: its own plan distance 5
        // is the base's 5, the point lifted (and, in the section frame,
        // leaning back).
        let s = 5.0;
        close(
            section.point,
            Vec3::new(s, 0.0, 2.0 + grade * s) + lift,
            1e-8,
            &format!("{frame:?}"),
        );
    }
    // Beside a banked curve the offset keeps the rolled lateral.
    let banked = Curve3::Banked(Banked3::new(
        Elevated3::new(line(), ElevationLaw::level(0.0)),
        CantLaw::new(vec![CantPiece::constant(100.0, 0.15)]),
        CantLaw::zero(100.0),
        1.5,
        BankConvention::TangentRotation,
    ));
    let rail = StationPiece::between(StationCurve::Three(&banked), 0.0, 100.0).unwrap();
    let basis = offset(&[rail], planar_section(0.75));
    assert_eq!(basis.convention(), DistanceConvention::PlanDistance);
    let own = basis.section_on(40.0, SeamSide::Outgoing).unwrap();
    let at = rail.section_on(40.0, SeamSide::Outgoing).unwrap();
    assert!(at.lateral.z.abs() > 0.05, "the base is rolled");
    close(own.point, at.point + 0.75 * at.lateral, 1e-9, "point");
    close(own.lateral, at.lateral, 1e-9, "the rolled lateral");
    close(own.up, at.up, 1e-9, "the rolled up");
}

fn planar_section(lateral: Scalar) -> OffsetLaw {
    linear(
        PathOffsets::new(lateral, 0.0, 0.0),
        PathOffsets::new(lateral, 0.0, 0.0),
    )
}

#[test]
fn a_collapse_a_cusp_and_a_self_crossing_are_refused_by_name() {
    let circle = circle(Point2::ZERO, 5.0);
    let arc = StationPiece::between_parameters(StationCurve::Two(&circle), 0.0, 0.5 * PI).unwrap();
    // To the centre, and beyond it.
    refused(offset_pieces(&[arc], planar(5.0)), "collapses");
    refused(offset_pieces(&[arc], planar(7.0)), "collapses");
    // A widening that reaches the centre part way: a cusp.
    refused(
        offset_pieces(
            &[arc],
            linear(PathOffsets::default(), PathOffsets::new(6.0, 0.0, 0.0)),
        ),
        "cusp",
    );
    // A longitudinal law running back faster than the base runs ahead.
    let line = line();
    let straight = StationPiece::between(StationCurve::Two(&line), 0.0, 10.0).unwrap();
    refused(
        offset_pieces(
            &[straight],
            linear(PathOffsets::default(), PathOffsets::new(0.0, 0.0, -20.0)),
        ),
        "cusp",
    );
    // A tightening spiral offset outwards by more than its turns' spacing
    // crosses its own earlier turn.
    let spiral = Curve2::Intrinsic(Intrinsic2::new(
        Frame2 {
            origin: Point2::ZERO,
            x: Vec2::X,
            y: Vec2::Y,
        },
        CurvatureLaw::clothoid(0.0, 1.0, 30.0),
        30.0,
    ));
    let whole = StationPiece::whole(StationCurve::Two(&spiral)).unwrap();
    refused(offset_pieces(&[whole], planar(-3.0)), "crosses itself");
}

#[test]
fn an_offset_across_a_corner_is_refused_and_one_through_it_is_a_seam() {
    let corner = Curve2::Polyline(Polyline2 {
        points: vec![
            Point2::ZERO,
            Point2::new(10.0, 0.0),
            Point2::new(10.0, 10.0),
        ],
        closed: false,
    });
    let whole = StationPiece::whole(StationCurve::Two(&corner)).unwrap();
    // One offset piece per span between the polyline's vertices.
    let pieces = offset_pieces(&[whole], planar(1.0)).unwrap();
    assert_eq!(pieces.len(), 2);
    refused(CompositeBasis::new(pieces), "an offset across a corner");
    // An offset that vanishes at the corner passes through it: a seam.
    let through = offset(
        &[
            StationPiece::between(StationCurve::Two(&corner), 0.0, 10.0).unwrap(),
            StationPiece::between(StationCurve::Two(&corner), 10.0, 20.0).unwrap(),
        ],
        linear(
            PathOffsets::new(0.0, 0.0, 0.0),
            PathOffsets::new(0.0, 2.0, 0.0),
        ),
    );
    let seams = through.seams().unwrap();
    assert_eq!(seams.len(), 1);
    let incoming = through
        .section_on(seams[0].distance, SeamSide::Incoming)
        .unwrap();
    let outgoing = through
        .section_on(seams[0].distance, SeamSide::Outgoing)
        .unwrap();
    close(incoming.point, outgoing.point, 1e-9, "one point");
    assert!(incoming.tangent.x > 0.99 && outgoing.tangent.y > 0.99);
    // A base span crossing a seam of its curve is refused: split it.
    refused(
        StationOffset::new(whole, planar(1.0)),
        "crosses a seam of its curve",
    );
}

#[test]
fn a_curve_path_with_offset_pieces_reads_as_the_composite_bitwise() {
    let line = line();
    let circle = circle(Point2::new(10.0, 5.0), 5.0);
    let base = [
        StationPiece::between(StationCurve::Two(&line), 0.0, 10.0).unwrap(),
        StationPiece::between_parameters(StationCurve::Two(&circle), 1.5 * PI, 2.0 * PI).unwrap(),
    ];
    let e = ReferenceCurveEvaluator::new();
    for law in [
        planar(1.0),
        linear(
            PathOffsets::new(1.0, 0.0, 0.0),
            PathOffsets::new(2.0, 0.5, 0.0),
        ),
    ] {
        let path: CurvePath = offset_pieces(&base, law)
            .unwrap()
            .into_iter()
            .map(Into::into)
            .collect();
        assert!(path
            .pieces()
            .iter()
            .all(|piece| matches!(piece.curve, PathCurve::Offset(_))));
        let composite = CompositeBasis::from_path(&path).unwrap();
        assert_eq!(composite.path(), path, "the path round-trips");
        assert_eq!(
            e.path_distance_convention(&path),
            DistanceConvention::ArcLength3d
        );
        let length = composite.length();
        let joint = composite.seams().unwrap()[0].distance;
        for d in [0.0, 3.0, joint, joint + 1.0, length] {
            let at = CurveMeasure::Distance(d);
            for side in SIDES {
                let station = composite.section_on(d, side).unwrap();
                assert_eq!(e.path_point_at_on(&path, at, side).unwrap(), station.point);
                assert_eq!(
                    e.path_tangent_at_on(&path, at, side).unwrap(),
                    station.tangent
                );
                assert_eq!(
                    e.path_frame_at_on(&path, at, side).unwrap(),
                    station.frame()
                );
                assert_eq!(
                    e.path_frame_is_exact_at(&path, at, side),
                    composite.frame_is_exact_at(d, side)
                );
            }
        }
    }
}
