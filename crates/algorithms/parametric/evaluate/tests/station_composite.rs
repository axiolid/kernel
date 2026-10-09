//! Stations along a composite basis (#285): pieces of curves laid end to
//! end, joints read by the seam rule. Expected frames are closed forms and
//! stations on the individual pieces at the offset distance, never the
//! composite code under test.

use axiolid_core::{Frame2, Point2, Point3, Scalar, Transform3, Vec2, Vec3};
use axiolid_curve::{
    Circle2, Curve2, Curve3, Elevated3, ElevationLaw, Line2, Line3, Polyline2, SeamSide,
};
use axiolid_evaluate::station::{
    station_section2, CompositeBasis, DistanceConvention, SectionFrame, StationCurve, StationPiece,
};

const EPS: Scalar = 1e-12;
const PI: Scalar = core::f64::consts::PI;

fn close3(actual: Vec3, expected: Vec3, eps: Scalar, what: &str) {
    assert!(
        (actual - expected).abs().max_element() <= eps,
        "{what}: {actual:?} != {expected:?} (off by {:e})",
        (actual - expected).abs().max_element()
    );
}

fn same(a: &SectionFrame, b: &SectionFrame, eps: Scalar, what: &str) {
    close3(a.point, b.point, eps, &format!("{what}: point"));
    close3(a.tangent, b.tangent, eps, &format!("{what}: tangent"));
    close3(a.lateral, b.lateral, eps, &format!("{what}: lateral"));
    close3(a.up, b.up, eps, &format!("{what}: up"));
}

fn refused<T: core::fmt::Debug>(result: Result<T, axiolid_contracts::GeomError>, needle: &str) {
    match result {
        Err(error) => assert!(
            error.to_string().contains(needle),
            "expected a refusal naming {needle:?}, got {error}"
        ),
        Ok(value) => panic!("expected a refusal naming {needle:?}, got {value:?}"),
    }
}

/// An L: 3 along `+x`, then 4 along `+y`.
fn ell() -> Curve2 {
    Curve2::Polyline(Polyline2 {
        points: vec![Point2::ZERO, Point2::new(3.0, 0.0), Point2::new(3.0, 4.0)],
        closed: false,
    })
}

/// `+x` from the origin.
fn line() -> Curve2 {
    Curve2::Line(Line2 {
        origin: Point2::ZERO,
        direction: Vec2::X,
    })
}

/// Radius 5 about `(15, 0)`: at angle `pi` it passes `(10, 0)` heading
/// `-y`, so a line along `+x` ending there turns right by a quarter.
fn circle() -> Curve2 {
    Curve2::Circle(Circle2 {
        frame: Frame2 {
            origin: Point2::new(15.0, 0.0),
            x: Vec2::X,
            y: Vec2::Y,
        },
        radius: 5.0,
    })
}

/// The closed-form frame on that circle at angle `a`, counter-clockwise.
fn on_circle(a: Scalar) -> (Point3, Vec3, Vec3) {
    let (s, c) = a.sin_cos();
    let t = Vec3::new(-s, c, 0.0);
    (
        Point3::new(15.0 + 5.0 * c, 5.0 * s, 0.0),
        t,
        Vec3::new(-t.y, t.x, 0.0),
    )
}

/// The line over `[0, 10]`, then the quarter arc from angle `pi` to
/// `3 pi / 2`: a right-angle corner at distance 10, length `10 + 5 pi / 2`.
fn line_and_arc<'c>(line: &'c Curve2, circle: &'c Curve2) -> CompositeBasis<'c> {
    CompositeBasis::new(vec![
        StationPiece::between(StationCurve::Two(line), 0.0, 10.0).unwrap(),
        StationPiece::between_parameters(StationCurve::Two(circle), PI, 1.5 * PI).unwrap(),
    ])
    .unwrap()
}

#[test]
fn stations_on_a_line_and_an_arc_agree_with_the_pieces_at_the_offset_distance() {
    let (line, circle) = (line(), circle());
    let composite = line_and_arc(&line, &circle);
    assert!((composite.length() - (10.0 + 2.5 * PI)).abs() <= EPS);
    assert_eq!(composite.convention(), DistanceConvention::ArcLength3d);

    // Inside the line: the line's own station, and the closed form.
    for d in [0.0, 4.0, 9.5] {
        let on = composite.section_on(d, SeamSide::Outgoing).unwrap();
        same(&on, &station_section2(&line, d).unwrap(), EPS, "line piece");
        close3(on.point, Point3::new(d, 0.0, 0.0), EPS, "line point");
        close3(on.lateral, Vec3::Y, EPS, "line lateral");
    }
    // Inside the arc: the circle's station at r * pi + (d - 10), and the
    // closed form at angle pi + (d - 10) / r.
    for into in [0.5, 2.0, 7.0] {
        let d = 10.0 + into;
        let on = composite.section_on(d, SeamSide::Incoming).unwrap();
        same(
            &on,
            &station_section2(&circle, 5.0 * PI + into).unwrap(),
            1e-9,
            "arc piece",
        );
        let (p, t, l) = on_circle(PI + into / 5.0);
        close3(on.point, p, 1e-9, "arc point");
        close3(on.tangent, t, 1e-9, "arc tangent");
        close3(on.lateral, l, 1e-9, "arc lateral");
        close3(on.up, Vec3::Z, EPS, "arc up");
    }

    // On the joint, and within the tolerance on either side of it: the
    // outgoing side is the arc's start (heading -y, lateral +x), the
    // incoming one the line's end (heading +x, lateral +y).
    for d in [10.0, 10.0 - 1e-12, 10.0 + 1e-12] {
        let outgoing = composite.section_on(d, SeamSide::Outgoing).unwrap();
        let incoming = composite.section_on(d, SeamSide::Incoming).unwrap();
        close3(
            outgoing.point,
            Point3::new(10.0, 0.0, 0.0),
            1e-9,
            "joint point",
        );
        close3(
            incoming.point,
            Point3::new(10.0, 0.0, 0.0),
            EPS,
            "joint point",
        );
        close3(outgoing.tangent, -Vec3::Y, 1e-12, "outgoing tangent");
        close3(outgoing.lateral, Vec3::X, 1e-12, "outgoing lateral");
        close3(incoming.tangent, Vec3::X, EPS, "incoming tangent");
        close3(incoming.lateral, Vec3::Y, EPS, "incoming lateral");
        same(
            &incoming,
            &station_section2(&line, 10.0).unwrap(),
            EPS,
            "incoming is the line's end",
        );
    }
    // Just off the joint the side does not matter.
    let off = 10.0 + 1e-6;
    same(
        &composite.section_on(off, SeamSide::Incoming).unwrap(),
        &composite.section_on(off, SeamSide::Outgoing).unwrap(),
        EPS,
        "off the joint",
    );
    // At the start and the end one piece exists, and both sides read it.
    let end = composite.length();
    for d in [0.0, end] {
        same(
            &composite.section_on(d, SeamSide::Incoming).unwrap(),
            &composite.section_on(d, SeamSide::Outgoing).unwrap(),
            EPS,
            "an end",
        );
    }
    let (p, t, _) = on_circle(1.5 * PI);
    let last = composite.section_on(end, SeamSide::Incoming).unwrap();
    close3(last.point, p, 1e-9, "end point");
    close3(last.tangent, t, 1e-9, "end tangent");

    // Beyond the length, by name.
    refused(
        composite.section_on(end + 1e-3, SeamSide::Outgoing),
        "beyond the curve's length",
    );
}

#[test]
fn a_reversed_piece_is_read_from_its_end_with_tangent_and_lateral_negated() {
    // The arc from (10, 0) clockwise to (15, 5): the circle's span from
    // pi / 2 to pi, traversed backwards. Heading +y at the joint.
    let (line, circle) = (line(), circle());
    let composite = CompositeBasis::new(vec![
        StationPiece::between(StationCurve::Two(&line), 0.0, 10.0).unwrap(),
        StationPiece::between_parameters(StationCurve::Two(&circle), 0.5 * PI, PI)
            .unwrap()
            .reversed(),
    ])
    .unwrap();
    for into in [0.0, 1.0, 2.5 * PI] {
        let on = composite
            .section_on(10.0 + into, SeamSide::Outgoing)
            .unwrap();
        let (p, t, l) = on_circle(PI - into / 5.0);
        close3(on.point, p, 1e-9, "reversed point");
        close3(on.tangent, -t, 1e-9, "reversed tangent");
        close3(on.lateral, -l, 1e-9, "reversed lateral");
        close3(on.up, Vec3::Z, EPS, "reversed up");
    }
    // The joint's outgoing side is the reversed arc's start, heading +y.
    let joint = composite.section_on(10.0, SeamSide::Outgoing).unwrap();
    close3(joint.tangent, Vec3::Y, 1e-12, "turns left");
}

#[test]
fn a_gap_a_reversed_piece_and_mixed_conventions_are_refused_by_name() {
    let (line, circle) = (line(), circle());
    let short = StationPiece::between(StationCurve::Two(&line), 0.0, 9.5).unwrap();
    let arc = StationPiece::between_parameters(StationCurve::Two(&circle), PI, 1.5 * PI).unwrap();
    refused(CompositeBasis::new(vec![short, arc]), "a gap of 0.5");
    // The arc from (15, 5) to (10, 0), not declared reversed: its END meets
    // the line.
    let backwards =
        StationPiece::between_parameters(StationCurve::Two(&circle), 0.5 * PI, PI).unwrap();
    let full = StationPiece::between(StationCurve::Two(&line), 0.0, 10.0).unwrap();
    refused(
        CompositeBasis::new(vec![full, backwards]),
        "a reversed piece",
    );
    // The line itself reversed against an arc that starts at its far end.
    refused(
        CompositeBasis::new(vec![full.reversed(), arc]),
        "a reversed piece",
    );
    refused(CompositeBasis::new(Vec::new()), "no pieces");

    // An elevated piece measures plan distance, a 3D line arc length.
    let elevated = Curve3::Elevated(Elevated3::new(
        Curve2::Line(Line2 {
            origin: Point2::new(10.0, 0.0),
            direction: Vec2::X,
        }),
        ElevationLaw::constant_grade(0.0, 0.05),
    ));
    let flat = Curve3::Line(Line3 {
        origin: Point3::ZERO,
        direction: Vec3::X,
    });
    refused(
        CompositeBasis::new(vec![
            StationPiece::between(StationCurve::Three(&flat), 0.0, 10.0).unwrap(),
            StationPiece::between(StationCurve::Three(&elevated), 0.0, 5.0).unwrap(),
        ]),
        "measure stations differently",
    );
}

#[test]
fn a_composite_of_elevated_pieces_is_measured_in_plan_distance() {
    // Two grades end to end on a straight plan: 5% over [0, 10], then -2%
    // from (10, 0, 0.5). A station at 15 is plan distance 15: x = 15, not
    // the 3D arc length.
    let rise = Curve3::Elevated(Elevated3::new(
        line(),
        ElevationLaw::constant_grade(0.0, 0.05),
    ));
    let fall = Curve3::Elevated(Elevated3::new(
        Curve2::Line(Line2 {
            origin: Point2::new(10.0, 0.0),
            direction: Vec2::X,
        }),
        ElevationLaw::constant_grade(0.5, -0.02),
    ));
    let composite = CompositeBasis::new(vec![
        StationPiece::between(StationCurve::Three(&rise), 0.0, 10.0).unwrap(),
        StationPiece::between(StationCurve::Three(&fall), 0.0, 20.0).unwrap(),
    ])
    .unwrap();
    assert_eq!(composite.convention(), DistanceConvention::PlanDistance);
    assert!((composite.length() - 30.0).abs() <= EPS);
    let on = composite.section_on(15.0, SeamSide::Outgoing).unwrap();
    close3(
        on.point,
        Point3::new(15.0, 0.0, 0.4),
        EPS,
        "plan distance 15",
    );
    // The grade break at the joint: incoming 5%, outgoing -2%.
    let k = |g: Scalar| (1.0 + g * g).sqrt();
    let incoming = composite.section_on(10.0, SeamSide::Incoming).unwrap();
    let outgoing = composite.section_on(10.0, SeamSide::Outgoing).unwrap();
    close3(
        incoming.tangent,
        Vec3::new(1.0, 0.0, 0.05) / k(0.05),
        EPS,
        "incoming grade",
    );
    close3(
        outgoing.tangent,
        Vec3::new(1.0, 0.0, -0.02) / k(0.02),
        EPS,
        "outgoing grade",
    );
}

#[test]
fn seams_are_the_joints_and_the_pieces_own_and_exactness_follows_the_pieces() {
    // 2 along +x to the origin, a polyline (a corner at its 3), then the
    // line on from its end: joints at 2 and 9, the vertex at 5.
    let lead = Curve2::Line(Line2 {
        origin: Point2::new(-2.0, 0.0),
        direction: Vec2::X,
    });
    let ell = ell();
    let up = Curve2::Line(Line2 {
        origin: Point2::new(3.0, 4.0),
        direction: Vec2::new(0.0, 2.0),
    });
    let composite = CompositeBasis::new(vec![
        StationPiece::between(StationCurve::Two(&lead), 0.0, 2.0).unwrap(),
        StationPiece::whole(StationCurve::Two(&ell)).unwrap(),
        StationPiece::whole(StationCurve::Two(&up)).unwrap(),
    ])
    .unwrap();
    // An untrimmed line is its domain [0, 1]: length 2.
    assert!((composite.length() - 11.0).abs() <= EPS);
    let seams = composite.exact_seams().unwrap();
    let at: Vec<Scalar> = seams.iter().map(|s| s.distance).collect();
    assert_eq!(at, vec![2.0, 5.0, 9.0]);
    assert!(seams.iter().all(|s| s.exact && !s.smooth));
    // The polyline's vertex inside the composite reads its sides.
    let incoming = composite.section_on(5.0, SeamSide::Incoming).unwrap();
    close3(incoming.tangent, Vec3::X, EPS, "vertex incoming");
    let outgoing = composite.section_on(5.0, SeamSide::Outgoing).unwrap();
    close3(outgoing.tangent, Vec3::Y, EPS, "vertex outgoing");
    close3(
        outgoing.point,
        Point3::new(3.0, 0.0, 0.0),
        EPS,
        "vertex point",
    );

    // Exactness: lines only are exact; after an arc nothing is.
    let (line, circle) = (line(), circle());
    let mixed = line_and_arc(&line, &circle);
    assert!(mixed.frame_is_exact_at(4.0, SeamSide::Outgoing));
    assert!(mixed.frame_is_exact_at(10.0, SeamSide::Incoming));
    assert!(!mixed.frame_is_exact_at(10.0, SeamSide::Outgoing));
    assert!(!mixed.frame_is_exact_at(12.0, SeamSide::Outgoing));
    assert!(composite.frame_is_exact_at(1.0, SeamSide::Outgoing));
    assert!(!composite.frame_is_exact_at(3.0, SeamSide::Outgoing));
    // A line after an arc is not exact: the arc's length places it.
    let tail = Curve2::Line(Line2 {
        origin: Point2::new(15.0, -5.0),
        direction: Vec2::X,
    });
    let arc_first = CompositeBasis::new(vec![
        StationPiece::between_parameters(StationCurve::Two(&circle), PI, 1.5 * PI).unwrap(),
        StationPiece::between(StationCurve::Two(&tail), 0.0, 5.0).unwrap(),
    ])
    .unwrap();
    assert!(!arc_first.frame_is_exact_at(2.5 * PI + 1.0, SeamSide::Outgoing));
    // An inexactly placed line is not exact either.
    let placed = StationPiece::between(StationCurve::Two(&line), 0.0, 1.0)
        .unwrap()
        .placed(
            Transform3::from_translation(Vec3::new(0.0, 0.0, 1.0)),
            false,
        );
    assert!(!placed.frame_is_exact());
    assert!(StationPiece::between(StationCurve::Two(&line), 0.0, 1.0)
        .unwrap()
        .frame_is_exact());

    // A trim of the composite keeps the pieces between two distances.
    let clipped = mixed.pieces_between(5.0, 12.0).unwrap();
    let trimmed = CompositeBasis::new(clipped).unwrap();
    assert!((trimmed.length() - 7.0).abs() <= 1e-12);
    same(
        &trimmed.section_on(6.0, SeamSide::Outgoing).unwrap(),
        &mixed.section_on(11.0, SeamSide::Outgoing).unwrap(),
        1e-12,
        "trimmed station",
    );
}

#[test]
fn a_placed_piece_carries_its_frame_or_reads_its_own_where_tilted() {
    let line = line();
    let piece = StationPiece::between(StationCurve::Two(&line), 0.0, 4.0).unwrap();
    // Upright: a quarter turn about +Z and a lift, carried as is.
    let turn = Transform3::from_mat3_translation(
        axiolid_core::Mat3::from_rotation_z(0.5 * PI),
        Vec3::new(1.0, 2.0, 3.0),
    );
    let upright = piece.placed(turn, true);
    let on = upright.section_on(1.0, SeamSide::Outgoing).unwrap();
    close3(on.point, Point3::new(1.0, 3.0, 3.0), EPS, "carried point");
    close3(on.tangent, Vec3::Y, EPS, "carried tangent");
    close3(on.lateral, -Vec3::X, EPS, "carried lateral");
    // Tilted about x by 30 degrees: the line still runs +x, so its own
    // reference-up frame is upright, not the tilted carried one.
    let tilt = Transform3::from_rotation_x(PI / 6.0);
    let on = piece
        .placed(tilt, true)
        .section_on(1.0, SeamSide::Outgoing)
        .unwrap();
    close3(on.tangent, Vec3::X, EPS, "tilted tangent");
    close3(on.lateral, Vec3::Y, EPS, "own lateral");
    close3(on.up, Vec3::Z, EPS, "own up");
    // A plan-measured curve so placed is refused by name.
    let grade = Curve3::Elevated(Elevated3::new(
        line,
        ElevationLaw::constant_grade(0.0, 0.01),
    ));
    refused(
        StationPiece::between(StationCurve::Three(&grade), 0.0, 4.0)
            .unwrap()
            .placed(tilt, true)
            .section_on(1.0, SeamSide::Outgoing),
        "tilts +Z",
    );
}

#[test]
fn a_piece_is_read_from_inside_at_its_ends_and_reversed_at_its_seams() {
    // The L from its vertex on: at the start, both sides read +y, never
    // the trimmed-away first segment.
    let ell = ell();
    let from_vertex = CompositeBasis::new(vec![StationPiece::between(
        StationCurve::Two(&ell),
        3.0,
        7.0,
    )
    .unwrap()])
    .unwrap();
    for side in [SeamSide::Incoming, SeamSide::Outgoing] {
        let start = from_vertex.section_on(0.0, side).unwrap();
        close3(start.tangent, Vec3::Y, EPS, "trimmed at the vertex");
        close3(start.point, Point3::new(3.0, 0.0, 0.0), EPS, "vertex point");
    }
    // The L up to its vertex: at the end, both sides read +x, never the
    // segment trimmed away after it.
    let to_vertex = CompositeBasis::new(vec![StationPiece::between(
        StationCurve::Two(&ell),
        0.0,
        3.0,
    )
    .unwrap()])
    .unwrap();
    for side in [SeamSide::Incoming, SeamSide::Outgoing] {
        let end = to_vertex.section_on(3.0, side).unwrap();
        close3(end.tangent, Vec3::X, EPS, "trimmed to the vertex");
    }
    // The whole L backwards: from (3, 4) down to the vertex at 4, then to
    // the origin. Incoming at the vertex is heading -y, outgoing -x.
    let backwards = CompositeBasis::new(vec![StationPiece::whole(StationCurve::Two(&ell))
        .unwrap()
        .reversed()])
    .unwrap();
    let incoming = backwards.section_on(4.0, SeamSide::Incoming).unwrap();
    let outgoing = backwards.section_on(4.0, SeamSide::Outgoing).unwrap();
    close3(incoming.tangent, -Vec3::Y, EPS, "reversed incoming");
    close3(outgoing.tangent, -Vec3::X, EPS, "reversed outgoing");
    close3(outgoing.lateral, -Vec3::Y, EPS, "reversed lateral");
    let seams = backwards.exact_seams().unwrap();
    assert_eq!(seams.len(), 1);
    assert!((seams[0].distance - 4.0).abs() <= EPS);
}

#[test]
fn a_trim_of_a_reversed_piece_keeps_its_span() {
    let (line, circle) = (line(), circle());
    let composite = CompositeBasis::new(vec![
        StationPiece::between(StationCurve::Two(&line), 0.0, 10.0).unwrap(),
        StationPiece::between_parameters(StationCurve::Two(&circle), 0.5 * PI, PI)
            .unwrap()
            .reversed(),
    ])
    .unwrap();
    let trimmed = CompositeBasis::new(composite.pieces_between(5.0, 12.0).unwrap()).unwrap();
    for d in [1.0, 6.0, 7.0] {
        same(
            &trimmed.section_on(d, SeamSide::Outgoing).unwrap(),
            &composite.section_on(d + 5.0, SeamSide::Outgoing).unwrap(),
            1e-12,
            "trimmed reversed piece",
        );
    }
}

#[test]
fn a_circle_trimmed_across_its_parameter_seam_is_one_piece() {
    // Radius 2 about the origin from 3 pi / 2 to 5 pi / 2: up the right
    // side from (0, -2) to (0, 2), through angle 0 at distance pi.
    let circle = Curve2::Circle(Circle2 {
        frame: Frame2 {
            origin: Point2::ZERO,
            x: Vec2::X,
            y: Vec2::Y,
        },
        radius: 2.0,
    });
    let composite = CompositeBasis::new(vec![StationPiece::between_parameters(
        StationCurve::Two(&circle),
        1.5 * PI,
        2.5 * PI,
    )
    .unwrap()])
    .unwrap();
    assert!((composite.length() - 2.0 * PI).abs() <= 1e-12);
    for (d, angle) in [
        (0.5, 1.5 * PI + 0.25),
        (PI, 2.0 * PI),
        (1.5 * PI, 2.25 * PI),
    ] {
        let on = composite.section_on(d, SeamSide::Outgoing).unwrap();
        let (s, c) = Scalar::sin_cos(angle);
        close3(on.point, Point3::new(2.0 * c, 2.0 * s, 0.0), 1e-9, "point");
        close3(on.tangent, Vec3::new(-s, c, 0.0), 1e-9, "tangent");
    }
}

#[test]
fn a_joint_after_an_ellipse_is_not_exact() {
    // A quarter ellipse from (2, 0) to (0, 1), then on along -x.
    let ellipse = Curve2::Ellipse(axiolid_curve::Ellipse2 {
        frame: Frame2 {
            origin: Point2::ZERO,
            x: Vec2::X,
            y: Vec2::Y,
        },
        semi_axis_x: 2.0,
        semi_axis_y: 1.0,
    });
    let back = Curve2::Line(Line2 {
        origin: Point2::new(0.0, 1.0),
        direction: -Vec2::X,
    });
    let composite = CompositeBasis::new(vec![
        StationPiece::between_parameters(StationCurve::Two(&ellipse), 0.0, 0.5 * PI).unwrap(),
        StationPiece::between(StationCurve::Two(&back), 0.0, 1.0).unwrap(),
    ])
    .unwrap();
    let seams = composite.seams().unwrap();
    assert_eq!(seams.len(), 1);
    assert!(!seams[0].exact);
    refused(composite.exact_seams(), "arc length is a quadrature");
}
