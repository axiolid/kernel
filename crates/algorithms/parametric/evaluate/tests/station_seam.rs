//! Stations on seams (#263): which piece a station on a seam reads, where
//! the seams are, and the mitre between two pieces. Expected frames come
//! from closed forms and from plain station evaluation just off the seam,
//! never from the seam code under test.

use axiolid_contracts::GeomError;
use axiolid_core::{Frame2, Point2, Point3, Scalar, Vec2, Vec3};
use axiolid_curve::{
    BSplineCurve, BankConvention, Banked3, CantLaw, CantPiece, Chain2, ChainPiece2, CurvatureLaw,
    Curve2, Curve3, Elevated3, ElevationLaw, KnotSpec, Line2, Polyline2, Polyline3, SeamSide,
};
use axiolid_evaluate::station::{
    exact_station_seams2, exact_station_seams3, station_seams2, station_section2,
    station_section2_on, station_section3, station_section3_on, Mitre, SectionFrame,
    MITRE_TOLERANCE, SEAM_TANGENT_TOLERANCE,
};

const EPS: Scalar = 1e-12;

fn close3(actual: Vec3, expected: Vec3, eps: Scalar, what: &str) {
    assert!(
        (actual - expected).abs().max_element() <= eps,
        "{what}: {actual:?} != {expected:?} (off by {:e})",
        (actual - expected).abs().max_element()
    );
}

fn same_axes(a: &SectionFrame, b: &SectionFrame, eps: Scalar, what: &str) {
    close3(a.tangent, b.tangent, eps, &format!("{what}: tangent"));
    close3(a.lateral, b.lateral, eps, &format!("{what}: lateral"));
    close3(a.up, b.up, eps, &format!("{what}: up"));
}

fn refused(error: GeomError, needle: &str) {
    assert!(
        error.to_string().contains(needle),
        "expected a refusal naming {needle:?}, got {error}"
    );
}

fn ell2() -> Curve2 {
    Curve2::Polyline(Polyline2 {
        points: vec![
            Point2::new(0.0, 0.0),
            Point2::new(3.0, 0.0),
            Point2::new(3.0, 4.0),
        ],
        closed: false,
    })
}

#[test]
fn a_polyline_vertex_reads_the_side_it_is_asked_for() {
    let curve = ell2();
    let corner = Point3::new(3.0, 0.0, 0.0);
    // On the vertex, and within the arc-length tolerance on either side.
    for s in [3.0, 3.0 - 1e-13, 3.0 + 1e-13] {
        let outgoing = station_section2_on(&curve, s, SeamSide::Outgoing).unwrap();
        let incoming = station_section2_on(&curve, s, SeamSide::Incoming).unwrap();
        close3(
            outgoing.point,
            corner,
            0.0,
            "outgoing reads the vertex itself",
        );
        close3(
            incoming.point,
            corner,
            0.0,
            "incoming reads the vertex itself",
        );
        close3(outgoing.tangent, Vec3::Y, EPS, "outgoing: the second leg");
        close3(outgoing.lateral, -Vec3::X, EPS, "outgoing lateral");
        close3(incoming.tangent, Vec3::X, EPS, "incoming: the first leg");
        close3(incoming.lateral, Vec3::Y, EPS, "incoming lateral");
        // The default is the outgoing piece.
        assert_eq!(station_section2(&curve, s).unwrap(), outgoing);
    }
    // Off the seam the side does not matter.
    for s in [1.0, 2.9, 3.1, 7.0, 0.0] {
        assert_eq!(
            station_section2_on(&curve, s, SeamSide::Incoming).unwrap(),
            station_section2(&curve, s).unwrap(),
            "at {s}"
        );
    }
}

#[test]
fn each_side_of_a_seam_is_the_limit_of_its_own_piece() {
    // Independent evaluation: plain stations just before and after the
    // seam, read by the ordinary path, frame the two sides.
    let h = 1e-6;
    let three = Curve3::Polyline(Polyline3 {
        points: vec![
            Point3::new(0.0, 0.0, 0.0),
            Point3::new(4.0, 0.0, 3.0),
            Point3::new(4.0, 6.0, 3.0),
        ],
        closed: false,
    });
    let elevated = Curve3::Elevated(grade_break());
    for (curve, seam, what) in [
        (&three, 5.0, "3D polyline"),
        (&elevated, 50.0, "grade break"),
    ] {
        let before = station_section3(curve, seam - h).unwrap();
        let after = station_section3(curve, seam + h).unwrap();
        let incoming = station_section3_on(curve, seam, SeamSide::Incoming).unwrap();
        let outgoing = station_section3_on(curve, seam, SeamSide::Outgoing).unwrap();
        same_axes(&incoming, &before, 1e-9, &format!("{what}: incoming"));
        same_axes(&outgoing, &after, 1e-9, &format!("{what}: outgoing"));
        close3(
            incoming.point,
            outgoing.point,
            1e-12,
            &format!("{what}: one point"),
        );
        assert!(
            incoming.tangent.dot(outgoing.tangent) < 1.0 - 1e-6,
            "{what}: the seam turns"
        );
    }
}

fn grade_break() -> Elevated3 {
    Elevated3::new(
        Curve2::Line(Line2 {
            origin: Point2::ZERO,
            direction: Vec2::X,
        }),
        ElevationLaw::Piecewise {
            breaks: vec![50.0],
            laws: vec![
                ElevationLaw::constant_grade(0.0, 0.02),
                ElevationLaw::constant_grade(1.0, -0.01),
            ],
        },
    )
}

#[test]
fn a_grade_break_reads_the_grade_of_the_selected_piece() {
    let curve = Curve3::Elevated(grade_break());
    let frame = |grade: Scalar| {
        let k = grade.hypot(1.0);
        (
            Vec3::new(1.0 / k, 0.0, grade / k),
            Vec3::new(-grade / k, 0.0, 1.0 / k),
        )
    };
    for (side, grade) in [(SeamSide::Incoming, 0.02), (SeamSide::Outgoing, -0.01)] {
        let section = station_section3_on(&curve, 50.0, side).unwrap();
        let (tangent, up) = frame(grade);
        close3(section.point, Point3::new(50.0, 0.0, 1.0), EPS, "point");
        close3(section.tangent, tangent, EPS, "tangent");
        close3(section.lateral, Vec3::Y, EPS, "lateral");
        close3(section.up, up, EPS, "up leans with the selected grade");
    }
}

#[test]
fn a_cant_jump_reads_the_roll_of_the_selected_piece() {
    // Level and straight: 150 mm cant for 100 m, then none. Rotation about
    // the tangent, b = 1.5: the incoming lateral rises by 0.15 / 1.5.
    let curve = Curve3::Banked(Banked3::new(
        Elevated3::new(
            Curve2::Line(Line2 {
                origin: Point2::ZERO,
                direction: Vec2::X,
            }),
            ElevationLaw::level(0.0),
        ),
        CantLaw::new(vec![
            CantPiece::constant(100.0, 0.15),
            CantPiece::constant(50.0, 0.0),
        ]),
        CantLaw::zero(150.0),
        1.5,
        BankConvention::TangentRotation,
    ));
    let incoming = station_section3_on(&curve, 100.0, SeamSide::Incoming).unwrap();
    let outgoing = station_section3_on(&curve, 100.0, SeamSide::Outgoing).unwrap();
    close3(incoming.tangent, Vec3::X, EPS, "tangent");
    close3(outgoing.tangent, Vec3::X, EPS, "tangent");
    close3(
        incoming.lateral,
        Vec3::new(0.0, 0.99_f64.sqrt(), 0.1),
        EPS,
        "rolled by the 150 mm cant",
    );
    close3(outgoing.lateral, Vec3::Y, EPS, "level after the seam");
    assert_eq!(station_section3(&curve, 100.0).unwrap(), outgoing);
    // The seam is listed, not smooth: the roll jumps although the tangent
    // does not, so no mitre is cut.
    let seams = exact_station_seams3(&curve).unwrap();
    assert_eq!(seams.len(), 1);
    assert_eq!(seams[0].distance, 100.0);
    assert!(!seams[0].smooth && seams[0].exact);
    assert!(Mitre::between(&incoming, &outgoing).unwrap().is_none());
}

/// Degree 2 with a double interior knot: the curve passes through the
/// middle control point with a corner there, `+x` arriving and `+y`
/// leaving. The control points before it are collinear, so the arc length
/// to the corner is exactly 2.
fn cornered_spline() -> Curve2 {
    Curve2::BSpline(BSplineCurve {
        degree: 2,
        control_points: vec![
            Point2::new(0.0, 0.0),
            Point2::new(1.0, 0.0),
            Point2::new(2.0, 0.0),
            Point2::new(2.0, 1.0),
            Point2::new(2.0, 2.0),
        ],
        knots: vec![0.0, 1.0, 2.0],
        multiplicities: vec![3, 2, 3],
        weights: None,
        closed: false,
        self_intersect: None,
        knot_spec: KnotSpec::Unspecified,
    })
}

#[test]
fn a_spline_corner_knot_reads_either_span() {
    let curve = cornered_spline();
    let seams = station_seams2(&curve).unwrap();
    assert_eq!(seams.len(), 1);
    assert_eq!(seams[0].parameter, 1.0);
    assert!((seams[0].distance - 2.0).abs() <= 1e-12);
    assert!(!seams[0].exact, "a quadrature, not stored data");
    refused(
        exact_station_seams2(&curve).unwrap_err(),
        "its arc length to the knot is a quadrature",
    );
    let incoming = station_section2_on(&curve, 2.0, SeamSide::Incoming).unwrap();
    let outgoing = station_section2_on(&curve, 2.0, SeamSide::Outgoing).unwrap();
    close3(incoming.point, Point3::new(2.0, 0.0, 0.0), EPS, "corner");
    close3(outgoing.point, Point3::new(2.0, 0.0, 0.0), EPS, "corner");
    close3(incoming.tangent, Vec3::X, EPS, "arriving along +x");
    close3(outgoing.tangent, Vec3::Y, EPS, "leaving along +y");
    // And each is the limit of its own span.
    let h = 1e-7;
    same_axes(
        &incoming,
        &station_section2(&curve, 2.0 - h).unwrap(),
        1e-6,
        "incoming limit",
    );
    same_axes(
        &outgoing,
        &station_section2(&curve, 2.0 + h).unwrap(),
        1e-6,
        "outgoing limit",
    );
}

#[test]
fn seams_are_read_from_stored_data() {
    // A polyline: interior vertices at the running sum of the legs, a
    // repeated vertex one seam.
    let polyline = Curve2::Polyline(Polyline2 {
        points: vec![
            Point2::new(0.0, 0.0),
            Point2::new(3.0, 0.0),
            Point2::new(3.0, 0.0),
            Point2::new(3.0, 4.0),
            Point2::new(0.0, 4.0),
        ],
        closed: false,
    });
    let seams = exact_station_seams2(&polyline).unwrap();
    let at: Vec<(Scalar, Scalar)> = seams.iter().map(|s| (s.distance, s.parameter)).collect();
    assert_eq!(at, vec![(3.0, 2.0), (7.0, 3.0)]);
    assert!(seams.iter().all(|s| !s.smooth && s.exact));
    let incoming = station_section2_on(&polyline, 3.0, SeamSide::Incoming).unwrap();
    close3(incoming.tangent, Vec3::X, EPS, "past the repeated vertex");
    // A chain: its joins and its pieces' curvature seams, smooth.
    let chain = Curve2::Chain(Chain2::new(
        Frame2 {
            origin: Point2::ZERO,
            x: Vec2::X,
            y: Vec2::Y,
        },
        vec![
            ChainPiece2::Intrinsic {
                curvature: CurvatureLaw::circular(0.0),
                length: 10.0,
            },
            ChainPiece2::Intrinsic {
                curvature: CurvatureLaw::Piecewise {
                    breaks: vec![5.0],
                    laws: vec![CurvatureLaw::circular(0.01), CurvatureLaw::circular(0.02)],
                },
                length: 20.0,
            },
        ],
    ));
    let seams = exact_station_seams2(&chain).unwrap();
    let at: Vec<Scalar> = seams.iter().map(|s| s.distance).collect();
    assert_eq!(at, vec![10.0, 15.0]);
    assert!(seams.iter().all(|s| s.smooth && s.exact));
    // Smooth seams are read like any other station: both sides agree.
    assert_eq!(
        station_section2_on(&chain, 10.0, SeamSide::Incoming).unwrap(),
        station_section2(&chain, 10.0).unwrap()
    );
    // An elevated curve on that chain: the plan's seams and the profile's.
    let elevated = Curve3::Elevated(Elevated3::new(
        chain,
        ElevationLaw::Piecewise {
            breaks: vec![15.0, 25.0],
            laws: vec![
                ElevationLaw::level(0.0),
                ElevationLaw::constant_grade(0.0, 0.01),
                ElevationLaw::level(0.1),
            ],
        },
    ));
    let seams = exact_station_seams3(&elevated).unwrap();
    let at: Vec<(Scalar, bool)> = seams.iter().map(|s| (s.distance, s.smooth)).collect();
    assert_eq!(at, vec![(10.0, true), (15.0, false), (25.0, false)]);
    // A line has none; an unmeasurable family is refused as stations are.
    let line = Curve2::Line(Line2 {
        origin: Point2::ZERO,
        direction: Vec2::X,
    });
    assert!(station_seams2(&line).unwrap().is_empty());
}

#[test]
fn a_mitre_bisects_the_turn_and_meets_both_extrusions() {
    let curve = ell2();
    let incoming = station_section2_on(&curve, 3.0, SeamSide::Incoming).unwrap();
    let outgoing = station_section2_on(&curve, 3.0, SeamSide::Outgoing).unwrap();
    let mitre = Mitre::between(&incoming, &outgoing).unwrap().unwrap();
    let half = core::f64::consts::FRAC_1_SQRT_2;
    close3(mitre.normal, Vec3::new(half, half, 0.0), EPS, "bisector");
    // A point one unit to the left and half a unit up: on the first leg's
    // prism at (3 - 1, 1, 0.5), on the second's at the same point.
    let left_in = incoming.place(1.0, 0.5, 0.0);
    let left_out = outgoing.place(1.0, 0.5, 0.0);
    let placed = mitre.place(left_in, left_out, 0.0);
    close3(placed, Point3::new(2.0, 1.0, 0.5), EPS, "inner corner");
    assert!(((placed - incoming.point).dot(mitre.normal)).abs() <= EPS);
    // To the right: the outer corner.
    let placed = mitre.place(
        incoming.place(-1.0, 0.0, 0.0),
        outgoing.place(-1.0, 0.0, 0.0),
        0.0,
    );
    close3(placed, Point3::new(4.0, -1.0, 0.0), EPS, "outer corner");
    // A longitudinal offset moves the section along the bisector.
    let placed = mitre.place(incoming.point, outgoing.point, 2.0);
    close3(
        placed,
        Point3::new(3.0 + 2.0 * half, 2.0 * half, 0.0),
        EPS,
        "along the normal",
    );
}

#[test]
fn a_straight_seam_has_no_mitre_and_a_reversal_none_at_all() {
    let straight = Curve2::Polyline(Polyline2 {
        points: vec![
            Point2::new(0.0, 0.0),
            Point2::new(3.0, 0.0),
            Point2::new(5.0, 0.0),
        ],
        closed: false,
    });
    let incoming = station_section2_on(&straight, 3.0, SeamSide::Incoming).unwrap();
    let outgoing = station_section2_on(&straight, 3.0, SeamSide::Outgoing).unwrap();
    assert!(Mitre::between(&incoming, &outgoing).unwrap().is_none());
    // Within the tangent tolerance of straight: still no mitre.
    let mut nudged = outgoing;
    let angle = 0.5 * SEAM_TANGENT_TOLERANCE;
    nudged.tangent = Vec3::new(angle.cos(), angle.sin(), 0.0);
    assert!(Mitre::between(&incoming, &nudged).unwrap().is_none());
    // Turning back on itself.
    let back = Curve2::Polyline(Polyline2 {
        points: vec![
            Point2::new(0.0, 0.0),
            Point2::new(3.0, 0.0),
            Point2::new(1.0, 1e-9),
        ],
        closed: false,
    });
    let incoming = station_section2_on(&back, 3.0, SeamSide::Incoming).unwrap();
    let outgoing = station_section2_on(&back, 3.0, SeamSide::Outgoing).unwrap();
    assert!(0.5 * (incoming.tangent + outgoing.tangent).length() <= MITRE_TOLERANCE);
    refused(
        Mitre::between(&incoming, &outgoing).unwrap_err(),
        "turns back on itself",
    );
}
