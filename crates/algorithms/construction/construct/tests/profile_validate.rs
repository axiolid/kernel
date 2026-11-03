//! `ExactProfileValidator` against fixtures with known answers (#352).
//!
//! The validator judges the arc ring the exact path lowers a contour to, so
//! a fixture's verdict is checked two ways: against what the geometry
//! says, and (for the valid ones) against the exact extruder building it.

use axiolid_construct::extrude::extrude_profile_exact;
use axiolid_construct::profile_validate::{
    ContourRef, ExactProfileValidator, ProfileFinding, ProfileReport,
};
use axiolid_construct::section_lower::{circle_contour, rectangle_contour, section_contour};
use axiolid_contracts::GeomError;
use axiolid_core::{Frame2, Interval, Point2, Tolerance, Transform2, Vec2, Vec3};
use axiolid_curve::{Circle2, Curve2, Ellipse2, Line2, Polyline2};
use axiolid_overlay::OverlayError;
use axiolid_profile::{
    CenterLineProfile, CircleProfile, Contour, ContourProfile, EllipseProfile, Profile,
    ProfileSegment, RectangleProfile, SectionProfile, ValidateProfile,
};

const TAU: f64 = core::f64::consts::TAU;
const PI: f64 = core::f64::consts::PI;

fn tol() -> Tolerance {
    Tolerance::new(1e-6, 1e-6).expect("valid tolerance")
}

fn p(x: f64, y: f64) -> Point2 {
    Point2::new(x, y)
}

fn line(from: Point2, to: Point2) -> ProfileSegment {
    ProfileSegment {
        curve: Curve2::Line(Line2 {
            origin: from,
            direction: to - from,
        }),
        domain: Interval::UNIT,
        same_sense: true,
    }
}

fn frame(centre: Point2) -> Frame2 {
    Frame2 {
        origin: centre,
        x: Vec2::X,
        y: Vec2::Y,
    }
}

fn arc(centre: Point2, radius: f64, from: f64, to: f64) -> ProfileSegment {
    ProfileSegment {
        curve: Curve2::Circle(Circle2 {
            frame: frame(centre),
            radius,
        }),
        domain: Interval::new(from, to),
        same_sense: true,
    }
}

fn polygon(points: &[Point2]) -> Contour {
    Contour::new(
        (0..points.len())
            .map(|i| line(points[i], points[(i + 1) % points.len()]))
            .collect(),
    )
}

fn rect(x0: f64, y0: f64, x1: f64, y1: f64) -> Contour {
    polygon(&[p(x0, y0), p(x1, y0), p(x1, y1), p(x0, y1)])
}

/// A full circle as two half-turn segments.
fn disc(centre: Point2, radius: f64) -> Contour {
    Contour::new(vec![
        arc(centre, radius, 0.0, PI),
        arc(centre, radius, PI, TAU),
    ])
}

fn contours(outer: Contour, holes: Vec<Contour>) -> ContourProfile {
    ContourProfile { outer, holes }
}

fn validate(profile: &ContourProfile) -> ProfileReport {
    ExactProfileValidator.validate(profile, tol())
}

fn validate_profile(profile: &Profile, tolerance: Tolerance) -> ProfileReport {
    ExactProfileValidator.validate(profile, tolerance)
}

fn near(a: Point2, x: f64, y: f64) -> bool {
    (a.x - x).abs() < 1e-9 && (a.y - y).abs() < 1e-9
}

// -------------------------------------------------------------- closure

#[test]
fn an_open_contour_names_the_closing_segment() {
    let open = polygon_open(&[p(0.0, 0.0), p(2.0, 0.0), p(2.0, 2.0), p(0.5, 2.0)]);
    let report = validate(&contours(open, vec![]));
    match report.findings.as_slice() {
        [ProfileFinding::Gap {
            contour: ContourRef::Outer,
            segment: 2,
            gap,
        }] => assert!((gap - 0.5_f64.hypot(2.0)).abs() < 1e-12, "{gap}"),
        other => panic!("{other:?}"),
    }
}

/// Segments between consecutive points only: the last point does not close.
fn polygon_open(points: &[Point2]) -> Contour {
    Contour::new(points.windows(2).map(|w| line(w[0], w[1])).collect())
}

fn square_with_gap(gap: f64) -> Contour {
    // The second segment starts `gap` to the right of where the first ends.
    Contour::new(vec![
        line(p(0.0, 0.0), p(2.0, 0.0)),
        line(p(2.0 + gap, 0.0), p(2.0, 2.0)),
        line(p(2.0, 2.0), p(0.0, 2.0)),
        line(p(0.0, 2.0), p(0.0, 0.0)),
    ])
}

#[test]
fn joint_gaps_are_judged_against_the_tolerance() {
    let wide = validate(&contours(square_with_gap(2e-6), vec![]));
    assert!(
        matches!(
            wide.findings.as_slice(),
            [ProfileFinding::Gap {
                segment: 0,
                contour: ContourRef::Outer,
                ..
            }]
        ),
        "{wide:?}"
    );
    let narrow = validate(&contours(square_with_gap(5e-7), vec![]));
    assert!(narrow.is_valid(), "{narrow:?}");
    assert!((narrow.max_joint_gap - 5e-7).abs() < 1e-15, "{narrow:?}");
    // The closing joint is a joint too.
    let closing = Contour::new(vec![
        line(p(0.0, 0.0), p(2.0, 0.0)),
        line(p(2.0, 0.0), p(2.0, 2.0)),
        line(p(2.0, 2.0), p(0.0, 2.0)),
        line(p(0.0, 2.0), p(0.0, 5e-7)),
    ]);
    let closed = validate(&contours(closing, vec![]));
    assert!(closed.is_valid(), "{closed:?}");
    assert!((closed.max_joint_gap - 5e-7).abs() < 1e-15, "{closed:?}");
    let exact = validate(&contours(rect(0.0, 0.0, 2.0, 2.0), vec![]));
    assert!(exact.is_strict() && exact.max_joint_gap == 0.0, "{exact:?}");
}

fn rolled_sections() -> Vec<(&'static str, Profile)> {
    let section = |s| Profile::Section(s);
    vec![
        (
            "IPE 300 with root and toe radii",
            section(SectionProfile::I {
                depth: 0.3,
                width: 0.15,
                web_thickness: 0.0071,
                flange_thickness: 0.0107,
                fillet_radius: Some(0.015),
                flange_edge_radius: Some(0.0035),
                flange_slope: None,
            }),
        ),
        (
            "IPE 300 sharp",
            section(SectionProfile::I {
                depth: 0.3,
                width: 0.15,
                web_thickness: 0.0071,
                flange_thickness: 0.0107,
                fillet_radius: None,
                flange_edge_radius: None,
                flange_slope: None,
            }),
        ),
        (
            "T 100",
            section(SectionProfile::T {
                depth: 0.1,
                flange_width: 0.1,
                web_thickness: 0.011,
                flange_thickness: 0.011,
                fillet_radius: Some(0.011),
                flange_edge_radius: Some(0.0055),
                web_edge_radius: Some(0.003),
                web_slope: None,
                flange_slope: None,
            }),
        ),
        (
            "UPE 200",
            section(SectionProfile::U {
                depth: 0.2,
                flange_width: 0.08,
                web_thickness: 0.006,
                flange_thickness: 0.011,
                fillet_radius: Some(0.013),
                edge_radius: Some(0.0045),
                flange_slope: None,
            }),
        ),
        (
            "L 120x80x10",
            section(SectionProfile::L {
                depth: 0.12,
                width: Some(0.08),
                thickness: 0.01,
                fillet_radius: Some(0.011),
                edge_radius: Some(0.0055),
                leg_slope: None,
            }),
        ),
        (
            "rounded rectangle",
            Profile::Rectangle(RectangleProfile {
                x: 0.2,
                y: 0.1,
                thickness: None,
                outer_radius: Some(0.0158),
                inner_radius: None,
            }),
        ),
        (
            "RHS 200x100x6.3",
            Profile::Rectangle(RectangleProfile {
                x: 0.2,
                y: 0.1,
                thickness: Some(0.0063),
                outer_radius: Some(0.0158),
                inner_radius: Some(0.0095),
            }),
        ),
        (
            "CHS 168.3x7.1",
            Profile::Circle(CircleProfile {
                radius: 0.08415,
                thickness: Some(0.0071),
            }),
        ),
        (
            "circle",
            Profile::Circle(CircleProfile {
                radius: 1.5,
                thickness: None,
            }),
        ),
        (
            "rectangle",
            Profile::Rectangle(RectangleProfile {
                x: 2.0,
                y: 1.0,
                thickness: None,
                outer_radius: None,
                inner_radius: None,
            }),
        ),
    ]
}

#[test]
fn the_validator_agrees_with_the_exact_path_at_zero_tolerance() {
    for (name, profile) in rolled_sections() {
        let report = validate_profile(&profile, Tolerance::ZERO);
        assert!(report.is_valid(), "{name}: {report:?}");
        // Extrude the lowered contour: the exact path has no annular circle.
        let contour = match &profile {
            Profile::Section(section) => section_contour(section),
            Profile::Rectangle(rectangle) => rectangle_contour(rectangle),
            Profile::Circle(circle) => circle_contour(circle),
            other => panic!("{name}: not a parametric profile: {other:?}"),
        }
        .expect("a parametric profile lowers");
        extrude_profile_exact(&Profile::Contour(contour), Vec3::Z, 1.0, Tolerance::ZERO)
            .unwrap_or_else(|error| panic!("{name}: the exact path refuses it: {error}"));
    }
}

#[test]
fn derived_and_centre_line_profiles_are_valid() {
    let mirror = Transform2::from_cols(
        Vec2::new(2.0, 0.0),
        Vec2::new(0.0, -2.0),
        Vec2::new(3.0, -1.0),
    );
    let derived = Profile::Derived {
        basis: Box::new(Profile::Circle(CircleProfile {
            radius: 1.5,
            thickness: None,
        })),
        transform: mirror,
    };
    let report = validate_profile(&derived, tol());
    assert!(report.is_valid(), "{report:?}");
    extrude_profile_exact(&derived, Vec3::Z, 1.0, tol()).expect("the exact path builds it");

    let centre = Profile::CenterLine(CenterLineProfile::from_width(
        Contour::new(vec![line(p(0.0, 0.0), p(2.0, 0.0))]),
        0.4,
    ));
    let report = validate_profile(&centre, tol());
    assert!(report.is_valid(), "{report:?}");
    extrude_profile_exact(&centre, Vec3::Z, 1.0, tol()).expect("the exact path builds it");
}

// ----------------------------------------------------------------- arcs

/// Two circles of radius 2 and 1.6, 0.5 apart, cut a crescent: the outer
/// arc runs round the big circle's left, the inner back along the small one.
fn crescent() -> Contour {
    let (r1, r2, d) = (2.0_f64, 1.6_f64, 0.5_f64);
    let x = (r1 * r1 - r2 * r2 + d * d) / (2.0 * d);
    let y = (r1 * r1 - x * x).sqrt();
    let a = y.atan2(x);
    let b = y.atan2(x - d);
    Contour::new(vec![
        arc(p(0.0, 0.0), r1, a, TAU - a),
        ProfileSegment {
            same_sense: false,
            ..arc(p(d, 0.0), r2, b, TAU - b)
        },
    ])
}

#[test]
fn a_crescent_is_valid_and_strict() {
    let profile = contours(crescent(), vec![]);
    let report = validate(&profile);
    assert!(report.is_strict(), "{report:?}");
    assert!(report.max_joint_gap < 1e-14, "{report:?}");
    extrude_profile_exact(&Profile::Contour(profile), Vec3::Z, 1.0, tol())
        .expect("the exact path builds it");
}

#[test]
fn a_chord_crossing_an_arc_names_both_segments() {
    // Upper half circle, then a line leaving the circle through the arc.
    let outer = Contour::new(vec![
        arc(p(0.0, 0.0), 1.0, 0.0, PI),
        line(p(-1.0, 0.0), p(1.5, 1.0)),
        line(p(1.5, 1.0), p(1.0, 0.0)),
    ]);
    let report = validate(&contours(outer, vec![]));
    let found = report.findings.iter().find_map(|f| match f {
        ProfileFinding::SelfContact {
            contour: ContourRef::Outer,
            segments,
            at,
        } => Some((*segments, *at)),
        _ => None,
    });
    let (segments, at) = found.unwrap_or_else(|| panic!("{report:?}"));
    assert_eq!(segments, [0, 1]);
    // t = 5 / 7.25 along the line.
    let t = 5.0 / 7.25;
    assert!(near(at, -1.0 + 2.5 * t, t), "{at:?}");
}

#[test]
fn a_contact_on_a_split_arc_names_the_arc_not_a_sub_edge() {
    // Three quarters of a turn lowers to three sub-arcs; a chord cuts the
    // third quarter and the contact is reported on segment 0.
    let outer = Contour::new(vec![
        arc(p(0.0, 0.0), 1.0, 0.0, 1.5 * PI),
        line(p(0.0, -1.0), p(-2.0, 1.5)),
        line(p(-2.0, 1.5), p(1.0, 0.0)),
    ]);
    let report = validate(&contours(outer, vec![]));
    let segments: Vec<[usize; 2]> = report
        .findings
        .iter()
        .filter_map(|f| match f {
            ProfileFinding::SelfContact { segments, .. } => Some(*segments),
            _ => None,
        })
        .collect();
    assert!(!segments.is_empty(), "{report:?}");
    for pair in segments {
        assert!(pair.iter().all(|&s| s <= 2) && pair[0] == 0, "{pair:?}");
    }
}

// ------------------------------------------------------------- families

#[test]
fn unsupported_segment_families_are_named() {
    let ellipse = ProfileSegment {
        curve: Curve2::Ellipse(Ellipse2 {
            frame: frame(p(1.0, 1.0)),
            semi_axis_x: 1.0,
            semi_axis_y: 0.5,
        }),
        domain: Interval::new(0.0, PI),
        same_sense: true,
    };
    let outer = Contour::new(vec![
        line(p(0.0, 0.0), p(1.0, 0.0)),
        ellipse,
        line(p(0.0, 1.0), p(0.0, 0.0)),
    ]);
    let report = validate(&contours(outer, vec![]));
    assert_eq!(
        report.findings,
        vec![ProfileFinding::UnsupportedSegment {
            contour: ContourRef::Outer,
            segment: 1,
            kind: "elliptical contour segment",
        }]
    );

    let polyline = ProfileSegment {
        curve: Curve2::Polyline(Polyline2 {
            points: vec![p(1.0, 0.0), p(1.0, 1.0)],
            closed: false,
        }),
        domain: Interval::UNIT,
        same_sense: true,
    };
    let outer = Contour::new(vec![line(p(0.0, 0.0), p(1.0, 0.0)), polyline]);
    let report = validate(&contours(outer, vec![]));
    assert!(
        matches!(
            report.findings.as_slice(),
            [ProfileFinding::UnsupportedSegment {
                segment: 1,
                kind: "polyline contour segment",
                ..
            }]
        ),
        "{report:?}"
    );
}

#[test]
fn profiles_that_do_not_lower_are_reported_whole() {
    let cases: Vec<(&str, Profile)> = vec![
        (
            "ellipse",
            Profile::Ellipse(EllipseProfile {
                semi_axis_x: 2.0,
                semi_axis_y: 1.0,
            }),
        ),
        (
            "composite",
            Profile::Composite(vec![Profile::Circle(CircleProfile {
                radius: 1.0,
                thickness: None,
            })]),
        ),
        (
            "negative rectangle",
            Profile::Rectangle(RectangleProfile {
                x: -1.0,
                y: 1.0,
                thickness: None,
                outer_radius: None,
                inner_radius: None,
            }),
        ),
    ];
    for (name, profile) in cases {
        let report = validate_profile(&profile, tol());
        assert!(
            matches!(report.findings.as_slice(), [ProfileFinding::NotLowered(_)]),
            "{name}: {report:?}"
        );
        assert!(report.contacts.is_empty());
    }
    let rectangle = validate_profile(
        &Profile::Rectangle(RectangleProfile {
            x: -1.0,
            y: 1.0,
            thickness: None,
            outer_radius: None,
            inner_radius: None,
        }),
        tol(),
    );
    assert!(matches!(
        rectangle.findings.as_slice(),
        [ProfileFinding::NotLowered(GeomError::InvalidInput(_))]
    ));
    let ellipse = validate_profile(
        &Profile::Ellipse(EllipseProfile {
            semi_axis_x: 2.0,
            semi_axis_y: 1.0,
        }),
        tol(),
    );
    assert!(matches!(
        ellipse.findings.as_slice(),
        [ProfileFinding::NotLowered(
            GeomError::UnsupportedInput { .. }
        )]
    ));
}

#[test]
fn a_contour_with_one_segment_is_refused() {
    let report = validate(&contours(
        Contour::new(vec![line(p(0.0, 0.0), p(1.0, 0.0))]),
        vec![],
    ));
    assert_eq!(
        report.findings,
        vec![ProfileFinding::TooFewSegments {
            contour: ContourRef::Outer,
            count: 1
        }]
    );
}

// ---------------------------------------------------------------- holes

#[test]
fn holes_are_nested_by_contour() {
    let outer = || rect(0.0, 0.0, 10.0, 10.0);

    let outside = validate(&contours(outer(), vec![rect(20.0, 20.0, 22.0, 22.0)]));
    assert_eq!(
        outside.findings,
        vec![ProfileFinding::HoleOutsideOuter { hole: 0 }]
    );

    let crossing = validate(&contours(
        outer(),
        vec![rect(1.0, 1.0, 5.0, 5.0), rect(3.0, 3.0, 8.0, 8.0)],
    ));
    assert_eq!(
        crossing.findings,
        vec![ProfileFinding::HolesOverlap { holes: [0, 1] }]
    );

    let nested = validate(&contours(
        outer(),
        vec![
            rect(7.0, 7.0, 8.0, 8.0),
            rect(1.0, 1.0, 5.0, 5.0),
            rect(2.0, 2.0, 3.0, 3.0),
        ],
    ));
    assert_eq!(
        nested.findings,
        vec![ProfileFinding::HolesOverlap { holes: [1, 2] }]
    );

    let sharing = validate(&contours(outer(), vec![rect(2.0, 0.0, 5.0, 5.0)]));
    assert_eq!(
        sharing.findings,
        vec![ProfileFinding::BoundaryOverlap {
            contours: [ContourRef::Outer, ContourRef::Hole(0)]
        }]
    );
}

#[test]
fn a_hole_touching_the_outer_is_a_contact() {
    let hole = polygon(&[p(8.0, 4.0), p(10.0, 5.0), p(8.0, 6.0)]);
    let report = validate(&contours(rect(0.0, 0.0, 10.0, 10.0), vec![hole]));
    assert!(report.is_valid() && !report.is_strict(), "{report:?}");
    assert_eq!(report.contacts.len(), 1);
    assert_eq!(
        report.contacts[0].contours,
        [ContourRef::Outer, ContourRef::Hole(0)]
    );
    assert!(near(report.contacts[0].at, 10.0, 5.0));
}

#[test]
fn a_plate_with_a_round_hole_is_valid_and_builds() {
    let profile = contours(rect(-2.0, -2.0, 2.0, 2.0), vec![disc(p(0.0, 0.0), 1.0)]);
    let report = validate(&profile);
    assert!(report.is_strict(), "{report:?}");
    extrude_profile_exact(&Profile::Contour(profile), Vec3::Z, 1.0, tol())
        .expect("the exact path builds it");
}

// ----------------------------------------------------------- degenerate

#[test]
fn degenerate_segments_and_rings_are_typed() {
    // A sub-tolerance radius: lowers, then the ring check refuses it.
    let tiny = validate(&contours(disc(p(0.0, 0.0), 0.75e-6), vec![]));
    assert!(
        matches!(
            tiny.findings.first(),
            Some(ProfileFinding::DegenerateRing {
                contour: ContourRef::Outer,
                error: OverlayError::ZeroRadiusArc,
            })
        ),
        "{tiny:?}"
    );

    let flat = ProfileSegment {
        curve: Curve2::Circle(Circle2 {
            frame: Frame2 {
                origin: p(0.0, 0.0),
                x: Vec2::X,
                y: Vec2::X,
            },
            radius: 1.0,
        }),
        domain: Interval::new(0.0, 1.0),
        same_sense: true,
    };
    let report = validate(&contours(
        Contour::new(vec![flat, line(p(1.0, 0.0), p(0.0, 0.0))]),
        vec![],
    ));
    assert!(
        matches!(
            report.findings.as_slice(),
            [ProfileFinding::DegenerateSegment {
                segment: 0,
                reason: "circular profile segment has a degenerate frame",
                ..
            }]
        ),
        "{report:?}"
    );

    let empty = arc(p(0.0, 0.0), 1.0, 1.0, 1.0);
    let report = validate(&contours(
        Contour::new(vec![empty, line(p(0.5, 0.8), p(0.0, 0.0))]),
        vec![],
    ));
    assert!(
        matches!(
            report.findings.as_slice(),
            [ProfileFinding::DegenerateSegment { segment: 0, .. }]
        ),
        "{report:?}"
    );

    let overshoot = arc(p(0.0, 0.0), 1.0, 0.0, 7.0);
    let report = validate(&contours(
        Contour::new(vec![overshoot, line(p(0.75, 0.66), p(1.0, 0.0))]),
        vec![],
    ));
    assert!(
        matches!(
            report.findings.as_slice(),
            [ProfileFinding::UnsupportedSegment {
                segment: 0,
                kind: "circular profile segment sweeping more than a whole turn",
                ..
            }]
        ),
        "{report:?}"
    );
}

#[test]
fn a_broken_contour_stops_the_nesting_check() {
    // The hole lies outside, but the outer does not lower: only the
    // outer's finding is reported.
    let open = polygon_open(&[p(0.0, 0.0), p(2.0, 0.0), p(2.0, 2.0), p(0.5, 2.0)]);
    let report = validate(&contours(open, vec![rect(20.0, 20.0, 22.0, 22.0)]));
    assert!(
        matches!(report.findings.as_slice(), [ProfileFinding::Gap { .. }]),
        "{report:?}"
    );
    // A hole that does not lower is named by its index.
    let report = validate(&contours(
        rect(0.0, 0.0, 10.0, 10.0),
        vec![
            rect(1.0, 1.0, 2.0, 2.0),
            polygon_open(&[p(4.0, 4.0), p(5.0, 4.0), p(5.0, 5.0)]),
        ],
    ));
    assert!(
        matches!(
            report.findings.as_slice(),
            [ProfileFinding::Gap {
                contour: ContourRef::Hole(1),
                ..
            }]
        ),
        "{report:?}"
    );
}

#[test]
fn self_contact_is_reported_per_contour_with_segment_indices() {
    // A bow tie in hole 1 only.
    let bow = polygon(&[p(1.0, 1.0), p(3.0, 3.0), p(3.0, 1.0), p(1.0, 4.0)]);
    let outer = rect(-10.0, -10.0, 0.0, 0.0);
    // The bow tie sits in the second hole's slot, outside the outer: only
    // the self-contact is reported, nesting is skipped.
    let report = validate(&contours(outer, vec![rect(-3.0, -3.0, -1.0, -1.0), bow]));
    match report.findings.as_slice() {
        [ProfileFinding::SelfContact {
            contour: ContourRef::Hole(1),
            segments: [0, 2],
            ..
        }] => {}
        other => panic!("{other:?}"),
    }
}

#[test]
fn reports_repeat() {
    let profile = contours(crescent(), vec![]);
    assert_eq!(validate(&profile), validate(&profile));
    let broken = contours(rect(0.0, 0.0, 10.0, 10.0), vec![rect(8.0, 8.0, 12.0, 12.0)]);
    assert_eq!(validate(&broken), validate(&broken));
}

#[test]
fn a_hole_names_its_own_segments_across_a_split_arc() {
    // Hole 1 is the arc-and-chord contour of the chord test, moved into the
    // plate: the semicircle lowers to two ring edges, so ring edge and
    // segment indices differ. Hole 0 is a plain square, whose indices agree.
    let moved = Contour::new(vec![
        arc(p(5.0, 5.0), 1.0, 0.0, PI),
        line(p(4.0, 5.0), p(6.5, 6.0)),
        line(p(6.5, 6.0), p(6.0, 5.0)),
    ]);
    let report = validate(&contours(
        rect(0.0, 0.0, 10.0, 10.0),
        vec![rect(1.0, 1.0, 2.0, 2.0), moved],
    ));
    match report.findings.as_slice() {
        [ProfileFinding::SelfContact {
            contour: ContourRef::Hole(1),
            segments: [0, 1],
            ..
        }] => {}
        other => panic!("{other:?}"),
    }
}
