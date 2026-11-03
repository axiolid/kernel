//! `extrude_profile_exact` refuses profiles the exact validator flags (#401).
//!
//! Each refusal fixture builds a solid that is not a valid manifold if the
//! extruder is not stopped first. The valid fixtures check the refusal does
//! not reach profiles the validator accepts.

use axiolid_construct::extrude::extrude_profile_exact;
use axiolid_construct::profile_validate::ExactProfileValidator;
use axiolid_contracts::GeomError;
use axiolid_core::{Frame2, Interval, Point2, Tolerance, Vec2, Vec3};
use axiolid_curve::{Circle2, Curve2, Line2};
use axiolid_profile::{Contour, ContourProfile, Profile, ProfileSegment, ValidateProfile};

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

fn arc(centre: Point2, radius: f64, from: f64, to: f64) -> ProfileSegment {
    ProfileSegment {
        curve: Curve2::Circle(Circle2 {
            frame: Frame2 {
                origin: centre,
                x: Vec2::X,
                y: Vec2::Y,
            },
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

fn disc(centre: Point2, radius: f64) -> Contour {
    Contour::new(vec![
        arc(centre, radius, 0.0, PI),
        arc(centre, radius, PI, TAU),
    ])
}

fn profile(outer: Contour, holes: Vec<Contour>) -> Profile {
    Profile::Contour(ContourProfile { outer, holes })
}

fn extrude(profile: &Profile) -> Result<(), GeomError> {
    extrude_profile_exact(profile, Vec3::Z, 1.0, tol()).map(|_| ())
}

/// The validator flags the fixture and the extruder refuses it as
/// `InvalidInput` naming the finding.
fn assert_invalid(profile: &Profile, named: &str) {
    assert!(
        !ExactProfileValidator.validate(profile, tol()).is_valid(),
        "the fixture must be flagged by the validator"
    );
    match extrude(profile) {
        Err(GeomError::InvalidInput(message)) => {
            assert!(message.contains(named), "{message}");
        }
        other => panic!("expected InvalidInput naming {named}, got {other:?}"),
    }
}

#[test]
fn a_chord_crossing_an_arc_is_refused() {
    let outer = Contour::new(vec![
        arc(p(0.0, 0.0), 1.0, 0.0, PI),
        line(p(-1.0, 0.0), p(1.5, 1.0)),
        line(p(1.5, 1.0), p(1.0, 0.0)),
    ]);
    assert_invalid(&profile(outer, vec![]), "SelfContact");
}

#[test]
fn a_bow_tie_is_refused() {
    let outer = polygon(&[p(0.0, 0.0), p(2.0, 2.0), p(2.0, 0.0), p(0.0, 2.0)]);
    assert_invalid(&profile(outer, vec![]), "SelfContact");
}

#[test]
fn a_pinch_is_refused() {
    // The vertex (2, 0) touches the non-adjacent edge (0, 0)-(4, 0).
    let outer = polygon(&[
        p(0.0, 0.0),
        p(4.0, 0.0),
        p(4.0, 4.0),
        p(2.0, 0.0),
        p(0.0, 4.0),
    ]);
    assert_invalid(&profile(outer, vec![]), "SelfContact");
}

#[test]
fn a_hole_outside_the_outer_contour_is_refused() {
    let plate = profile(
        rect(0.0, 0.0, 10.0, 10.0),
        vec![rect(20.0, 20.0, 22.0, 22.0)],
    );
    assert_invalid(&plate, "HoleOutsideOuter");
}

#[test]
fn two_crossing_holes_are_refused() {
    let plate = profile(
        rect(0.0, 0.0, 10.0, 10.0),
        vec![rect(1.0, 1.0, 5.0, 5.0), rect(3.0, 3.0, 8.0, 8.0)],
    );
    assert_invalid(&plate, "HolesOverlap");
}

#[test]
fn a_hole_crossing_the_outer_contour_is_refused() {
    let plate = profile(rect(0.0, 0.0, 10.0, 10.0), vec![rect(8.0, 4.0, 12.0, 6.0)]);
    assert_invalid(&plate, "HoleOutsideOuter");
}

#[test]
fn a_hole_crossing_an_arc_outer_is_refused() {
    let plate = profile(disc(p(0.0, 0.0), 2.0), vec![rect(1.0, -0.5, 3.0, 0.5)]);
    assert!(!ExactProfileValidator.validate(&plate, tol()).is_valid());
    assert!(matches!(extrude(&plate), Err(GeomError::InvalidInput(_))));
}

#[test]
fn a_refusal_holds_below_the_profile_plane() {
    let bow_tie = profile(
        polygon(&[p(0.0, 0.0), p(2.0, 2.0), p(2.0, 0.0), p(0.0, 2.0)]),
        vec![],
    );
    let result = extrude_profile_exact(&bow_tie, Vec3::Z, -1.0, tol());
    assert!(
        matches!(result, Err(GeomError::InvalidInput(_))),
        "{result:?}"
    );
}

#[test]
fn valid_profiles_still_extrude() {
    let plates = [
        profile(rect(0.0, 0.0, 4.0, 3.0), vec![]),
        profile(rect(-2.0, -2.0, 2.0, 2.0), vec![disc(p(0.0, 0.0), 1.0)]),
        profile(
            disc(p(0.0, 0.0), 3.0),
            vec![rect(-1.0, -1.0, 0.0, 1.0), disc(p(1.5, 0.0), 0.5)],
        ),
    ];
    for plate in &plates {
        assert!(ExactProfileValidator.validate(plate, tol()).is_valid());
        extrude(plate).expect("a valid profile extrudes");
    }
}

#[test]
fn a_tangent_contact_the_validator_accepts_still_extrudes() {
    // A disc hole tangent to the outer disc from inside: one isolated
    // contact, valid under the OGC rule.
    let tangent = profile(disc(p(0.0, 0.0), 2.0), vec![disc(p(1.0, 0.0), 1.0)]);
    let report = ExactProfileValidator.validate(&tangent, tol());
    assert!(report.is_valid() && !report.is_strict(), "{report:?}");
    extrude(&tangent).expect("an accepted contact extrudes unchanged");
}

fn composite(member: Profile) -> Profile {
    Profile::Composite(vec![member])
}

#[test]
fn a_composite_member_that_is_not_a_region_is_refused() {
    let pinch = polygon(&[
        p(0.0, 0.0),
        p(4.0, 0.0),
        p(4.0, 4.0),
        p(2.0, 0.0),
        p(0.0, 4.0),
    ]);
    let members = [
        (
            profile(
                rect(0.0, 0.0, 10.0, 10.0),
                vec![rect(20.0, 20.0, 22.0, 22.0)],
            ),
            "HoleOutsideOuter",
        ),
        (
            profile(rect(0.0, 0.0, 10.0, 10.0), vec![rect(8.0, 4.0, 12.0, 6.0)]),
            "HoleOutsideOuter",
        ),
        (
            profile(
                rect(0.0, 0.0, 10.0, 10.0),
                vec![rect(1.0, 1.0, 5.0, 5.0), rect(3.0, 3.0, 8.0, 8.0)],
            ),
            "HolesOverlap",
        ),
        (profile(pinch, vec![]), "SelfContact"),
    ];
    for (member, named) in members {
        match extrude(&composite(member)) {
            Err(GeomError::InvalidInput(message)) => {
                assert!(message.contains(named), "{named}: {message}");
            }
            other => panic!("expected InvalidInput naming {named}, got {other:?}"),
        }
    }
}

#[test]
fn a_valid_composite_still_extrudes() {
    let members = vec![
        profile(rect(0.0, 0.0, 4.0, 3.0), vec![]),
        profile(rect(2.0, 1.0, 6.0, 5.0), vec![disc(p(4.0, 3.0), 0.5)]),
    ];
    extrude(&Profile::Composite(members)).expect("a valid composite extrudes");
}

#[test]
fn refusal_classes_are_kept() {
    // An ellipse hole stays unsupported, an empty arc and a ring along one
    // line stay degenerate: the lowering and orientation refuse them as
    // before, whatever the region check would say.
    let ellipse = ProfileSegment {
        curve: Curve2::Ellipse(axiolid_curve::Ellipse2 {
            frame: Frame2 {
                origin: p(5.0, 5.0),
                x: Vec2::X,
                y: Vec2::Y,
            },
            semi_axis_x: 2.0,
            semi_axis_y: 1.0,
        }),
        domain: Interval::new(0.0, PI),
        same_sense: true,
    };
    let ellipse_hole = Contour::new(vec![ellipse, line(p(3.0, 5.0), p(7.0, 5.0))]);
    assert!(matches!(
        extrude(&profile(rect(0.0, 0.0, 10.0, 10.0), vec![ellipse_hole])),
        Err(GeomError::UnsupportedInput { .. })
    ));
    let empty_arc = Contour::new(vec![
        arc(p(5.0, 5.0), 1.0, 1.0, 1.0),
        line(p(5.0, 5.0), p(6.0, 6.0)),
    ]);
    assert!(matches!(
        extrude(&profile(rect(0.0, 0.0, 10.0, 10.0), vec![empty_arc])),
        Err(GeomError::Degenerate(_))
    ));
    let flat = polygon(&[p(0.0, 0.0), p(1.0, 0.0), p(2.0, 0.0)]);
    assert!(matches!(
        extrude(&profile(flat, vec![])),
        Err(GeomError::Degenerate(_))
    ));
    let flat_hole = polygon(&[p(1.0, 1.0), p(2.0, 1.0), p(3.0, 1.0)]);
    assert!(matches!(
        extrude(&profile(rect(0.0, 0.0, 10.0, 10.0), vec![flat_hole])),
        Err(GeomError::Degenerate(_))
    ));
}

#[test]
fn edges_and_arcs_below_the_tolerance_still_extrude() {
    // A vertex 5e-7 short of a corner and a corner arc of radius 5e-7, at
    // tolerance 1e-6: the exact solid is valid and extruded before the check.
    let short = polygon(&[
        p(0.0, 0.0),
        p(4.0, 0.0),
        p(4.0, 3.0),
        p(4.0 - 5e-7, 3.0),
        p(0.0, 3.0),
    ]);
    extrude(&profile(short, vec![])).expect("a sub-tolerance edge extrudes");
    let r = 5e-7;
    let fillet = Contour::new(vec![
        line(p(0.0, 0.0), p(4.0, 0.0)),
        line(p(4.0, 0.0), p(4.0, 3.0 - r)),
        arc(p(4.0 - r, 3.0 - r), r, 0.0, PI / 2.0),
        line(p(4.0 - r, 3.0), p(0.0, 3.0)),
        line(p(0.0, 3.0), p(0.0, 0.0)),
    ]);
    extrude(&profile(fillet, vec![])).expect("a sub-tolerance arc extrudes");
}

#[test]
fn contacts_that_cut_the_interior_apart_are_refused() {
    let diamond = polygon(&[p(0.0, 5.0), p(5.0, 0.0), p(10.0, 5.0), p(5.0, 10.0)]);
    assert_invalid(
        &profile(rect(0.0, 0.0, 10.0, 10.0), vec![diamond]),
        "DisconnectedInterior",
    );
}
