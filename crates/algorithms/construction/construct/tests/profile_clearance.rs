//! Opt-in clearance of `ExactProfileValidator::validate_with` (#388).
//!
//! The overlay decides which lowered edges are closer than the distance; the
//! validator names the segments an author wrote, merges the sub-arcs of one
//! arc into one finding, and leaves the default report unchanged.

use axiolid_construct::extrude::extrude_profile_exact;
use axiolid_construct::profile_validate::{ContourRef, ExactProfileValidator, ProfileFinding};
use axiolid_core::{Frame2, Interval, Point2, Tolerance, Vec2, Vec3};
use axiolid_curve::{Circle2, Curve2, Line2};
use axiolid_overlay::ArcCheckOptions;
use axiolid_profile::{Contour, ContourProfile, Profile, ProfileSegment, ValidateProfile};

const TAU: f64 = core::f64::consts::TAU;
const PI: f64 = core::f64::consts::PI;
const TOL: f64 = 1e-3;

fn tol() -> Tolerance {
    Tolerance::new(TOL, TOL).expect("valid tolerance")
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

fn rect(x0: f64, y0: f64, x1: f64, y1: f64) -> Contour {
    let points = [p(x0, y0), p(x1, y0), p(x1, y1), p(x0, y1)];
    Contour::new(
        (0..4)
            .map(|i| line(points[i], points[(i + 1) % 4]))
            .collect(),
    )
}

/// A full circle as two half-turn segments: 0 the upper half, 1 the lower.
fn disc(centre: Point2, radius: f64) -> Contour {
    Contour::new(vec![
        arc(centre, radius, 0.0, PI),
        arc(centre, radius, PI, TAU),
    ])
}

fn plate(holes: Vec<Contour>) -> ContourProfile {
    ContourProfile {
        outer: rect(0.0, 0.0, 10.0, 10.0),
        holes,
    }
}

fn clearance() -> ArcCheckOptions {
    ArcCheckOptions::new().with_clearance(tol().linear())
}

fn clearances(profile: &ContourProfile) -> Vec<([(ContourRef, usize); 2], bool)> {
    ExactProfileValidator
        .validate_contours_with(profile, tol(), clearance())
        .findings
        .into_iter()
        .filter_map(|finding| match finding {
            ProfileFinding::Clearance {
                segments,
                undecided,
                ..
            } => Some((segments, undecided)),
            _ => None,
        })
        .collect()
}

const OUTER: ContourRef = ContourRef::Outer;
const HOLE: ContourRef = ContourRef::Hole(0);

#[test]
fn a_hole_near_the_outer_contour_names_the_segments() {
    let close = plate(vec![rect(0.5 * TOL, 1.0, 5.0, 4.0)]);
    assert_eq!(
        clearances(&close),
        vec![
            ([(OUTER, 3), (HOLE, 0)], false),
            ([(OUTER, 3), (HOLE, 2)], false),
            ([(OUTER, 3), (HOLE, 3)], false),
        ]
    );
    let apart = plate(vec![rect(2.0 * TOL, 1.0, 5.0, 4.0)]);
    let report = ExactProfileValidator.validate_contours_with(&apart, tol(), clearance());
    assert!(report.is_valid(), "{report:?}");
}

#[test]
fn the_sub_arcs_of_one_arc_make_one_finding() {
    // The lower half turn lowers to two sub-arcs; both rise within half the
    // tolerance of the bottom side only at the bottom point, but the pair
    // is one segment pair either way.
    let close = plate(vec![disc(p(5.0, 0.5 + 0.5 * TOL), 0.5)]);
    assert_eq!(clearances(&close), vec![([(OUTER, 0), (HOLE, 1)], false)]);
    let apart = plate(vec![disc(p(5.0, 0.5 + 2.0 * TOL), 0.5)]);
    assert!(clearances(&apart).is_empty());
}

#[test]
fn a_lower_bound_is_below_the_distance() {
    let close = plate(vec![rect(0.5 * TOL, 1.0, 5.0, 4.0)]);
    let report = ExactProfileValidator.validate_contours_with(&close, tol(), clearance());
    assert!(!report.is_valid());
    for finding in &report.findings {
        let ProfileFinding::Clearance { lower_bound, .. } = finding else {
            panic!("only clearance findings: {finding:?}");
        };
        assert!(
            (0.4 * TOL..0.5 * TOL).contains(lower_bound),
            "{lower_bound}"
        );
    }
}

#[test]
fn the_default_report_is_unchanged_and_extrusion_ignores_clearance() {
    let close = plate(vec![rect(0.5 * TOL, 1.0, 5.0, 4.0)]);
    let plain = ExactProfileValidator.validate(&close, tol());
    assert!(plain.is_valid(), "{plain:?}");
    assert_eq!(
        ExactProfileValidator.validate_contours_with(&close, tol(), ArcCheckOptions::new()),
        plain
    );
    let profile = Profile::Contour(close);
    assert_eq!(
        ExactProfileValidator.validate_with(&profile, tol(), ArcCheckOptions::new()),
        plain
    );
    extrude_profile_exact(&profile, Vec3::Z, 1.0, tol()).expect("clearance is not extruded on");
}

#[test]
fn a_merged_bound_holds_for_every_sub_arc() {
    // The disc's lower half lowers to sub-arcs of different distance to the
    // bottom side; the segment pair's bound must not exceed the nearest.
    let g = 1e-4;
    let c = p(5.0, 1.0 + g);
    let hole = Contour::new(vec![
        arc(c, 1.0, 1.5 * PI, 2.5 * PI),
        arc(c, 1.0, 0.5 * PI, 1.5 * PI),
    ]);
    let report = ExactProfileValidator.validate_contours_with(
        &plate(vec![hole]),
        tol(),
        ArcCheckOptions::new().with_clearance(1.2),
    );
    let bounds: Vec<(ContourRef, f64)> = report
        .findings
        .iter()
        .filter_map(|finding| match finding {
            ProfileFinding::Clearance {
                segments,
                lower_bound,
                ..
            } if segments[0] == (OUTER, 0) => Some((segments[1].0, *lower_bound)),
            _ => None,
        })
        .collect();
    assert!(!bounds.is_empty(), "{report:?}");
    for (_, bound) in bounds {
        assert!(bound <= g, "bound {bound} exceeds the distance {g}");
    }
}

#[test]
fn a_profile_is_judged_after_it_lowers() {
    let profile = Profile::Contour(plate(vec![rect(0.5 * TOL, 1.0, 5.0, 4.0)]));
    let report = ExactProfileValidator.validate_with(&profile, tol(), clearance());
    assert_eq!(report.findings.len(), 3, "{report:?}");
}
