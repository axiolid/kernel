//! Contour and arc extrusion with holes (ADR 0058).
//!
//! Volumes are checked against closed forms chosen so a DROPPED hole is
//! visible: the outer area alone differs from the correct answer by exactly
//! the hole area.

use axiolid_brep_audit::geometric_audit;
use axiolid_construct::contour_lower::{
    arc_ring_signed_area, contour_to_arc_ring, orient_arc_ring,
};
use axiolid_construct::extrude::extrude_profile_exact;
use axiolid_core::{Frame2, Interval, Point2, Tolerance, Vec2, Vec3};
use axiolid_curve::{Circle2, Curve2, Line2};
use axiolid_profile::{Contour, ContourProfile, Profile, ProfileSegment};
use axiolid_surface::Surface;

const DEPTH: f64 = 2.0;

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

/// A full circle as four quarter arcs.
///
/// Quarters rather than halves: ADR 0053 refuses a segment sweeping half a
/// turn or more, because `bulge = tan(sweep/4)` stops determining the arc
/// from its chord there.
fn circle(centre: Point2, radius: f64) -> Contour {
    let frame = Frame2 {
        origin: centre,
        x: Vec2::X,
        y: Vec2::Y,
    };
    let quarter = core::f64::consts::FRAC_PI_2;
    Contour::new(
        (0..4)
            .map(|index| ProfileSegment {
                curve: Curve2::Circle(Circle2 { frame, radius }),
                domain: Interval::new(quarter * index as f64, quarter * (index + 1) as f64),
                same_sense: true,
            })
            .collect(),
    )
}

fn square(half: f64) -> Contour {
    Contour::new(vec![
        line(Point2::new(-half, -half), Point2::new(half, -half)),
        line(Point2::new(half, -half), Point2::new(half, half)),
        line(Point2::new(half, half), Point2::new(-half, half)),
        line(Point2::new(-half, half), Point2::new(-half, -half)),
    ])
}

fn extrude(outer: Contour, holes: Vec<Contour>) -> axiolid_brep::ExactBRep {
    let profile = Profile::Contour(ContourProfile { outer, holes });
    let solid = extrude_profile_exact(&profile, Vec3::Z, DEPTH, Tolerance::METRE)
        .expect("a contour with holes extrudes");
    let health = geometric_audit(&solid, Tolerance::METRE);
    assert!(
        health.is_consistent(),
        "solid must audit clean, found {:?}",
        health.defects()
    );
    solid
}

#[test]
fn a_plate_with_a_round_hole_loses_exactly_the_hole_volume() {
    // 4x4 plate, radius-1 hole. A dropped hole would read 16*DEPTH; a hole
    // wound the wrong way would read 16+pi rather than 16-pi.
    let radius = 1.0;
    let solid = extrude(square(2.0), vec![circle(Point2::new(0.0, 0.0), radius)]);

    // `exact_properties` is planar-only, so the area comes from the exact
    // rings via the signed-area routine the extruder itself uses to decide
    // winding. Outer is positive, a correctly wound hole is negative, so the
    // sum is the net section -- and a dropped or mis-wound hole changes it.
    let outer_area = arc_ring_signed_area(
        &contour_to_arc_ring(&square(2.0), Tolerance::METRE).expect("outer lowers"),
    );
    let hole_area = arc_ring_signed_area(
        &orient_arc_ring(
            &contour_to_arc_ring(&circle(Point2::new(0.0, 0.0), radius), Tolerance::METRE)
                .expect("hole lowers"),
            false,
        )
        .expect("hole orients"),
    );
    let expected = 16.0 - core::f64::consts::PI * radius * radius;
    let got = outer_area + hole_area;
    assert!(
        (got - expected).abs() < 1e-9,
        "net section: expected {expected}, got {got}"
    );
    assert!(
        hole_area < 0.0,
        "a hole must wind clockwise, got {hole_area}"
    );

    // The hole wall is a genuine cylinder, not a fan of planes.
    let radii: Vec<f64> = solid
        .surfaces()
        .iter()
        .filter_map(|s| match s {
            Surface::Cylinder(c) => Some(c.radius),
            _ => None,
        })
        .collect();
    assert_eq!(radii.len(), 4, "four quarter walls, got {radii:?}");
    for value in &radii {
        assert!((value - radius).abs() < 1e-12, "got {value}");
    }
}

#[test]
fn a_hole_given_counter_clockwise_is_reoriented_not_rejected() {
    // An imported contour may hand holes in either winding. The extruder
    // orients them, so both inputs must build the SAME solid rather than one
    // building an inside-out passage.
    let radius = 0.8;
    let ccw = circle(Point2::new(0.0, 0.0), radius);
    let cw = Contour::new(reverse_segments(&ccw));

    let from_ccw = extrude(square(2.0), vec![ccw]);
    let from_cw = extrude(square(2.0), vec![cw]);

    // Same face count and same cylinder radii: the hole survived both ways.
    assert_eq!(
        from_ccw.topology().faces().len(),
        from_cw.topology().faces().len(),
        "both windings must give the same topology"
    );
    let walls = |solid: &axiolid_brep::ExactBRep| {
        solid
            .surfaces()
            .iter()
            .filter(|s| matches!(s, Surface::Cylinder(_)))
            .count()
    };
    assert_eq!(walls(&from_ccw), 4);
    assert_eq!(walls(&from_cw), 4);
}

/// Reverse a contour of arc segments, flipping each segment's sense.
fn reverse_segments(contour: &Contour) -> Vec<ProfileSegment> {
    contour
        .segments
        .iter()
        .rev()
        .map(|segment| ProfileSegment {
            curve: segment.curve.clone(),
            domain: segment.domain,
            same_sense: !segment.same_sense,
        })
        .collect()
}

#[test]
fn several_holes_each_contribute_their_own_walls() {
    let radius = 0.5;
    let solid = extrude(
        square(3.0),
        vec![
            circle(Point2::new(-1.5, 0.0), radius),
            circle(Point2::new(1.5, 0.0), radius),
        ],
    );
    let cylinders = solid
        .surfaces()
        .iter()
        .filter(|s| matches!(s, Surface::Cylinder(_)))
        .count();
    assert_eq!(cylinders, 8, "four quarter walls per hole, got {cylinders}");
}

#[test]
fn a_straight_outer_with_a_straight_hole_still_routes_to_the_polygon_path() {
    // All-straight sections keep their existing behaviour, including the
    // face naming the polygon path carries and the arc path does not.
    let solid = extrude(square(2.0), vec![square(1.0)]);
    let expected = (16.0 - 4.0) * DEPTH;
    let got = axiolid_measure::exact_properties(&solid, Tolerance::METRE)
        .map(|p| p.signed_volume)
        .expect("an all-planar prism is measurable");
    assert!(
        (got - expected).abs() < 1e-12,
        "expected {expected}, got {got}"
    );
}
