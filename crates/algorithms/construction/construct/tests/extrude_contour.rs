//! Extruding arbitrary exact contours (ADR 0053).
//!
//! The claim: a contour of lines and arcs extrudes to an exact solid with
//! genuine `Plane` and `Cylinder` walls, and anything not exactly
//! representable is refused rather than sampled into chords.

use axiolid_brep_audit::geometric_audit;
use axiolid_construct::extrude::extrude_profile_exact;
use axiolid_core::{Frame2, Interval, Point2, Tolerance, Vec2, Vec3};
use axiolid_curve::{Circle2, Curve2, Ellipse2, Line2};
use axiolid_profile::{Contour, ContourProfile, Profile, ProfileSegment};
use axiolid_surface::Surface;

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

/// A quarter arc on the unit-radius circle centred at `centre`.
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

fn profile(outer: Contour) -> Profile {
    Profile::Contour(ContourProfile {
        outer,
        holes: Vec::new(),
    })
}

#[test]
fn a_polygonal_contour_extrudes_to_planar_walls() {
    let contour = Contour::new(vec![
        line(Point2::new(0.0, 0.0), Point2::new(3.0, 0.0)),
        line(Point2::new(3.0, 0.0), Point2::new(3.0, 1.0)),
        line(Point2::new(3.0, 1.0), Point2::new(1.0, 1.0)),
        line(Point2::new(1.0, 1.0), Point2::new(1.0, 3.0)),
        line(Point2::new(1.0, 3.0), Point2::new(0.0, 3.0)),
        line(Point2::new(0.0, 3.0), Point2::new(0.0, 0.0)),
    ]);
    // An L-shape: six sides, which no rectangle or circle profile can carry.
    let solid = extrude_profile_exact(&profile(contour), Vec3::Z, 2.0, Tolerance::METRE)
        .expect("an L-shaped contour extrudes");

    let planes = solid
        .surfaces()
        .iter()
        .filter(|surface| matches!(surface, Surface::Plane(_)))
        .count();
    assert_eq!(planes, 8, "six walls plus two caps, got {planes}");

    let health = geometric_audit(&solid, Tolerance::METRE);
    assert!(
        health.is_consistent(),
        "contour solid must audit clean, found {:?}",
        health.defects()
    );
}

#[test]
fn an_arc_contour_extrudes_to_a_cylindrical_wall() {
    // A stadium/slot outline: two straight sides closed by two half-turns...
    // except a half turn is exactly what a single bulge cannot express, so
    // each end is built from two quarter arcs.
    let radius = 1.0;
    let right = Point2::new(2.0, 0.0);
    let left = Point2::new(-2.0, 0.0);
    let contour = Contour::new(vec![
        line(Point2::new(-2.0, -1.0), Point2::new(2.0, -1.0)),
        // right cap: -90 deg -> 0 -> +90 deg
        arc(right, radius, -std::f64::consts::FRAC_PI_2, 0.0),
        arc(right, radius, 0.0, std::f64::consts::FRAC_PI_2),
        line(Point2::new(2.0, 1.0), Point2::new(-2.0, 1.0)),
        // left cap: +90 deg -> 180 -> 270
        arc(
            left,
            radius,
            std::f64::consts::FRAC_PI_2,
            std::f64::consts::PI,
        ),
        arc(
            left,
            radius,
            std::f64::consts::PI,
            3.0 * std::f64::consts::FRAC_PI_2,
        ),
    ]);

    let solid = extrude_profile_exact(&profile(contour), Vec3::Z, 2.0, Tolerance::METRE)
        .expect("a stadium contour extrudes");

    let cylinders: Vec<f64> = solid
        .surfaces()
        .iter()
        .filter_map(|surface| match surface {
            Surface::Cylinder(cylinder) => Some(cylinder.radius),
            _ => None,
        })
        .collect();
    assert_eq!(
        cylinders.len(),
        4,
        "four quarter arcs must give four cylindrical walls, got {}",
        cylinders.len()
    );
    for found in &cylinders {
        assert!(
            (found - radius).abs() < 1e-12,
            "each wall must carry the contour's own radius, got {found}"
        );
    }

    let health = geometric_audit(&solid, Tolerance::METRE);
    assert!(
        health.is_consistent(),
        "arc contour solid must audit clean, found {:?}",
        health.defects()
    );
}

#[test]
fn an_inexact_curve_kind_is_refused_not_sampled() {
    // Tessellating an ellipse into chords would close, validate, and report a
    // plausible volume while not being the requested shape.
    let contour = Contour::new(vec![
        ProfileSegment {
            curve: Curve2::Ellipse(Ellipse2 {
                frame: Frame2 {
                    origin: Point2::new(0.0, 0.0),
                    x: Vec2::X,
                    y: Vec2::Y,
                },
                semi_axis_x: 2.0,
                semi_axis_y: 1.0,
            }),
            domain: Interval::new(0.0, std::f64::consts::PI),
            same_sense: true,
        },
        line(Point2::new(-2.0, 0.0), Point2::new(2.0, 0.0)),
    ]);
    let error = extrude_profile_exact(&profile(contour), Vec3::Z, 1.0, Tolerance::METRE)
        .expect_err("an elliptical segment is not exactly representable here");
    let text = format!("{error:?}");
    assert!(
        text.contains("elliptical"),
        "the refusal must name the curve kind, got {text}"
    );
}

#[test]
fn an_open_contour_is_refused_rather_than_closed_silently() {
    // Bridging the gap would change the profile the caller asked for.
    let contour = Contour::new(vec![
        line(Point2::new(0.0, 0.0), Point2::new(2.0, 0.0)),
        line(Point2::new(2.0, 0.0), Point2::new(2.0, 2.0)),
        line(Point2::new(2.0, 2.0), Point2::new(0.5, 2.0)),
    ]);
    let error = extrude_profile_exact(&profile(contour), Vec3::Z, 1.0, Tolerance::METRE)
        .expect_err("an open contour has no area");
    let text = format!("{error:?}");
    assert!(
        text.contains("does not close"),
        "the refusal must say the contour is open, got {text}"
    );
}

#[test]
fn a_contour_with_holes_now_builds_a_through_passage() {
    // ADR 0053 refused this; ADR 0058 implements it. The assertion flipped
    // from "must refuse" to a measured volume so the test states the new
    // capability instead of quietly disappearing.
    let outer = Contour::new(vec![
        line(Point2::new(0.0, 0.0), Point2::new(4.0, 0.0)),
        line(Point2::new(4.0, 0.0), Point2::new(4.0, 4.0)),
        line(Point2::new(4.0, 4.0), Point2::new(0.0, 4.0)),
        line(Point2::new(0.0, 4.0), Point2::new(0.0, 0.0)),
    ]);
    let hole = Contour::new(vec![
        line(Point2::new(1.0, 1.0), Point2::new(1.0, 2.0)),
        line(Point2::new(1.0, 2.0), Point2::new(2.0, 2.0)),
        line(Point2::new(2.0, 2.0), Point2::new(2.0, 1.0)),
        line(Point2::new(2.0, 1.0), Point2::new(1.0, 1.0)),
    ]);
    let value = Profile::Contour(ContourProfile {
        outer,
        holes: vec![hole],
    });
    let solid = extrude_profile_exact(&value, Vec3::Z, 1.0, Tolerance::METRE)
        .expect("a contour with holes extrudes");

    // 16 minus the 1x1 hole, times unit depth. A dropped hole reads 16.
    let volume = axiolid_measure::exact_properties(&solid, Tolerance::METRE)
        .map(|properties| properties.signed_volume)
        .expect("an all-planar prism is measurable");
    assert!((volume - 15.0).abs() < 1e-12, "expected 15, got {volume}");
}

/// A quarter arc on a LEFT-handed frame: y is -perp(x), so the parameter
/// runs clockwise in world orientation.
fn arc_left_handed(centre: Point2, radius: f64, from: f64, to: f64) -> ProfileSegment {
    ProfileSegment {
        curve: Curve2::Circle(Circle2 {
            frame: Frame2 {
                origin: centre,
                x: Vec2::X,
                y: Vec2::new(0.0, -1.0),
            },
            radius,
        }),
        domain: Interval::new(from, to),
        same_sense: true,
    }
}

#[test]
fn frame_handedness_decides_which_way_the_arc_turns() {
    // A left-handed frame runs the parameter backwards relative to the
    // plane, so the world sweep is the negation of the parameter sweep.
    // Ignoring that puts the arc on the wrong side of its chord -- the ring
    // still closes, so only the resulting geometry can tell.
    let radius = 1.0;
    let quarter = std::f64::consts::FRAC_PI_2;

    let centre_for = |left_handed: bool| {
        let segment = if left_handed {
            arc_left_handed(Point2::new(0.0, 0.0), radius, 0.0, quarter)
        } else {
            arc(Point2::new(0.0, 0.0), radius, 0.0, quarter)
        };
        // Close the quarter arc back to its start through the centre.
        let start = Point2::new(radius, 0.0);
        let end = if left_handed {
            Point2::new(0.0, -radius)
        } else {
            Point2::new(0.0, radius)
        };
        let contour = Contour::new(vec![
            segment,
            line(end, Point2::new(0.0, 0.0)),
            line(Point2::new(0.0, 0.0), start),
        ]);
        let solid = extrude_profile_exact(&profile(contour), Vec3::Z, 1.0, Tolerance::METRE)
            .expect("a quarter-pie contour extrudes");
        let health = geometric_audit(&solid, Tolerance::METRE);
        assert!(
            health.is_consistent(),
            "handedness {left_handed}: {:?}",
            health.defects()
        );
        solid
            .surfaces()
            .iter()
            .find_map(|surface| match surface {
                Surface::Cylinder(cylinder) => Some(cylinder.frame.origin),
                _ => None,
            })
            .expect("a cylindrical wall")
    };

    let right = centre_for(false);
    let left = centre_for(true);
    // Both arcs lie on the same circle centred at the origin, whichever way
    // the frame runs: that is the whole point of honouring handedness.
    assert!(
        right.x.hypot(right.y) < 1e-9,
        "right-handed arc centre should be the origin, got {right:?}"
    );
    assert!(
        left.x.hypot(left.y) < 1e-9,
        "left-handed arc centre should be the origin, got {left:?}"
    );
}

#[test]
fn a_half_turn_segment_is_refused_not_silently_bent() {
    // `bulge` is tan(sweep/4); at a half turn the tangent is 1 and at more
    // than that the chord no longer determines the arc. Accepting it would
    // silently build a different curve.
    let contour = Contour::new(vec![
        arc(Point2::new(0.0, 0.0), 1.0, 0.0, std::f64::consts::PI),
        line(Point2::new(-1.0, 0.0), Point2::new(1.0, 0.0)),
    ]);
    let error = extrude_profile_exact(&profile(contour), Vec3::Z, 1.0, Tolerance::METRE)
        .expect_err("a half turn cannot be one bulge");
    let text = format!("{error:?}");
    assert!(
        text.contains("half a turn"),
        "the refusal must name the half turn, got {text}"
    );
}
