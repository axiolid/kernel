//! Exact centre-line extrusion (ADR 0056).
//!
//! Areas are checked against closed forms computed independently, not
//! against another run of the offsetter.

use axiolid_brep_audit::geometric_audit;
use axiolid_construct::center_line_exact::center_line_contour;
use axiolid_construct::extrude::extrude_profile_exact;
use axiolid_core::{Frame2, Interval, Point2, Tolerance, Vec2, Vec3};
use axiolid_curve::{Circle2, Curve2, Ellipse2, Line2};
use axiolid_profile::{CenterLineProfile, Contour, Profile, ProfileSegment};
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

fn center_line(segments: Vec<ProfileSegment>, width: f64) -> Profile {
    Profile::CenterLine(CenterLineProfile::from_width(Contour::new(segments), width))
}

#[test]
fn a_straight_centre_line_extrudes_to_width_times_length() {
    let (length, width, depth) = (5.0, 0.4, 2.0);
    let profile = center_line(
        vec![line(Point2::new(0.0, 0.0), Point2::new(length, 0.0))],
        width,
    );
    let solid = extrude_profile_exact(&profile, Vec3::Z, depth, Tolerance::METRE)
        .expect("a straight centre line extrudes");

    let properties =
        axiolid_measure::exact_properties(&solid, Tolerance::METRE).expect("all faces planar");
    let expected = length * width * depth;
    assert!(
        (properties.signed_volume - expected).abs() < 1e-12,
        "expected {expected}, got {}",
        properties.signed_volume
    );
    assert!(geometric_audit(&solid, Tolerance::METRE).is_consistent());
}

#[test]
fn a_curved_centre_line_keeps_exact_arc_walls() {
    // A quarter-arc path of radius 2, width 0.4. The offsets are CONCENTRIC
    // arcs of radius 1.8 and 2.2, so the solid must carry two genuine
    // cylindrical walls -- a flattened offset would emit a fan of planes.
    let (radius, width) = (2.0, 0.4);
    let quarter = core::f64::consts::FRAC_PI_2;
    let profile = center_line(
        vec![arc(Point2::new(0.0, 0.0), radius, 0.0, quarter)],
        width,
    );
    let solid = extrude_profile_exact(&profile, Vec3::Z, 1.0, Tolerance::METRE)
        .expect("a curved centre line extrudes");

    let mut radii: Vec<f64> = solid
        .surfaces()
        .iter()
        .filter_map(|s| match s {
            Surface::Cylinder(c) => Some(c.radius),
            _ => None,
        })
        .collect();
    radii.sort_by(f64::total_cmp);
    assert_eq!(radii.len(), 2, "two concentric walls, got {radii:?}");
    assert!(
        (radii[0] - (radius - width / 2.0)).abs() < 1e-12,
        "inner wall should be 1.8, got {}",
        radii[0]
    );
    assert!(
        (radii[1] - (radius + width / 2.0)).abs() < 1e-12,
        "outer wall should be 2.2, got {}",
        radii[1]
    );

    let health = geometric_audit(&solid, Tolerance::METRE);
    assert!(
        health.is_consistent(),
        "offset arcs must agree with their pcurves, found {:?}",
        health.defects()
    );
}

#[test]
fn the_offset_side_follows_the_turn_direction() {
    // Measured: the left normal of a counter-clockwise arc points INWARD, so
    // the left offset shrinks the radius. Reversing the sweep must swap
    // which side grows -- but both radii must still be 1.8 and 2.2.
    let (radius, width) = (2.0, 0.4);
    let quarter = core::f64::consts::FRAC_PI_2;
    for (from, to) in [(0.0, quarter), (quarter, 0.0)] {
        let profile = center_line(vec![arc(Point2::new(0.0, 0.0), radius, from, to)], width);
        let solid = extrude_profile_exact(&profile, Vec3::Z, 1.0, Tolerance::METRE)
            .expect("either traversal extrudes");
        let mut radii: Vec<f64> = solid
            .surfaces()
            .iter()
            .filter_map(|s| match s {
                Surface::Cylinder(c) => Some(c.radius),
                _ => None,
            })
            .collect();
        radii.sort_by(f64::total_cmp);
        assert_eq!(radii.len(), 2, "sweep {from}->{to}: got {radii:?}");
        assert!(
            (radii[0] - 1.8).abs() < 1e-12 && (radii[1] - 2.2).abs() < 1e-12,
            "sweep {from}->{to} must still straddle the path, got {radii:?}"
        );
    }
}

#[test]
fn an_elliptical_path_segment_is_refused_not_sampled() {
    // Measured: offsetting a 3:1 ellipse by 0.3 gives a curve whose best-fit
    // ellipse residual is 0.074 -- it is simply not an ellipse. Emitting one
    // anyway, or chord-sampling it, would be a silent approximation.
    let profile = center_line(
        vec![ProfileSegment {
            curve: Curve2::Ellipse(Ellipse2 {
                frame: Frame2 {
                    origin: Point2::new(0.0, 0.0),
                    x: Vec2::X,
                    y: Vec2::Y,
                },
                semi_axis_x: 3.0,
                semi_axis_y: 1.0,
            }),
            domain: Interval::new(0.0, 1.0),
            same_sense: true,
        }],
        0.4,
    );
    let error = extrude_profile_exact(&profile, Vec3::Z, 1.0, Tolerance::METRE)
        .expect_err("an ellipse has no elliptical offset");
    assert!(
        format!("{error:?}").contains("elliptical"),
        "the refusal must name the curve kind, got {error:?}"
    );
}

#[test]
fn a_half_width_that_collapses_an_arc_is_refused() {
    // Half-width 2.5 against a path arc of radius 2: the inner offset would
    // have radius -0.5, which is not a curve. Emitting it would produce a
    // ring that turns itself inside out at the corner.
    let profile = center_line(
        vec![arc(
            Point2::new(0.0, 0.0),
            2.0,
            0.0,
            core::f64::consts::FRAC_PI_2,
        )],
        5.0,
    );
    let error = extrude_profile_exact(&profile, Vec3::Z, 1.0, Tolerance::METRE)
        .expect_err("the width swallows the arc");
    assert!(format!("{error:?}").contains("collapses"), "got {error:?}");
}

#[test]
fn a_disconnected_path_is_refused() {
    let profile = center_line(
        vec![
            line(Point2::new(0.0, 0.0), Point2::new(2.0, 0.0)),
            line(Point2::new(3.0, 0.0), Point2::new(5.0, 0.0)),
        ],
        0.4,
    );
    let error = center_line_contour(
        match &profile {
            Profile::CenterLine(c) => c,
            _ => unreachable!(),
        },
        Tolerance::METRE,
    )
    .expect_err("a gap is not a path");
    assert!(
        format!("{error:?}").contains("disconnected"),
        "got {error:?}"
    );
}

#[test]
fn a_closed_path_is_refused_because_it_denotes_an_annulus() {
    // A closed centre line encloses a hole; a single outer ring cannot say
    // that, and welding the ends would silently fill it.
    let profile = center_line(
        vec![
            line(Point2::new(0.0, 0.0), Point2::new(2.0, 0.0)),
            line(Point2::new(2.0, 0.0), Point2::new(2.0, 2.0)),
            line(Point2::new(2.0, 2.0), Point2::new(0.0, 2.0)),
            line(Point2::new(0.0, 2.0), Point2::new(0.0, 0.0)),
        ],
        0.2,
    );
    let error = extrude_profile_exact(&profile, Vec3::Z, 1.0, Tolerance::METRE)
        .expect_err("a closed centre line is an annulus");
    assert!(format!("{error:?}").contains("annulus"), "got {error:?}");
}

#[test]
fn a_non_positive_half_width_is_refused() {
    for width in [0.0, -1.0, f64::NAN] {
        let profile = center_line(
            vec![line(Point2::new(0.0, 0.0), Point2::new(2.0, 0.0))],
            width,
        );
        assert!(
            extrude_profile_exact(&profile, Vec3::Z, 1.0, Tolerance::METRE).is_err(),
            "width {width} must be refused"
        );
    }
}

#[test]
fn a_left_handed_arc_frame_still_offsets_to_the_right_side() {
    // `Circle2` stores x and y independently, so the frame can be
    // left-handed: the parameter then runs clockwise in world orientation.
    // The offset side depends on the WORLD turn, so dropping handedness
    // swaps which side grows -- and the ring still closes, so only the
    // resulting radii reveal it.
    let (radius, width) = (2.0, 0.4);
    let quarter = core::f64::consts::FRAC_PI_2;
    let profile = center_line(
        vec![ProfileSegment {
            curve: Curve2::Circle(Circle2 {
                frame: Frame2 {
                    origin: Point2::new(0.0, 0.0),
                    x: Vec2::X,
                    // y = -perp(x): a left-handed frame.
                    y: Vec2::new(0.0, -1.0),
                },
                radius,
            }),
            domain: Interval::new(0.0, quarter),
            same_sense: true,
        }],
        width,
    );
    let solid = extrude_profile_exact(&profile, Vec3::Z, 1.0, Tolerance::METRE)
        .expect("a left-handed arc path extrudes");

    let mut radii: Vec<f64> = solid
        .surfaces()
        .iter()
        .filter_map(|s| match s {
            Surface::Cylinder(c) => Some(c.radius),
            _ => None,
        })
        .collect();
    radii.sort_by(f64::total_cmp);
    assert_eq!(radii.len(), 2, "two walls, got {radii:?}");
    // The offsets must straddle the path radius whichever way the frame runs.
    assert!(
        (radii[0] - (radius - width / 2.0)).abs() < 1e-12
            && (radii[1] - (radius + width / 2.0)).abs() < 1e-12,
        "expected 1.8 and 2.2, got {radii:?}"
    );
    assert!(geometric_audit(&solid, Tolerance::METRE).is_consistent());
}
