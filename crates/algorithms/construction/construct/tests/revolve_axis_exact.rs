//! Exact full-turn revolution of a profile touching the axis (#379).
//!
//! An edge lying on the axis sweeps no face, an edge ending on it exactly
//! perpendicularly sweeps a disc cap, and any other straight edge ending on
//! it sweeps a cone up to its apex: a solid cylinder, a stepped shaft, a
//! cone, a bolt, a spinning top. Volumes are closed forms from the inputs
//! (`pi r^2 h`, a sum of those); every solid must audit clean. An arc ending
//! on the axis with its circle centred on it sweeps a sphere up to its pole
//! (a ball, a domed column); its centre is the sphere's to rounding, which
//! the quarter turns of an arc's parameter must not disturb. A profile
//! crossing the axis or touching it at an isolated vertex stays refused by
//! name, as does an edge ending on the axis within the tolerance of
//! perpendicular to it, and an arc whose circle is not centred on the axis (a
//! spindle or horn torus).

use axiolid_brep::ExactBRep;
use axiolid_brep_audit::geometric_audit;
use axiolid_construct::revolve_contour::revolve_arc_ring;
use axiolid_construct::revolve_exact::revolve_profile_exact;
use axiolid_contracts::{GeomError, GeomResult};
use axiolid_core::{Frame2, Interval, Point2, Point3, Tolerance, Vec2, Vec3};
use axiolid_curve::{Circle2, Curve2, Line2};
use axiolid_overlay::{ArcRing, ArcVertex};
use axiolid_profile::{Contour, ContourProfile, Profile, ProfileSegment, RectangleProfile};
use axiolid_surface::Surface;

const PI: f64 = std::f64::consts::PI;
const TAU: f64 = std::f64::consts::TAU;

fn polygon(points: &[(f64, f64)]) -> Profile {
    let points: Vec<Point2> = points.iter().map(|&(x, y)| Point2::new(x, y)).collect();
    let segments = (0..points.len())
        .map(|i| {
            let (from, to) = (points[i], points[(i + 1) % points.len()]);
            ProfileSegment {
                curve: Curve2::Line(Line2 {
                    origin: from,
                    direction: to - from,
                }),
                domain: Interval::UNIT,
                same_sense: true,
            }
        })
        .collect();
    Profile::Contour(ContourProfile {
        outer: Contour::new(segments),
        holes: Vec::new(),
    })
}

fn try_revolve(profile: &Profile, axis_x: f64) -> GeomResult<ExactBRep> {
    revolve_profile_exact(
        profile,
        Point3::new(axis_x, 0.0, 0.0),
        Vec3::Y,
        TAU,
        Tolerance::METRE,
    )
}

fn revolve(profile: &Profile, axis_x: f64) -> ExactBRep {
    let solid = try_revolve(profile, axis_x).expect("a profile touching the axis revolves");
    let health = geometric_audit(&solid, Tolerance::METRE);
    assert!(health.is_consistent(), "{:?}", health.defects());
    let topology = axiolid_topology::audit_brep(solid.topology());
    assert!(topology.is_closed_manifold(), "{topology:?}");
    solid
}

fn volume(solid: &ExactBRep) -> f64 {
    axiolid_measure::exact_properties(solid, Tolerance::METRE)
        .expect("measurable")
        .signed_volume
}

fn area(solid: &ExactBRep) -> f64 {
    axiolid_measure::exact_properties(solid, Tolerance::METRE)
        .expect("measurable")
        .area
}

fn close(got: f64, want: f64) {
    assert!(
        (got - want).abs() <= 1e-10 * want.abs(),
        "expected {want}, got {got}"
    );
}

fn count(solid: &ExactBRep, of: fn(&Surface) -> bool) -> usize {
    solid.surfaces().iter().filter(|s| of(s)).count()
}

fn named(result: GeomResult<ExactBRep>, wanted: &str) {
    match result {
        Err(GeomError::UnsupportedInput { input, .. }) => {
            assert!(input.contains(wanted), "got {input:?}, wanted {wanted:?}");
        }
        other => panic!("expected a refusal naming {wanted:?}, got {other:?}"),
    }
}

#[test]
fn a_solid_cylinder_has_the_volume_pi_r_squared_h() {
    let (r, h) = (1.5, 4.0);
    let solid = revolve(&polygon(&[(0.0, 0.0), (r, 0.0), (r, h), (0.0, h)]), 0.0);
    close(volume(&solid), PI * r * r * h);
    assert_eq!(count(&solid, |s| matches!(s, Surface::Cylinder(_))), 1);
    assert_eq!(count(&solid, |s| matches!(s, Surface::Plane(_))), 2);
}

#[test]
fn a_cylinder_wound_clockwise_or_on_the_negative_side_is_the_same_solid() {
    let (r, h) = (1.5, 4.0);
    let clockwise = revolve(&polygon(&[(0.0, h), (r, h), (r, 0.0), (0.0, 0.0)]), 0.0);
    close(volume(&clockwise), PI * r * r * h);
    let mirrored = revolve(&polygon(&[(0.0, 0.0), (0.0, h), (-r, h), (-r, 0.0)]), 0.0);
    close(volume(&mirrored), PI * r * r * h);
}

#[test]
fn a_ring_on_the_negative_side_of_the_axis_is_outward() {
    // x in [-3, -1] about the axis: the tube of the mirrored side.
    let solid = revolve(
        &polygon(&[(-3.0, 0.0), (-1.0, 0.0), (-1.0, 2.0), (-3.0, 2.0)]),
        0.0,
    );
    close(volume(&solid), PI * (9.0 - 1.0) * 2.0);
}

#[test]
fn a_rectangle_with_an_edge_on_the_offset_axis_is_a_solid_cylinder() {
    // The rectangle spans x in [-1, 1]; the axis through x = 1 runs along its
    // right edge, so the cylinder has radius 2 and height 3.
    let profile = Profile::Rectangle(RectangleProfile {
        x: 2.0,
        y: 3.0,
        thickness: None,
        outer_radius: None,
        inner_radius: None,
    });
    let solid = revolve(&profile, 1.0);
    close(volume(&solid), PI * 4.0 * 3.0);
}

#[test]
fn a_stepped_shaft_sums_its_cylinders() {
    // Radius 2 over [0, 1], radius 1 over [1, 3]: the axis edge runs the
    // whole length and two discs cap it.
    let solid = revolve(
        &polygon(&[
            (0.0, 0.0),
            (2.0, 0.0),
            (2.0, 1.0),
            (1.0, 1.0),
            (1.0, 3.0),
            (0.0, 3.0),
        ]),
        0.0,
    );
    close(volume(&solid), PI * (4.0 * 1.0 + 1.0 * 2.0));
    assert_eq!(count(&solid, |s| matches!(s, Surface::Cylinder(_))), 2);
    // Two discs and the annular shoulder.
    assert_eq!(count(&solid, |s| matches!(s, Surface::Plane(_))), 3);
}

#[test]
fn a_profile_crossing_the_axis_is_refused_by_name() {
    named(
        try_revolve(
            &polygon(&[(-1.0, 0.0), (1.0, 0.0), (1.0, 2.0), (-1.0, 2.0)]),
            0.0,
        ),
        "crossing the axis",
    );
}

#[test]
fn a_vertex_touching_the_axis_alone_is_refused_by_name() {
    // A diamond whose left corner is the only contact: an isolated vertex,
    // the pinch point of two cones.
    named(
        try_revolve(
            &polygon(&[(1.0, 0.0), (2.0, 1.0), (1.0, 2.0), (0.0, 1.0)]),
            0.0,
        ),
        "isolated vertex",
    );
}

/// The audited solid of a cone profile, with its volume and area checked.
fn solid_of(
    what: &str,
    points: &[(f64, f64)],
    axis_x: f64,
    volume_want: f64,
    area_want: f64,
) -> ExactBRep {
    let solid = revolve(&polygon(points), axis_x);
    let (v, a) = (volume(&solid), area(&solid));
    assert!(
        (v - volume_want).abs() <= 1e-10 * volume_want
            && (a - area_want).abs() <= 1e-10 * area_want,
        "{what}: volume {v} want {volume_want}, area {a} want {area_want}"
    );
    solid
}

#[test]
fn a_cone_with_its_apex_up_is_a_third_of_its_cylinder() {
    // Radius 1.5, height 2: V = pi r^2 h / 3, A = pi r (r + slant), slant 2.5.
    let solid = solid_of(
        "cone",
        &[(0.0, 0.0), (1.5, 0.0), (0.0, 2.0)],
        0.0,
        4.71238898038469,
        18.84955592153876,
    );
    close(volume(&solid), 1.5 * PI);
    assert_eq!(count(&solid, |s| matches!(s, Surface::Cone(_))), 1);
    assert_eq!(count(&solid, |s| matches!(s, Surface::Plane(_))), 1);
    // The cone is framed at its apex with radius exactly zero, so the apex
    // is at v = 0 and a consumer closing the pole gap meets no rounding.
    let cone = solid
        .surfaces()
        .iter()
        .find_map(|surface| match surface {
            Surface::Cone(cone) => Some(cone),
            _ => None,
        })
        .expect("a cone");
    assert_eq!((cone.radius, cone.frame.origin.y), (0.0, 2.0));
}

#[test]
fn a_cone_with_its_apex_down_is_the_same_solid() {
    let solid = solid_of(
        "cone",
        &[(0.0, 0.0), (1.5, 2.0), (0.0, 2.0)],
        0.0,
        4.71238898038469,
        18.84955592153876,
    );
    assert_eq!(count(&solid, |s| matches!(s, Surface::Cone(_))), 1);
    assert_eq!(count(&solid, |s| matches!(s, Surface::Plane(_))), 1);
}

#[test]
fn a_cone_on_the_negative_side_is_mirrored() {
    solid_of(
        "mirrored cone",
        &[(0.0, 0.0), (0.0, 2.0), (-1.5, 0.0)],
        0.0,
        1.5 * PI,
        18.84955592153876,
    );
}

#[test]
fn a_cone_about_an_offset_axis_is_the_same_cone() {
    solid_of(
        "offset cone",
        &[(2.5, 0.0), (4.0, 0.0), (2.5, 2.0)],
        2.5,
        1.5 * PI,
        18.84955592153876,
    );
}

#[test]
fn a_bolt_is_a_cylinder_with_a_conical_point() {
    // Radius 0.5 over [0, 3] with a point reaching down to -0.4.
    let solid = revolve(
        &polygon(&[(0.0, -0.4), (0.5, 0.0), (0.5, 3.0), (0.0, 3.0)]),
        0.0,
    );
    close(volume(&solid), 2.4609142453120048);
    assert_eq!(count(&solid, |s| matches!(s, Surface::Cone(_))), 1);
    assert_eq!(count(&solid, |s| matches!(s, Surface::Cylinder(_))), 1);
}

#[test]
fn a_spinning_top_has_two_cones_and_two_apexes() {
    let solid = solid_of(
        "spinning top",
        &[(0.0, 0.0), (1.0, 1.0), (0.0, 2.0)],
        0.0,
        2.0943951023931953,
        8.885765876316732,
    );
    assert_eq!(count(&solid, |s| matches!(s, Surface::Cone(_))), 2);
}

#[test]
fn a_frustum_ending_in_a_cone_has_a_disc_and_two_cones() {
    solid_of(
        "frustum to apex",
        &[(0.0, 0.0), (1.0, 0.0), (0.5, 1.0), (0.0, 2.0)],
        0.0,
        2.0943951023931953,
        10.166407384630519,
    );
}

#[test]
fn a_cavity_inside_a_spinning_top_is_a_void_not_a_second_solid() {
    // The top's range along the axis comes from its apexes, off every
    // station off the axis: the cavity nests only when it counts them.
    let solid = revolve(
        &polygon(&[
            (0.0, 0.0),
            (1.0, 1.0),
            (0.0, 2.0),
            (0.0, 1.5),
            (0.2, 1.5),
            (0.2, 0.5),
            (0.0, 0.5),
        ]),
        0.0,
    );
    close(volume(&solid), 1.9687313962496034);
    assert_eq!(solid.topology().solids().len(), 1);
    assert_eq!(solid.topology().solids()[0].voids.len(), 1);
    // The top is the solid and the cavity its void, not the other way round:
    // the outer shell holds both cones.
    let topology = solid.topology();
    let outer = &topology.shells()[topology.solids()[0].outer.index()];
    let cones = outer
        .faces
        .iter()
        .filter(|(face, _)| {
            let surface = topology.faces()[face.index()].surface.expect("a surface");
            matches!(solid.surfaces()[surface.index()], Surface::Cone(_))
        })
        .count();
    assert_eq!(cones, 2);
}

#[test]
fn a_vertex_within_tolerance_of_the_axis_but_not_on_it_is_refused() {
    named(
        try_revolve(
            &polygon(&[(1e-9, 0.0), (1.0, 0.0), (1.0, 1.0), (1e-9, 1.0)]),
            0.0,
        ),
        "within tolerance of the axis",
    );
}

#[test]
fn an_edge_a_hair_off_perpendicular_ending_on_the_axis_is_refused() {
    // The top edge rises by 1e-12 towards the axis: not a disc, exactly, and
    // not a cone to the tolerance.
    named(
        try_revolve(
            &polygon(&[(0.0, 0.0), (1.0, 0.0), (1.0, 2.0), (0.0, 2.0 + 1e-12)]),
            0.0,
        ),
        "within the tolerance of perpendicular",
    );
}

#[test]
fn an_arc_tangent_to_the_axis_is_refused_as_a_horn_torus() {
    // A quarter arc of the circle about (1, 0) starting at (0, 0) at angle pi,
    // where its x is exactly 0, then a straight run back along the top and
    // the axis: the arc's start is on the axis, but its circle touches the
    // axis there instead of being centred on it.
    let line = |from: Point2, to: Point2| ProfileSegment {
        curve: Curve2::Line(Line2 {
            origin: from,
            direction: to - from,
        }),
        domain: Interval::UNIT,
        same_sense: true,
    };
    let arc = ProfileSegment {
        curve: Curve2::Circle(Circle2 {
            frame: Frame2 {
                origin: Point2::new(1.0, 0.0),
                x: Vec2::X,
                y: Vec2::Y,
            },
            radius: 1.0,
        }),
        domain: Interval::new(PI, 1.5 * PI),
        same_sense: true,
    };
    let start = Point2::new(1.0 + PI.cos(), PI.sin());
    let end = Point2::new(1.0 + (1.5 * PI).cos(), (1.5 * PI).sin());
    assert_eq!(start.x, 0.0);
    let top = Point2::new(1.0, 1.0);
    let axis_top = Point2::new(0.0, 1.0);
    named(
        try_revolve(
            &Profile::Contour(ContourProfile {
                outer: Contour::new(vec![
                    arc,
                    line(end, top),
                    line(top, axis_top),
                    line(axis_top, start),
                ]),
                holes: Vec::new(),
            }),
            0.0,
        ),
        "spindle or horn torus",
    );
}

#[test]
fn an_arc_across_the_top_of_a_cylinder_is_not_a_disc() {
    // The top edge runs from (1, 1) to the axis at (0, 1), at the same
    // height, yet it is an arc, not a perpendicular line.
    let vertex = |x: f64, y: f64, bulge: f64| ArcVertex {
        point: Point2::new(x, y),
        bulge,
    };
    let ring = ArcRing::new(vec![
        vertex(0.0, 0.0, 0.0),
        vertex(1.0, 0.0, 0.0),
        vertex(1.0, 1.0, 0.1),
        vertex(0.0, 1.0, 0.0),
    ]);
    named(
        revolve_arc_ring(&ring, Point3::ZERO, Tolerance::METRE),
        "spindle or horn torus",
    );
}

fn cylinder_polygon() -> Profile {
    polygon(&[(0.0, 0.0), (1.0, 0.0), (1.0, 4.0), (0.0, 4.0)])
}

#[test]
fn a_tilted_axis_within_the_tolerance_is_refused_for_a_touching_section() {
    // (0, 4) is 4e-3 off the real axis, which the sign against the vertical
    // line through the origin cannot see.
    for tilt in [1e-12, 1e-3, 1.4e-3] {
        named(
            revolve_profile_exact(
                &cylinder_polygon(),
                Point3::ZERO,
                Vec3::new(tilt, 1.0, 0.0),
                TAU,
                Tolerance::METRE,
            ),
            "not exactly the profile's y axis",
        );
    }
    named(
        revolve_profile_exact(
            &cylinder_polygon(),
            Point3::ZERO,
            Vec3::new(0.0, 1.0, 1e-3),
            TAU,
            Tolerance::METRE,
        ),
        "not exactly the profile's y axis",
    );
}

#[test]
fn an_axis_off_the_plane_by_less_than_the_tolerance_is_refused_for_a_touching_section() {
    named(
        revolve_profile_exact(
            &cylinder_polygon(),
            Point3::new(0.0, 0.0, 1e-7),
            Vec3::Y,
            TAU,
            Tolerance::METRE,
        ),
        "not exactly the profile's y axis",
    );
}

#[test]
fn an_axis_pointing_down_is_exact_too() {
    let solid = revolve_profile_exact(
        &cylinder_polygon(),
        Point3::ZERO,
        Vec3::NEG_Y,
        TAU,
        Tolerance::METRE,
    )
    .expect("a downward axis is exactly y");
    close(volume(&solid), PI * 4.0);
}

#[test]
fn a_cavity_along_the_axis_is_a_void_shell() {
    // A cylinder of radius 2 and height 3 with a disc cavity of radius 1
    // between heights 1 and 2: two contacts with the axis.
    let solid = revolve(
        &polygon(&[
            (0.0, 0.0),
            (2.0, 0.0),
            (2.0, 3.0),
            (0.0, 3.0),
            (0.0, 2.0),
            (1.0, 2.0),
            (1.0, 1.0),
            (0.0, 1.0),
        ]),
        0.0,
    );
    // The same body from its parts: the outer cylinder, and the hole
    // revolved on its own, joined as a void as `revolve_region` does.
    let outer = revolve(
        &polygon(&[(0.0, 0.0), (2.0, 0.0), (2.0, 3.0), (0.0, 3.0)]),
        0.0,
    );
    let hole = revolve(
        &polygon(&[(0.0, 1.0), (1.0, 1.0), (1.0, 2.0), (0.0, 2.0)]),
        0.0,
    );
    close(volume(&solid), volume(&outer) - volume(&hole));
    close(volume(&solid), PI * (4.0 * 3.0 - 1.0));
    assert_eq!(solid.topology().solids().len(), 1);
    assert_eq!(solid.topology().solids()[0].voids.len(), 1);
    assert_eq!(
        solid.topology().faces().len(),
        outer.topology().faces().len() + hole.topology().faces().len()
    );
}

#[test]
fn two_cylinders_joined_along_the_axis_are_two_solids() {
    // Radius 1 over [0, 1] and radius 2 over [2, 3], joined by an axis edge
    // of no width: two separate solids.
    let solid = revolve(
        &polygon(&[
            (0.0, 0.0),
            (1.0, 0.0),
            (1.0, 1.0),
            (0.0, 1.0),
            (0.0, 2.0),
            (2.0, 2.0),
            (2.0, 3.0),
            (0.0, 3.0),
        ]),
        0.0,
    );
    assert_eq!(solid.topology().solids().len(), 2);
    close(volume(&solid), PI * (1.0 + 4.0));
}

#[test]
fn a_circle_tangent_to_the_axis_is_refused_as_a_horn_torus() {
    named(
        try_revolve(
            &Profile::Circle(axiolid_profile::CircleProfile {
                radius: 1.0,
                thickness: None,
            }),
            -1.0,
        ),
        "horn torus",
    );
}

fn line_segment(from: Point2, to: Point2) -> ProfileSegment {
    ProfileSegment {
        curve: Curve2::Line(Line2 {
            origin: from,
            direction: to - from,
        }),
        domain: Interval::UNIT,
        same_sense: true,
    }
}

fn arc_segment(centre: Point2, radius: f64, from: f64, to: f64) -> ProfileSegment {
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

fn contour_of(segments: Vec<ProfileSegment>) -> Profile {
    Profile::Contour(ContourProfile {
        outer: Contour::new(segments),
        holes: Vec::new(),
    })
}

/// The half disc of radius `r` about the axis at `axis_x`: a quarter turn
/// of the arc's parameter either side of its centre, closed along the axis.
fn half_disc(radius: f64, axis_x: f64) -> Profile {
    contour_of(vec![
        arc_segment(Point2::new(axis_x, 0.0), radius, -PI / 2.0, PI / 2.0),
        line_segment(Point2::new(axis_x, radius), Point2::new(axis_x, -radius)),
    ])
}

#[test]
fn a_half_disc_revolves_to_a_ball_at_any_radius_and_offset() {
    // The arc's ends are `cos` and `sin` of a quarter turn: exactly on the
    // axis only because the lowering evaluates a quadrant exactly.
    for (radius, axis_x) in [
        (1.0, 0.0),
        (0.37, 0.0),
        (1000.0, 0.0),
        (0.37, 12.5),
        (2.0, -3.0),
    ] {
        let ball = revolve(&half_disc(radius, axis_x), axis_x);
        let want = 4.0 / 3.0 * PI * radius.powi(3);
        let got = volume(&ball);
        assert!(
            (got - want).abs() <= 1e-10 * want,
            "radius {radius} at {axis_x}: volume {got}, want {want}"
        );
        let want = 4.0 * PI * radius * radius;
        let got = area(&ball);
        assert!(
            (got - want).abs() <= 1e-10 * want,
            "radius {radius} at {axis_x}: area {got}, want {want}"
        );
        assert!(ball
            .surfaces()
            .iter()
            .all(|s| matches!(s, Surface::Sphere(_))));
    }
}

#[test]
fn a_hemisphere_on_a_cylinder_is_a_domed_column() {
    let (r, h) = (0.3, 1.2);
    let column = contour_of(vec![
        line_segment(Point2::new(0.0, 0.0), Point2::new(r, 0.0)),
        line_segment(Point2::new(r, 0.0), Point2::new(r, h)),
        arc_segment(Point2::new(0.0, h), r, 0.0, PI / 2.0),
        line_segment(Point2::new(0.0, h + r), Point2::new(0.0, 0.0)),
    ]);
    let solid = revolve(&column, 0.0);
    close(volume(&solid), 0.3958406743523139);
    close(area(&solid), 3.1101767270538954);
    assert_eq!(count(&solid, |s| matches!(s, Surface::Sphere(_))), 1);
    assert_eq!(count(&solid, |s| matches!(s, Surface::Cylinder(_))), 1);
    assert_eq!(count(&solid, |s| matches!(s, Surface::Plane(_))), 1);
}

#[test]
fn a_ball_with_a_cavity_along_the_axis_has_one_void() {
    let ball = contour_of(vec![
        arc_segment(Point2::ZERO, 2.0, -PI / 2.0, PI / 2.0),
        line_segment(Point2::new(0.0, 2.0), Point2::new(0.0, 0.5)),
        line_segment(Point2::new(0.0, 0.5), Point2::new(0.5, 0.5)),
        line_segment(Point2::new(0.5, 0.5), Point2::new(0.5, -0.5)),
        line_segment(Point2::new(0.5, -0.5), Point2::new(0.0, -0.5)),
        line_segment(Point2::new(0.0, -0.5), Point2::new(0.0, -2.0)),
    ]);
    let solid = revolve(&ball, 0.0);
    close(volume(&solid), 32.72492347489368);
    assert_eq!(solid.topology().solids().len(), 1);
    assert_eq!(solid.topology().solids()[0].voids.len(), 1);
}

/// `(0, 0) -> (1, 1)` by an arc of the given bulge, then back along the top
/// and the axis.
fn arc_from_the_axis(bulge: f64) -> ArcRing {
    ArcRing::new(vec![
        ArcVertex {
            point: Point2::new(0.0, 0.0),
            bulge,
        },
        ArcVertex {
            point: Point2::new(1.0, 1.0),
            bulge: 0.0,
        },
        ArcVertex {
            point: Point2::new(0.0, 1.0),
            bulge: 0.0,
        },
    ])
}

#[test]
fn a_ring_with_a_bulge_that_is_a_quarter_circle_is_a_lower_hemisphere() {
    // `tan(pi / 8)` is the bulge of a quarter arc, from (0, 0) to (1, 1)
    // about (0, 1): a bulge cannot encode it exactly, so its computed centre
    // is off the axis by rounding, which the gate absorbs.
    let solid = revolve_arc_ring(
        &arc_from_the_axis((PI / 8.0).tan()),
        Point3::ZERO,
        Tolerance::METRE,
    )
    .expect("a lower hemisphere");
    let health = geometric_audit(&solid, Tolerance::METRE);
    assert!(health.is_consistent(), "{:?}", health.defects());
    close(volume(&solid), 2.0943951023931953);
    close(area(&solid), 9.42477796076938);
}

#[test]
fn an_arc_whose_circle_is_not_centred_on_the_axis_is_a_spindle_or_horn_torus() {
    named(
        revolve_arc_ring(&arc_from_the_axis(0.2), Point3::ZERO, Tolerance::METRE),
        "spindle or horn torus",
    );
}

#[test]
fn an_arc_running_across_the_axis_is_refused() {
    named(
        revolve_arc_ring(
            &arc_from_the_axis(-(1.0 + 2.0_f64.sqrt())),
            Point3::ZERO,
            Tolerance::METRE,
        ),
        "crossing the axis",
    );
}

#[test]
fn an_arc_with_both_ends_on_the_axis_is_refused() {
    let ring = ArcRing::new(vec![
        ArcVertex {
            point: Point2::new(0.0, -1.0),
            bulge: 0.3,
        },
        ArcVertex {
            point: Point2::new(0.0, 1.0),
            bulge: 0.0,
        },
        ArcVertex {
            point: Point2::new(0.0, 0.0),
            bulge: 0.0,
        },
    ]);
    named(
        revolve_arc_ring(&ring, Point3::ZERO, Tolerance::METRE),
        "both ends on the axis",
    );
}

#[test]
fn an_arc_ending_on_the_axis_with_no_axis_edge_is_a_pinch_point() {
    // The arc's pole has no edge along the axis: two cones' worth of
    // surface meeting at a point.
    let ring = ArcRing::new(vec![
        ArcVertex {
            point: Point2::new(0.0, 0.0),
            bulge: (PI / 8.0).tan(),
        },
        ArcVertex {
            point: Point2::new(1.0, 1.0),
            bulge: 0.0,
        },
        ArcVertex {
            point: Point2::new(2.0, 0.0),
            bulge: 0.0,
        },
        ArcVertex {
            point: Point2::new(1.0, -1.0),
            bulge: 0.0,
        },
    ]);
    named(
        revolve_arc_ring(&ring, Point3::ZERO, Tolerance::METRE),
        "isolated vertex where an arc ends",
    );
}

#[test]
fn a_solid_inside_a_cavity_along_the_axis_is_refused_by_name() {
    // Three nested pieces: the outer cylinder, its cavity and an island in it.
    named(
        try_revolve(
            &polygon(&[
                (0.0, 0.0),
                (3.0, 0.0),
                (3.0, 5.0),
                (0.0, 5.0),
                (0.0, 4.0),
                (2.0, 4.0),
                (2.0, 1.0),
                (0.0, 1.0),
                (0.0, 2.0),
                (1.0, 2.0),
                (1.0, 3.0),
                (0.0, 3.0),
            ]),
            0.0,
        ),
        "nested inside a cavity",
    );
}
