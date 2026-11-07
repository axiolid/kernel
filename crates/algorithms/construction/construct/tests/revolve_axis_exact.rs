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
//!
//! The axis may lie anywhere in the profile plane (#450): it is built exactly
//! as given, in the axis's own frame, never snapped onto +y. An axis with a
//! component out of the plane is refused by name.

use axiolid_brep::ExactBRep;
use axiolid_brep_audit::geometric_audit;
use axiolid_construct::revolve_contour::revolve_arc_ring;
use axiolid_construct::revolve_exact::revolve_profile_exact;
use axiolid_contracts::{GeomError, GeomResult};
use axiolid_core::{Frame2, Interval, Point2, Point3, Tolerance, Vec2, Vec3};
use axiolid_curve::{Circle2, Curve2, Curve3, Line2};
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
fn a_tilted_axis_through_a_vertex_puts_the_section_across_it() {
    // (0, 4) is left of the axis tilted by `tilt` while (1, 0) is right of it:
    // the section crosses the axis it is revolved about, as given (it used
    // to be snapped onto y and accepted).
    for tilt in [1e-12, 1e-3, 1.4e-3] {
        named(
            revolve_profile_exact(
                &cylinder_polygon(),
                Point3::ZERO,
                Vec3::new(tilt, 1.0, 0.0),
                TAU,
                Tolerance::METRE,
            ),
            "crossing the axis",
        );
    }
}

#[test]
fn an_axis_out_of_the_profile_plane_is_refused_by_name() {
    named(
        revolve_profile_exact(
            &cylinder_polygon(),
            Point3::ZERO,
            Vec3::new(0.0, 1.0, 1e-3),
            TAU,
            Tolerance::METRE,
        ),
        "out of the profile plane",
    );
    named(
        revolve_profile_exact(
            &cylinder_polygon(),
            Point3::new(0.0, 0.0, 1e-7),
            Vec3::Y,
            TAU,
            Tolerance::METRE,
        ),
        "out of the profile plane",
    );
    // A partial turn and the public ring entry refuse it too.
    named(
        revolve_profile_exact(
            &cylinder_polygon(),
            Point3::new(0.0, 0.0, 1e-7),
            Vec3::Y,
            1.0,
            Tolerance::METRE,
        ),
        "out of the profile plane",
    );
    named(
        revolve_arc_ring(
            &ArcRing::from_points(&[
                Point2::new(1.0, 0.0),
                Point2::new(2.0, 0.0),
                Point2::new(2.0, 1.0),
            ]),
            Point3::new(0.0, 0.0, 1e-7),
            Tolerance::METRE,
        ),
        "out of the profile plane",
    );
}

#[test]
fn a_zero_or_non_finite_axis_is_invalid_input() {
    for direction in [Vec3::ZERO, Vec3::new(f64::NAN, 1.0, 0.0)] {
        match revolve_profile_exact(
            &cylinder_polygon(),
            Point3::ZERO,
            direction,
            TAU,
            Tolerance::METRE,
        ) {
            Err(GeomError::InvalidInput(text)) => assert!(text.contains("finite non-zero")),
            other => panic!("expected invalid input, got {other:?}"),
        }
    }
}

/// `profile` revolved through `angle` about the in-plane axis through
/// `origin` along `direction`, audited.
fn revolve_about(profile: &Profile, origin: Point3, direction: Vec3, angle: f64) -> ExactBRep {
    let solid = revolve_profile_exact(profile, origin, direction, angle, Tolerance::METRE)
        .expect("an axis in the profile plane revolves");
    let health = geometric_audit(&solid, Tolerance::METRE);
    assert!(health.is_consistent(), "{:?}", health.defects());
    let topology = axiolid_topology::audit_brep(solid.topology());
    assert!(topology.is_closed_manifold(), "{topology:?}");
    solid
}

#[test]
fn a_tilted_axis_is_built_as_given_not_snapped_onto_y() {
    // x in [2, 3], y in [0, 1], about the axis through the origin tilted by
    // a milliradian: Pappus gives 2 pi (2.5 cos a - 0.5 sin a), where about +y
    // it would be 5 pi.
    let tilt: f64 = 1e-3;
    let axis = Vec3::new(tilt.sin(), tilt.cos(), 0.0);
    let section = polygon(&[(2.0, 0.0), (3.0, 0.0), (3.0, 1.0), (2.0, 1.0)]);
    let solid = revolve_about(&section, Point3::ZERO, axis, TAU);
    close(volume(&solid), 15.704813821837996);
    assert!((volume(&solid) - 5.0 * PI).abs() > 1e-3);
    // Every circle lies in a plane perpendicular to the axis, about it.
    let mut circles = 0;
    for curve in solid.curves3() {
        if let Curve3::Circle(circle) = curve {
            circles += 1;
            assert!(circle.frame.z.cross(axis).length() <= 1e-12);
            assert!(circle.frame.origin.cross(axis).length() <= 1e-12);
        }
    }
    assert!(circles >= 4);
    // The section stays where it was given, not mirrored across the axis:
    // each vertex is one of its corners to rounding.
    for vertex in solid.topology().vertices() {
        let corner = [(2.0, 0.0), (3.0, 0.0), (3.0, 1.0), (2.0, 1.0)]
            .iter()
            .any(|&(x, y)| (vertex.position - Vec3::new(x, y, 0.0)).length() <= 1e-12);
        assert!(corner, "{:?}", vertex.position);
    }
    // A quarter turn about the tilted axis.
    let quarter = revolve_about(&section, Point3::ZERO, axis, PI / 2.0);
    close(volume(&quarter), 3.926203455459499);
}

#[test]
fn an_axis_along_x_is_the_profiles_own_axis() {
    // y in [2, 3], x in [0, 1], about +x: 2 pi * 2.5 * 1.
    let section = polygon(&[(0.0, 2.0), (1.0, 2.0), (1.0, 3.0), (0.0, 3.0)]);
    let solid = revolve_about(&section, Point3::ZERO, Vec3::X, TAU);
    close(volume(&solid), 15.707963267948966);
    // The turn is about the x axis itself: every vertex is a corner of the
    // section or its mirror image across the axis (a full turn is the same
    // solid either way), in the profile plane.
    for vertex in solid.topology().vertices() {
        let p = vertex.position;
        assert!(
            (2.0 - 1e-12..=3.0 + 1e-12).contains(&p.y.abs()) && p.z.abs() <= 1e-12,
            "{p:?}"
        );
    }
    // A quarter turn by the right-hand rule about +x carries +y towards +z.
    let quarter = revolve_about(&section, Point3::ZERO, Vec3::X, PI / 2.0);
    close(volume(&quarter), 3.9269908169872414);
    for vertex in quarter.topology().vertices() {
        assert!(vertex.position.z >= -1e-12, "{:?}", vertex.position);
    }
    // And a downward axis keeps +y's sense: -x turns the other way.
    let against = revolve_about(&section, Point3::ZERO, Vec3::NEG_X, PI / 2.0);
    for vertex in against.topology().vertices() {
        assert!(vertex.position.z <= 1e-12, "{:?}", vertex.position);
    }
}

#[test]
fn a_cone_about_the_x_axis_has_its_apex_on_it() {
    let cone = polygon(&[(0.0, 0.0), (2.0, 0.0), (0.0, 1.0)]);
    let solid = revolve_about(&cone, Point3::ZERO, Vec3::X, TAU);
    close(volume(&solid), 2.0943951023931953);
}

#[test]
fn a_cone_about_a_slanted_axis_zeroes_the_axis_and_keeps_the_disc_flat() {
    // The axis (3, 4) passes through (0, 0) and (3, 4), exactly; the edge
    // from (5, 2.5) to (3, 4) is exactly perpendicular to it. Mapped into the
    // axis's frame both facts must survive the rounding of the rotation.
    let cone = polygon(&[(0.0, 0.0), (5.0, 2.5), (3.0, 4.0)]);
    let solid = revolve_about(&cone, Point3::ZERO, Vec3::new(3.0, 4.0, 0.0), TAU);
    close(volume(&solid), 32.72492347489368);
}

#[test]
fn a_disc_exactly_perpendicular_to_a_slanted_axis_keeps_its_height() {
    // (4, -3) is exactly perpendicular to the axis (3, 4) from the origin,
    // yet rotated its height along the axis is -4.4e-16, not the origin's
    // zero: the edge is a disc only because the height of its end on the
    // axis is the other end's.
    let cone = polygon(&[(0.0, 0.0), (4.0, -3.0), (3.0, 4.0)]);
    let solid = revolve_about(&cone, Point3::ZERO, Vec3::new(3.0, 4.0, 0.0), TAU);
    // Radius 5, height 5.
    close(volume(&solid), 125.0 * PI / 3.0);
    assert_eq!(count(&solid, |s| matches!(s, Surface::Plane(_))), 1);
    assert_eq!(count(&solid, |s| matches!(s, Surface::Cone(_))), 1);
}

#[test]
fn a_vertex_the_rotation_cannot_tell_from_the_axis_is_refused() {
    // (2.9999999999999902, 3.999999999999987) is 4e-16 left of the axis
    // along (3, 4), exactly, yet the rotated coordinate rounds to zero: the
    // axis frame cannot keep its side, so it is refused, not read as on the
    // axis or across it.
    named(
        revolve_profile_exact(
            &polygon(&[
                (0.0, 0.0),
                (4.0, -3.0),
                (2.9999999999999902, 3.999999999999987),
            ]),
            Point3::ZERO,
            Vec3::new(3.0, 4.0, 0.0),
            TAU,
            Tolerance::METRE,
        ),
        "closer to the axis than the rounding of the axis frame",
    );
}

#[test]
fn an_edge_a_hair_off_perpendicular_to_a_slanted_axis_is_refused() {
    // (4, -3) is exactly perpendicular to the axis (3, 4) from the origin;
    // the next float pair (3.9999999999999827, -2.999999999999987) is off by
    // 4e-16, yet its height along the axis rounds to exactly zero: a disc by
    // the heights, a cone by the edge. Refused, not read as a disc.
    named(
        revolve_profile_exact(
            &polygon(&[
                (0.0, 0.0),
                (3.9999999999999827, -2.999999999999987),
                (3.0, 4.0),
            ]),
            Point3::ZERO,
            Vec3::new(3.0, 4.0, 0.0),
            TAU,
            Tolerance::METRE,
        ),
        "within the tolerance of perpendicular",
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
        // The radii are differences of numbers of the axis's magnitude, so
        // the rounding the centre gate reads is the axis's, not the ball's.
        (1e-3, 12.5),
        (0.3, 1000.0),
        (1e-3, 1.0),
        (0.1, 100.0),
        (0.01, 100.0),
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
fn a_small_ball_far_from_the_origin_is_read_at_the_axis_magnitude() {
    // Radii 0.1 about x = 1e6: each is 1e-10 off, the centre computed from
    // them likewise.
    let ring = ArcRing::new(vec![
        ArcVertex {
            point: Point2::new(1e6, -0.1),
            bulge: (PI / 8.0).tan(),
        },
        ArcVertex {
            point: Point2::new(1e6 + 0.1, 0.0),
            bulge: (PI / 8.0).tan(),
        },
        ArcVertex {
            point: Point2::new(1e6, 0.1),
            bulge: 0.0,
        },
    ]);
    let solid = revolve_arc_ring(&ring, Point3::new(1e6, 0.0, 0.0), Tolerance::METRE)
        .expect("a ball far from the origin");
    let want = 4.0 / 3.0 * PI * 1e-3;
    assert!(
        (volume(&solid) - want).abs() <= 1e-8 * want,
        "{}",
        volume(&solid)
    );
}

#[test]
fn a_ball_about_a_distant_x_axis_is_read_at_its_premapped_magnitude() {
    // Mapped into the axis's frame the coordinates are small, but they are
    // differences of numbers near 1000.
    let ball = contour_of(vec![
        arc_segment(Point2::new(1000.0, 0.0), 0.3, 0.0, PI),
        line_segment(Point2::new(999.7, 0.0), Point2::new(1000.3, 0.0)),
    ]);
    let solid = revolve_about(&ball, Point3::new(1000.0, 0.0, 0.0), Vec3::X, TAU);
    let want = 4.0 / 3.0 * PI * 0.027;
    assert!(
        (volume(&solid) - want).abs() <= 1e-9 * want,
        "{}",
        volume(&solid)
    );
}

#[test]
fn a_semicircle_between_two_points_of_the_axis_is_a_ball() {
    // A bulge of exactly one is a sphere's meridian, split at its midpoint
    // as a half turn of an arc is; any other bulge between them is a lemon.
    for bulge in [1.0, -1.0] {
        let ring = ArcRing::new(vec![
            ArcVertex {
                point: Point2::new(0.0, -1.0),
                bulge,
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
        let solid = revolve_arc_ring(&ring, Point3::ZERO, Tolerance::METRE)
            .expect("a semicircle is a sphere");
        let health = geometric_audit(&solid, Tolerance::METRE);
        assert!(health.is_consistent(), "{:?}", health.defects());
        close(volume(&solid), 4.0 / 3.0 * PI);
    }
}

#[test]
fn a_near_flat_cone_the_representation_cannot_carry_is_refused_at_zero_tolerance() {
    // A rise of 1e-12 over a radius of 1: `atan` then `tan` is off the slope
    // by 2e-4 of the radius, and the cone leaves its edge by 7e-5.
    for points in [
        vec![(0.0, 0.0), (1.0, 1e-12), (0.0, 1e-12)],
        vec![(1.0, 0.0), (2.0, 1e-12), (2.0, 1.0), (1.0, 1.0)],
    ] {
        named(
            revolve_profile_exact(
                &polygon(&points),
                Point3::ZERO,
                Vec3::Y,
                TAU,
                Tolerance::ZERO,
            ),
            "cannot carry its slope",
        );
    }
}

#[test]
fn an_ordinary_cone_and_frustum_are_exact_at_zero_tolerance() {
    for (points, want) in [
        (vec![(0.0, 0.0), (1.5, 0.0), (0.0, 2.0)], 1.5 * PI),
        (
            vec![(1.0, 0.0), (2.0, 1.0), (2.0, 2.0), (1.0, 2.0)],
            // Annulus of inner radius 1: the frustum 1..2 over [0, 1] and
            // the cylinder of radius 2 over [1, 2], less the bore.
            13.0 * PI / 3.0,
        ),
    ] {
        let solid = revolve_profile_exact(
            &polygon(&points),
            Point3::ZERO,
            Vec3::Y,
            TAU,
            Tolerance::ZERO,
        )
        .expect("exact at zero tolerance");
        let health = geometric_audit(&solid, Tolerance::METRE);
        assert!(health.is_consistent(), "{:?}", health.defects());
        close(volume(&solid), want);
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
