//! Exact full-turn revolution of a profile touching the axis (#379, step 1).
//!
//! An edge lying on the axis sweeps no face and an edge ending on it exactly
//! perpendicularly sweeps a disc cap: a solid cylinder, a stepped shaft.
//! Volumes are closed forms from the inputs (`pi r^2 h`, a sum of those);
//! every solid must audit clean. A profile crossing the axis, touching it at
//! an isolated vertex, or ending on it at an oblique edge stays refused by
//! name (steps 2 and 3 of #379).

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
    // the apex of two cones (step 2).
    named(
        try_revolve(
            &polygon(&[(1.0, 0.0), (2.0, 1.0), (1.0, 2.0), (0.0, 1.0)]),
            0.0,
        ),
        "isolated vertex",
    );
}

#[test]
fn an_oblique_edge_ending_on_the_axis_is_refused_by_name() {
    // A cone: its slant meets the axis at the apex, a pole (step 3).
    named(
        try_revolve(&polygon(&[(0.0, 0.0), (1.0, 0.0), (0.0, 1.0)]), 0.0),
        "oblique edge",
    );
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
    // The top edge rises by 1e-12 towards the axis: not a disc, exactly.
    named(
        try_revolve(
            &polygon(&[(0.0, 0.0), (1.0, 0.0), (1.0, 2.0), (0.0, 2.0 + 1e-12)]),
            0.0,
        ),
        "oblique edge",
    );
}

#[test]
fn an_arc_ending_on_the_axis_is_refused_as_a_pole() {
    // A quarter arc of the circle about (1, 0) starting at (0, 0) at angle pi,
    // where its x is exactly 0, then a straight run back along the top and
    // the axis: the arc's start is on the axis (step 3).
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
        "arc or an oblique edge",
    );
}

#[test]
fn an_arc_across_the_top_of_a_cylinder_is_not_a_disc() {
    // The top edge runs from (1, 1) to the axis at (0, 1), at the same
    // height, yet it is an arc, not a perpendicular line (step 3).
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
        "arc or an oblique edge",
    );
}
