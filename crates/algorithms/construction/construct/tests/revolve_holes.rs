//! Full-turn revolutions of sections with holes and of circles (#111).
//!
//! A hole swept a full turn encloses a ring-shaped cavity, carried as a void
//! shell; a circle sweeps a torus. Volumes are Pappus, `2 pi R A`, computed
//! from the inputs; every solid must audit clean and measure signed.

use axiolid_brep::ExactBRep;
use axiolid_brep_audit::geometric_audit;
use axiolid_construct::revolve_exact::revolve_profile_exact;
use axiolid_core::{Frame2, Interval, Point2, Point3, Tolerance, Vec2, Vec3};
use axiolid_curve::{Circle2, Curve2, Line2};
use axiolid_profile::{
    CircleProfile, Contour, ContourProfile, EllipseProfile, Profile, ProfileSegment,
};
use axiolid_surface::Surface;

const PI: f64 = std::f64::consts::PI;
const TAU: f64 = std::f64::consts::TAU;

fn revolve(profile: &Profile, axis_x: f64) -> ExactBRep {
    let solid = revolve_profile_exact(
        profile,
        Point3::new(axis_x, 0.0, 0.0),
        Vec3::Y,
        TAU,
        Tolerance::METRE,
    )
    .expect("revolves");
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

fn tori(solid: &ExactBRep) -> usize {
    solid
        .surfaces()
        .iter()
        .filter(|surface| matches!(surface, Surface::Torus(_)))
        .count()
}

#[test]
fn a_circle_revolves_into_a_torus() {
    // Profile centred on its own origin, axis 4 away: R = 4.
    let (r, major) = (1.25, 4.0);
    let solid = revolve(
        &Profile::Circle(CircleProfile {
            radius: r,
            thickness: None,
        }),
        -major,
    );
    assert_eq!(tori(&solid), 4, "four quarter arcs, four torus faces");
    close(volume(&solid), TAU * major * PI * r * r);
}

#[test]
fn a_hollow_circle_revolves_into_a_torus_with_a_toroidal_cavity() {
    let (r, t, major) = (1.25, 0.5, 4.0);
    let solid = revolve(
        &Profile::Circle(CircleProfile {
            radius: r,
            thickness: Some(t),
        }),
        -major,
    );
    assert_eq!(solid.topology().solids()[0].voids.len(), 1);
    let inner = r - t;
    close(volume(&solid), TAU * major * PI * (r * r - inner * inner));
}

#[test]
fn a_contour_with_two_holes_keeps_both_cavities() {
    // A 4 x 2 rectangle at x in [3, 7] with a square hole and a round hole.
    let line = |from: Point2, to: Point2| ProfileSegment {
        curve: Curve2::Line(Line2 {
            origin: from,
            direction: to - from,
        }),
        domain: Interval::UNIT,
        same_sense: true,
    };
    let polygon = |points: &[Point2]| {
        Contour::new(
            (0..points.len())
                .map(|i| line(points[i], points[(i + 1) % points.len()]))
                .collect(),
        )
    };
    let square = |x0: f64, y0: f64, x1: f64, y1: f64| {
        polygon(&[
            Point2::new(x0, y0),
            Point2::new(x1, y0),
            Point2::new(x1, y1),
            Point2::new(x0, y1),
        ])
    };
    let quarter = std::f64::consts::FRAC_PI_2;
    let (cx, cy, rr) = (6.0, 0.0, 0.5);
    let round = Contour::new(
        (0..4)
            .map(|index| ProfileSegment {
                curve: Curve2::Circle(Circle2 {
                    frame: Frame2 {
                        origin: Point2::new(cx, cy),
                        x: Vec2::X,
                        y: Vec2::Y,
                    },
                    radius: rr,
                }),
                domain: Interval::new(quarter * index as f64, quarter * (index + 1) as f64),
                same_sense: true,
            })
            .collect(),
    );
    let profile = Profile::Contour(ContourProfile {
        outer: square(3.0, -1.0, 7.0, 1.0),
        holes: vec![square(3.5, -0.5, 4.5, 0.5), round],
    });
    let solid = revolve(&profile, 0.0);
    assert_eq!(solid.topology().solids()[0].voids.len(), 2);
    // Pappus per piece: rectangle at R = 5, square hole at R = 4, round
    // hole at R = 6.
    let expected = TAU * (5.0 * 8.0 - 4.0 * 1.0 - 6.0 * PI * rr * rr);
    close(volume(&solid), expected);
}

#[test]
fn an_ellipse_revolution_is_still_refused_by_name() {
    // An ellipse swept about an axis in its plane is not a torus; no
    // surface family here carries it exactly.
    let result = revolve_profile_exact(
        &Profile::Ellipse(EllipseProfile {
            semi_axis_x: 1.0,
            semi_axis_y: 0.5,
        }),
        Point3::new(-4.0, 0.0, 0.0),
        Vec3::Y,
        TAU,
        Tolerance::METRE,
    );
    assert!(
        format!("{result:?}").contains("ellipse exact revolution"),
        "got {result:?}"
    );
}
