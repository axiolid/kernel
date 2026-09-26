//! `FaceDomain` on a periodic face with no pole (#125, #167).
//!
//! A cylinder wall covering every angle must contain a point whatever turn
//! its angle is stated in: inversion may return `-0.2` or `2 pi - 0.2`.
//! The whole-period shifts once had the wrong sign, so a negative angle was
//! tried only further away from the domain and read as outside -- a region
//! of the ring below classified outside the ring, and the boolean dropped it.

use axiolid_construct::revolve_exact::revolve_profile_exact;
use axiolid_core::{Point2, Point3, Tolerance, Vec3};
use axiolid_measure::FaceDomain;
use axiolid_profile::{Profile, RectangleProfile};
use axiolid_surface::Surface;

const TAU: f64 = std::f64::consts::TAU;

#[test]
fn a_wall_contains_its_points_at_every_turn_of_the_angle() {
    let ring = revolve_profile_exact(
        &Profile::Rectangle(RectangleProfile {
            x: 2.0,
            y: 3.0,
            thickness: None,
            outer_radius: None,
            inner_radius: None,
        }),
        Point3::new(-5.0, 0.0, 0.0),
        Vec3::Y,
        TAU,
        Tolerance::METRE,
    )
    .expect("a ring");
    let topology = ring.topology();
    let wall = (0..topology.faces().len())
        .find(|&i| {
            matches!(
                &ring.surfaces()[topology.faces()[i].surface.unwrap().index()],
                Surface::Cylinder(c) if (c.radius - 6.0).abs() < 1e-12
            )
        })
        .expect("the outer wall");
    let face = topology.face_id_at(wall).unwrap();
    let domain = FaceDomain::new(&ring, face, Tolerance::METRE)
        .unwrap()
        .expect("classifiable");
    for u in [-0.2, TAU - 0.2, 2.0 * TAU - 0.2, -TAU - 0.2, 0.3] {
        assert_eq!(
            domain.contains(Point2::new(u, 1.5)).unwrap(),
            Some(true),
            "angle {u}"
        );
        assert_eq!(
            domain.contains(Point2::new(u, 3.5)).unwrap(),
            Some(false),
            "angle {u} above"
        );
    }
}
