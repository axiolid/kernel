//! Solids whose curved faces wind round their surface with no seam edge
//! (a dome bounded by its rim alone, a can by its two rims), as files may
//! deliver them. The boolean gives them seams first; the results audit
//! clean and measure to closed forms.

mod common;

use axiolid_brep::ExactBRep;
use axiolid_brep_audit::geometric_audit;
use axiolid_brep_boolean::boolean;
use axiolid_construct::boolean_exact::{boolean_arc_prisms_exact, ArcPrism};
use axiolid_core::{BooleanOperator, Frame3, Point2, Point3, Tolerance, Vec3};
use axiolid_measure::exact_properties;
use axiolid_overlay::ArcRing;
use axiolid_surface::{Cylinder, Sphere, Surface};
use common::seamless;

const PI: f64 = std::f64::consts::PI;

fn tol() -> Tolerance {
    Tolerance::METRE
}

fn world(z: f64) -> Frame3 {
    Frame3 {
        origin: Point3::new(0.0, 0.0, z),
        x: Vec3::X,
        y: Vec3::Y,
        z: Vec3::Z,
    }
}

fn volume(brep: &ExactBRep) -> f64 {
    exact_properties(brep, tol())
        .expect("measurable")
        .signed_volume
}

fn run(a: &ExactBRep, b: &ExactBRep, op: BooleanOperator) -> ExactBRep {
    let result = boolean(a, b, op, tol()).unwrap_or_else(|e| panic!("{op:?}: {e}"));
    let health = geometric_audit(&result, tol());
    assert!(health.is_consistent(), "{op:?}: {:?}", health.defects());
    let topology = axiolid_topology::audit_brep(result.topology());
    assert!(topology.is_closed_manifold(), "{op:?}: {topology:?}");
    result
}

fn close(what: &str, got: f64, want: f64) {
    assert!(
        (got - want).abs() <= 1e-9 * want.abs().max(1.0),
        "{what}: expected {want}, got {got}"
    );
}

fn block(x: (f64, f64), y: (f64, f64), z: (f64, f64)) -> ExactBRep {
    let p = ArcPrism {
        section: ArcRing::from_points(&[
            Point2::new(x.0, y.0),
            Point2::new(x.1, y.0),
            Point2::new(x.1, y.1),
            Point2::new(x.0, y.1),
        ]),
        bottom: z.0,
        top: z.1,
    };
    boolean_arc_prisms_exact(&p, &p, BooleanOperator::Intersection, tol()).expect("a box")
}

fn check(a: &ExactBRep, b: &ExactBRep, intersection: f64) {
    let u = volume(&run(a, b, BooleanOperator::Union));
    let i = volume(&run(a, b, BooleanOperator::Intersection));
    let d = volume(&run(a, b, BooleanOperator::Difference));
    close("intersection", i, intersection);
    close("identity", u + i, volume(a) + volume(b));
    close("difference", d, volume(a) - i);
}

#[test]
fn a_seamless_dome_halved() {
    let r = 1.2;
    let dome = seamless(
        Surface::Sphere(Sphere {
            frame: world(0.0),
            radius: r,
        }),
        &[(0.0, r, 0.0)],
    );
    close("dome", volume(&dome), 2.0 / 3.0 * PI * r * r * r);
    let half = block((0.0, 3.0), (-3.0, 3.0), (-1.0, 3.0));
    check(&dome, &half, PI * r * r * r / 3.0);
}

#[test]
fn a_seamless_can_cut_by_a_slab() {
    let r = 0.8;
    let can = seamless(
        Surface::Cylinder(Cylinder {
            frame: world(0.0),
            radius: r,
        }),
        &[(0.0, r, 0.0), (2.0, r, 2.0)],
    );
    close("can", volume(&can), PI * r * r * 2.0);
    let slab = block((-0.3, 0.4), (-3.0, 3.0), (0.5, 5.0));
    // The slab's width across the can's circle, times its height 1.5.
    let chord = |x: f64| 0.5 * (x * (r * r - x * x).sqrt() + r * r * (x / r).asin());
    let area = 2.0 * (chord(0.4) - chord(-0.3));
    check(&can, &slab, area * 1.5);
}
